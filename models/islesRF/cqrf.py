from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import data.usgs_ob as usgs
import data.caflow_ob as caflow
import datetime as dt
import joblib
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
import os

"""
Module contain functions that can be used to train a suite of Random Forest models to estimate CQ relationships.
See examples/random_forest.py for a notebook version of this file

Peter Isles' paper on Random Forests for CQ Realtionships: https://doi.org/10.1016/j.watres.2023.120876

Nutrient Data is gathered from Lake Champlain Tributary Long-term Monitoring Project:
 	https://dec.vermont.gov/watershed/lakes-ponds/monitor/lake-champlain-long-term-monitoring-project
	 
Streamflow data is gathered from the USGS Daily Values service using data/usgs_ob.py module.

In the script portion of this moduel, we are building RFs top predict total phosphorus (TP, ug/L) and total nitrogen (TN, mg/L) for each of the following 5 USGS sites:
	- Missisquoi River
	- Mill River
	- Jewett Brook
	- Rock River
	- Pike River
"""

def to_metricQ(streamflow_dict, streamflow_colname='streamflow'):
	'''
	Converts an entire dictionary of streamflow Series (such as that returned by usgs_ob.get_data()) to 
	a dictionary of dataframes, where streamflow has been converted from cubic feet /s to cubic m / s.
	
	Args:
	-- streamflow_dict (dict) [req]: dictionary of streamflow series, in the nested format as returned by usgs_ob.get_data()
	-- streamflow_colname (str) [opt]: the name of the streamflow series (key used in variables param of usgs_ob.get_data())
	
	Returns:
	-- metric_q_fs (dict): dictionary with the same location-name keys as streamflow_dict, with values being DataFrames of converted streamflow
	'''
	metric_q_dfs = {}
	for location, q in streamflow_dict.items():
		# make a copy of the data and convert to cubic m / s
		metric_q = q[streamflow_colname].copy() * 0.0283168
		metric_q_df = pd.DataFrame(metric_q)
		metric_q_df.columns = [f'{streamflow_colname} (m³/s)']
		metric_q_dfs[location] = metric_q_df
	return metric_q_dfs

def add_features(data, streamflow_colname='streamflow (m³/s)'):
	'''
	Function to add features to a daily discharge dataframe for a USGS gauge. Features added are the same as those used by Isles, 2023 (https://doi.org/10.1016/j.watres.2023.120876)
	Exact specifications for how features were engineered can be found in the code contained in the supplmentary material of said paper.

	Args:
	-- data (pd.DataFrame or pd.Series) [req]: daily discharge dataframe/series for a given USGS gauge. Should have datetime index and single streamflow column
	-- streamflow_colname (str) [opt]: name of the streamflow column

	Returns:
	features_df: a new dataframe containing the features needed for RF training
	'''
	df_features = pd.DataFrame(data).copy()
	# Natural log of Q
	df_features['ln_q'] = np.log(df_features[streamflow_colname])
	# Mean Q over days t-7 through t-1; calculated with logQ
	df_features['ln_q_7d_rolling_ave'] = df_features['ln_q'].rolling(window=dt.timedelta(days=7), min_periods=7).mean().shift(1)
	# Mean Q over days t-30 through t-1; calculated with logQ
	df_features['ln_q_30d_rolling_ave'] = df_features['ln_q'].rolling(window=dt.timedelta(days=30), min_periods=30).mean().shift(1)
	# ∆Q1day; calculated with Q
	df_features['q_1d_diff'] = df_features[streamflow_colname].diff()
	# day of year
	df_features['dofY'] = df_features.index.day_of_year
	# water year
	df_features['w_year'] = df_features.index.year
	df_features.index.name = 'date'
	df_features.index = df_features.index.tz_localize(None)
	return df_features.drop(columns=streamflow_colname)

def preprocess_nutrients(df, standardize_units=True):
	'''
	Preprocessing for nutrient data (TP and TN), gathered from thee Lake Champlain Tributary Long-term Monitoring Project
	https://dec.vermont.gov/watershed/lakes-ponds/monitor/lake-champlain-long-term-monitoring-project

	Args:
	-- df (pd.DataFrame) [req]: the raw nutrient dataframe, directly read in from csv
	-- standardize_units (bool) [opt]: whether or not to convert units to AEM3D standards

	Returns:
	df_processed (pd.DataFrame): a new dataframe containing a datetime index and a single column of nutrient data
	'''
	df_processed = df.copy()
	df_processed = df_processed.set_index(pd.to_datetime(df_processed['VisitDate']))
	df_processed.index.name = 'date'
	nutrient =  df_processed['Test'].iloc[0]
	match nutrient:
		case 'Total Phosphorus' | 'Dissolved Phosphorus':
			if nutrient == 'Total Phosphorus': nut='TP'
			elif nutrient == 'Dissolved Phosphorus': nut='DP'
			if standardize_units:
				nut_col = f'{nut}_mg/L'
				# convert from ug/L to mg/L
				df_processed[nut_col] = df_processed['Result'] / 1000
			else:
				nut_col = f'{nut}_ug/L'
				df_processed[nut_col] = df_processed['Result']
		case 'Total Nitrogen':
			nut_col = 'TN_mg/L'
			# TN already in mg/L
			df_processed[nut_col] = df_processed['Result']
		case _: raise ValueError(f"Unknown nutrient detected: {df_processed['Test'].iloc[0]}")
	df_processed = df_processed.drop(columns=[c for c in df_processed.columns if c != nut_col])
	
	# some dataframes have duplicate dates (multiple samples on the same day)
	# to handle duplicates, let's calculate the mean for a duplicated day and use that value
	grouped = df_processed.groupby(df_processed.index).agg({nut_col: 'mean'})
	# Remove duplicate rows based on the index, keeping the first occurrence
	print(f"Averaging the following duplicated rows...")
	print(df_processed[df_processed.index.duplicated(keep='first')])
	df_processed = df_processed[~df_processed.index.duplicated(keep='first')]
	# Update the column values with the calculated mean for duplicated timestamps
	df_processed.update(grouped)
	
	return df_processed

def combine_data(nutrient_dict, q_dict):
	'''
	Create a new dictionary of combined nutrient data and discharge features for a given location

	Args:
	-- nutrient_dict (dict) [req]: dictionary of nutrient data. Keys should be same location names as q_dict, values should be preprocessed nutrient data
	-- q_dict (dict) [req]: dictionary of discharge features dataframes, keyed by location.

	Returns:
	combined_dict (dict): a new dictionary of combined data, keys are the same as those of q_dict
	'''
	combined_dict = {}
	for loc, q_df in q_dict.items():
		nutrient_df = nutrient_dict[loc]
		combined_df = q_df.join(nutrient_df, how='inner', sort=True, validate='one_to_one')
		combined_dict[loc] = combined_df
	return combined_dict

def train_rf(data, nutrient_col, test_data_size=None, n_trees=500):
	'''
	Train a single random forest CQ model for a given tributary.

	Args:
	-- data (pd.DataFrame) [req]: The input data containing discharge feature columns and the target nutrient column.
	-- nutrient_col (str) [req]: The name of target nutrient column. Ex. 'TP_mg/L'
	-- test_data_size (float) [opt]: The proportion of the data to be used as the test set (between 0.0 and 1.0). Default is None (no split)
	-- n_trees (int) [opt]: The number of trees to be used in the Random Forest model. Default is 500 trees

	Returns:
	-- rf (RandomForestRegressor): The trained Random Forest model.
	-- x_test (pd.DataFrame): The test set features.
	-- y_test (pd.Series): The test set target values.
	'''
	# grab all features for X, nutrient column for Y
	x = data[[c for c in data.columns if c != nutrient_col]]
	y = data[nutrient_col]

	if test_data_size:
		# Split the data into training and testing sets
		x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_data_size, random_state=42)
	else:
		# if no train-test split, use all data
		x_train, y_train = x, y
		x_test, y_test = x, y

	# Create and train the Random Forest model
	rf = RandomForestRegressor(n_estimators=n_trees, random_state=42)
	rf.fit(x_train, y_train)

	return rf, x_test, y_test

def build_nutrient_models(model_data, nutrient, model_dir, test_data_size=None, save=False):
	'''
	Build and save a suite of random forest models for a given CQ relationship at multiple tributaries.

	Args:
	-- model_data (dict) [req]: A dictionary where keys are tributary names and values are dataframes containing complete model data (features and target) for that location.
	-- nutrient (str) [req]: The name of the nutrient column in the dataframes to be used as the target variable. Ex. 'TN_mg/L'
	-- model_dir (str) [req]: The directory where the trained models should be saved.
	-- test_data_size (float or None) [opt]: size of test data for split as a proportion between 0-1. Defaults to None (no train-test split).
	-- save (bool) [opt]: whether or not to save models.

	Returns:
	-- models (dict): A dictionary where keys are location names and values are dictionaries containing the trained model and the test data (`x_test`, `y_test`).
	'''
	models = {}
	for loc, data in model_data.items():
		# train model
		model, x_test, y_test = train_rf(data, nutrient, test_data_size)
		# saving the model with joblib package
		nutrient_name = nutrient.split('_')[0]
		model_path = os.path.join(model_dir, nutrient_name)
		model_fname = f'{loc}_{nutrient_name}.joblib'
		# make model save path if it doesn't exist
		if save:
			if not os.path.exists(model_path):
				os.makedirs(model_path)
				print(f'Model path created: {model_path}')
			# save the model
			joblib.dump(model, os.path.join(model_path, model_fname))
			print(f'Model saved to: {os.path.join(model_path, model_fname)}')
		models[loc] = {'model':model, 'x_test':x_test, 'y_test':y_test}
	return models


def feature_importances(model_dict):
	'''
	Calculate and print the feature importances for a given model.

	Args:
	-- model_dict (dict) [req]: A dictionary containing the trained model and test data, as returned by build_nutrient_models()[location]. Expected keys are:

	Returns:
	None
	'''
	# Calculate and print feature importances
	feature_importances = model_dict['model'].feature_importances_
	for feature, importance in zip(model_dict['x_test'].columns, feature_importances):
		print(f'Feature: {feature}, Importance: {importance}')

def plot_scatter_with_metrics(y_true, y_pred, nutrient, location, ax):
	# Calculate R-squared (R²) and Mean Squared Error (MSE)
	r2 = r2_score(y_true, y_pred)
	mse = mean_squared_error(y_true, y_pred)

	nutrient_name = nutrient.split('_')[0]
	unit = nutrient.split('_')[-1]
	
	# Create a scatter plot
	ax.scatter(y_true, y_pred, color='blue')

	# Add perfect predictions line
	ax.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], color='orange', linestyle='--', linewidth=2)

	# Add text for R-squared and MSE
	ax.text(0.30, 0.85, f'R² = {r2:.4f}', transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='right', color='r')
	ax.text(0.30, 0.75, f'MSE = {mse:.4f}', transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='right', color='r')

	# Add labels and title
	ax.set_xlabel(f'Actual {nutrient_name} ({unit})')
	ax.set_ylabel(f'Predicted {nutrient_name} ({unit})')
	ax.set_title(f'Random Forest {nutrient_name} Model for {location}')

def make_figure(tp_models, tn_models, n_row, n_col):
	# Create a grid of subplots
	fig, axs = plt.subplots(n_row, n_col, figsize=(14, 18))

	for i, (loc, model_dict) in enumerate(tp_models.items()):
		# Predict on the test set
		y_pred = model_dict['model'].predict(model_dict['x_test'])
		plot_scatter_with_metrics(model_dict['y_test'], y_pred, nutrient='TP_mg/L', location=loc, ax=axs[i,0])

	for i, (loc, model_dict) in enumerate(tn_models.items()):
		# Predict on the test set
		y_pred = model_dict['model'].predict(model_dict['x_test'])
		plot_scatter_with_metrics(model_dict['y_test'], y_pred, nutrient='TN_mg/L', location=loc, ax=axs[i,1])
	
	# Adjust the layout
	plt.subplots_adjust(hspace=0.5, wspace=0.3)  # Adjust these values as needed

	return fig

def plot_ts(complete_q_ts, complete_samples, reach, nutrient, save=False, plot_dir = "/users/n/b/nbeckage/ciroh/workspaces/notebooks/FEE/randForest/plots/"):
	match nutrient:
		case 'TP_mg/L':
			element = 'phosphorus'
		case 'TN_mg/L':
			element = 'nitrogen'

	fig = plt.figure(figsize=(12,6))
	ax = fig.add_subplot(1,1,1)

	ax.plot(complete_q_ts.index, complete_q_ts['preds'], color='orange', label='RandomForest', zorder=1)
	ax.scatter(complete_samples.index, complete_samples[nutrient], color='red', label="Samples")

	ax.grid(visible=True, which='both')

	ax.xaxis.set_major_locator(mdates.YearLocator())
	ax.xaxis.sbet_major_formatter(FuncFormatter(lambda x, pos: f'{mdates.num2date(x).year}' if mdates.num2date(x).year % 2 == 0 else ''))

	ax.set_xlabel('Date')
	ax.set_ylabel(f'Total {element.capitalize()} (mg/L)')
	ax.set_title(f'Total {element} RF predictions vs observationss for {reach.capitalize()}')
	ax.legend()
	if save:
		fig.savefig(os.path.join(plot_dir, f"{reach}_{element}_plot.png"))
	return fig

def plot_ts_interactive(complete_q_ts, complete_samples, reach, nutrient, save=False, plot_dir="/users/n/b/nbeckage/ciroh/workspaces/notebooks/FEE/randForest/plots/"):
	# TP, DP, or TN
	nut_abrv = nutrient.split("_")[0]
	
	match nut_abrv:
		case 'TP':
			full_nutrient = 'total phosphorus'
		case 'DP':
			full_nutrient = 'dissolved phosphorus'
		case 'TN':
			full_nutrient = 'total nitrogen'
		
	fig = go.Figure()

	fig.add_trace(go.Scatter(
		x=complete_q_ts.index,
		y=complete_q_ts['preds'],
		mode='lines',
		name='RandomForest',
		line=dict(color='orange'),
	))

	fig.add_trace(go.Scatter(
		x=complete_samples.index,
		y=complete_samples[nutrient],
		mode='markers',
		name='Samples',
		line=dict(color='red'),
	))

	# Customize the layout
	fig.update_layout(
		title=f'{full_nutrient.capitalize()} RF predictions vs observations for {reach.capitalize()}',
		xaxis_title='Date',
		yaxis_title=f'{full_nutrient} (mg/L)',
		xaxis=dict(
			tickmode='array',
			tickvals=pd.date_range(start=complete_samples.index.min(), end=complete_samples.index.max(), freq='YS'),  # Year start frequency for ticks
			ticktext=[str(year) if year % 2 == 0 else '' for year in range(complete_samples.index.min().year, complete_samples.index.max().year + 1)]  # Label every other year
		),
		showlegend=True
	)

	# Add gridlines
	fig.update_xaxes(showgrid=True)
	fig.update_yaxes(showgrid=True)

	# Show the plot
	fig.show()

	# Save the interactive figure as an HTML file
	if save:
		fig.write_html(os.path.join(plot_dir, f"{reach}_{nut_abrv}_interactive.html"))

def save_data(data_dict, save_dir, suffix):
	"""
	Dumps a dictionary of Pandas DataFrames or Series into a directory of CSV's.

	Args:
	-- data_dict (dict) [req]: Dictionary where keys are the names of the data (i.e. location, reach name, etc), and the values are the data, as Pandas DataFrames or Series.
	-- save_dir (str) [req]: Directory in which data will be saved.
	-- sufficx (str) [req]: File name suffic for CSV's, such as 'streamflow', 'nitro', etc.

	Returns:
	None
	"""
	for name, df in data_dict.items():
		save_path = os.path.join(save_dir, f"{name}_{suffix}.csv")
		print(f'Data saved at: {save_path}')
		df.to_csv(save_path)


##### MAIN SCRIPT #####

def main():
	# start and end date inferred from Peter's code
	start_dt = dt.datetime(1990, 10, 1, tzinfo=dt.timezone.utc)
	end_dt = dt.datetime(2023, 8, 31, tzinfo=dt.timezone.utc)

	# USGS gauges we want discharge for
	usgs_gauges = {'missisquoi':'04294000',
			'mill':'04292750',
			'jewett':'04292810'}
			# 'rock':'04294140',
			# 'pike':'04294300'}

	# Use Canadian data for rock and pike
	ca_guages = {'pike':'030424',
			  	 'rock':'030425'}

	# adding a day to end_dt because get_data()s are end date exclusive
	usgs_discharge = usgs.get_data(start_dt, end_dt+dt.timedelta(days=1), locations=usgs_gauges, variables={'discharge':'00060'}, service='dv')
	ca_discharge = caflow.get_data(start_dt, end_dt+dt.timedelta(days=1), locations=ca_guages, variables={'discharge':'Débit (m³/s)'}, service='dv')

	# convert USGS streamflow to metric
	metric_q_dfs = to_metricQ(usgs_discharge, streamflow_colname='discharge')

	# create a flattened single-layer dictionary of streamflow (necessary for combining dictionaries)
	ca_discharge_flat = {loc:loc_dict['discharge'] for loc, loc_dict in ca_discharge.items()}

	# now combine canadian and usgs data into one dictionary
	all_discharge = ca_discharge_flat | metric_q_dfs

	# add features to each df
	features_dfs = {location:add_features(discharge, streamflow_colname='discharge (m³/s)').dropna() for location, discharge in all_discharge.items()}

	# define directories where nutrient data are stored
	tp_dir = "/users/n/b/nbeckage/ciroh/workspaces/notebooks/FEE/randForest/nutrient_data/TP"
	tn_dir = "/users/n/b/nbeckage/ciroh/workspaces/notebooks/FEE/randForest/nutrient_data/TN"

	tp = {file.split('_')[0]:pd.read_csv(os.path.join(tp_dir, file)) for file in os.listdir(tp_dir)}
	tn = {file.split('_')[0]:pd.read_csv(os.path.join(tn_dir, file)) for file in os.listdir(tn_dir)}
	tp_processed = {location:preprocess_nutrients(tp_df) for location, tp_df in tp.items()}
	tn_processed = {location:preprocess_nutrients(tn_df) for location, tn_df in tn.items()}

	# combine nutrient and streamflow features
	full_tp = combine_data(tp_processed, features_dfs)
	full_tn = combine_data(tn_processed, features_dfs)

	# build and save models for nitrogren and phosphorus
	tp_models = build_nutrient_models(full_tp,
									  nutrient='TP_mg/L',
									  model_dir='/users/n/b/nbeckage/ciroh/workspaces/notebooks/FEE/randForest/models/',
									  test_data_size=None,
									  save=True)
	tn_models = build_nutrient_models(full_tn,
									  nutrient='TN_mg/L',
									  model_dir='/users/n/b/nbeckage/ciroh/workspaces/notebooks/FEE/randForest/models/',
									  test_data_size=None,
									  save=True)

if __name__ == "__main__":
	main()