
import json
from dataset import CIFilesDownloadProcess

if __name__ == "__main__":

    config_file_path = "config.json" # path to the configuration file

    with open(config_file_path, 'r') as config_file:
        config = json.load(config_file)

    app_key_path = config["app_key_path"]
    start_date = config["start_date"]
    end_date = config["end_date"]
    output_dir = config["output_dir"]
    json_path = config["json_path"]
    ci_script_path = config["ci_script_path"]
    convert_to_ci = config.get("convert_to_ci", False)
    remove_temp = config.get("remove_temp", True)

    d_down = CIFilesDownloadProcess(
        app_key_path=app_key_path,
        start_date=start_date,
        end_date=end_date,
        output_dir=output_dir,
        aoi_json_path=json_path,
        ci_python_script_path=ci_script_path,
        convert_to_ci=convert_to_ci,
        remove_temp=remove_temp
    )
    d_down.run_pipeline()

