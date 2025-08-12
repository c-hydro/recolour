
import os
import sys
import json
import rasterio
import shutil
from datetime import datetime, timedelta
from eccodes import *
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import box

## env information
## Step 1 : Create environment.yml
## Step 2: Create the Conda Environment
# conda env create -f environment.yml
## Step 3: Activate It
# source /root/envs/conda_default//miniconda/bin/activate /root/envs/conda_default//envs/recolour_default
## Step 4: Verify Dependencies
#python -c "import eccodes, geopandas, shapely; print('✅ Environment ready!')"

# Extra libraries
# conda install -c conda-forge eccodes python-eccodes
# conda install -c conda-forge rasterio

# === Load JSON config ===
def load_config(json_path):
    try:
        with open(json_path, 'r') as f:
            config = json.load(f)
        return config
    except Exception as e:
        print(f"❌ Failed to load config: {e}")
        sys.exit(1)

# === BUFR domain filtering ===
def file_intersects_domain(bufr_path, domain_union, lat_key, lon_key):
    try:
        with open(bufr_path, 'rb') as f:
            while True:
                try:
                    msg = codes_bufr_new_from_file(f)
                    if msg is None:
                        break
                    codes_set(msg, "unpack", 1)

                    lats = codes_get_array(msg, lat_key)
                    lons = codes_get_array(msg, lon_key)

                    for lat, lon in zip(lats, lons):
                        if -90 <= lat <= 90 and -180 <= lon <= 180:
                            if domain_union.intersects(Point(lon, lat)):
                                codes_release(msg)
                                return True

                    codes_release(msg)
                except EOFError:
                    break
                except Exception:
                    break
    except Exception:
        pass
    return False

def load_domain_geometry(path):
    if path.lower().endswith(".shp"):
        domain = gpd.read_file(path)
        return domain.geometry.unary_union
    elif path.lower().endswith((".tif", ".tiff")):
        with rasterio.open(path) as src:
            bounds = src.bounds
            return box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    else:
        print(f"❌ Unsupported domain file format: {path}")
        sys.exit(1)

# === Logging function ===
def log_message(message, log_file):
    print(message)
    with open(log_file, 'a') as log:
        log.write(message + '\n')

# === Main processing ===
def main(start_time_str, end_time_str, config_path):
    try:
        start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M")
        end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M")
    except ValueError:
        print("❌ Invalid time format. Use: 'YYYY-MM-DD HH:MM'")
        return

    config = load_config(config_path)
    input_root = config["INPUT_ROOT"]
    output_root = config["OUTPUT_ROOT"]
    domain_path = config["DOMAIN_PATH"]
    lat_key = config.get("LAT_KEY", "latitude")
    lon_key = config.get("LON_KEY", "longitude")
    log_file = config.get("LOG_FILE", "copied_files.log")

    # Clear log file
    with open(log_file, 'w') as log:
        log.write(f"# BUFR Copy Log: {start_time_str} to {end_time_str}\n")

    # Load shapefile
    log_message(f"🗺️  Loading domain from: {domain_path}", log_file)
    domain_union = load_domain_geometry(domain_path)

    log_message(f"📆 Filtering from {start_time} to {end_time}", log_file)
    log_message(f"📁 Input: {input_root}", log_file)
    log_message(f"📤 Output: {output_root}", log_file)

    current_time = start_time
    copied_count = 0

    while current_time < end_time:
        subdir = current_time.strftime("%Y/%m/%d/%H")
        input_folder = os.path.join(input_root, subdir)
        output_folder = os.path.join(output_root, subdir)

        if os.path.exists(input_folder):
            for fname in os.listdir(input_folder):
                if fname.endswith(".buf"):
                    full_input_path = os.path.join(input_folder, fname)
                    full_output_path = os.path.join(output_folder, fname)

                    if file_intersects_domain(full_input_path, domain_union, lat_key, lon_key):
                        os.makedirs(output_folder, exist_ok=True)
                        shutil.copy2(full_input_path, full_output_path)
                        copied_count += 1
                        log_message(f"✅ Copied: {full_input_path} → {full_output_path}", log_file)

        current_time += timedelta(hours=1)

    log_message(f"\n🎉 Done. Total files copied: {copied_count}", log_file)

# === CLI Entry ===
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage:")
        print("  python filter_bufr_by_time_and_region.py 'START_TIME' 'END_TIME' config.json")
        print("Format:")
        print("  START_TIME and END_TIME: 'YYYY-MM-DD HH:MM'")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])


