import re
from pathlib import Path
import numpy as np
from ascat.cell import RaggedArrayTs
from ascat.cell import OrthoMultiTimeseriesCell

class SwathProduct:
    from ascat.swath import Swath
    file_class = Swath

class AscatSwathProduct(SwathProduct):
    grid_name = None

    @classmethod
    def preprocess_(cls, ds):
        ds["location_id"].attrs["cf_role"] = "timeseries_id"
        ds.attrs["global_attributes_flag"] = 1
        ds.attrs["featureType"] = "point"
        ds.attrs["grid_mapping_name"] = cls.grid_name
        if "spacecraft" in ds.attrs:
            # Assumption: the spacecraft attribute is something like "metop-a"
            sat_id = {"a": 3, "b": 4, "c": 5}
            sat = ds.attrs["spacecraft"][-1].lower()
            ds["sat_id"] = ("obs",
                            np.repeat(np.int8(sat_id[sat]), ds["location_id"].size))
            del ds.attrs["spacecraft"]
        return ds

    @staticmethod
    def postprocess_(ds):
        for key, item in {
                "latitude": "lat",
                "longitude": "lon",
                "altitude": "alt"
        }.items():
            if key in ds:
                ds = ds.rename({key: item})
        if "altitude" not in ds:
            ds["alt"] = ("locations",
                         np.full_like(ds["lat"], fill_value=np.nan))
        return ds

# helper to configure product
def configure_product(source):

    folder = source["folder"]
    filename = source["filename"]

    # folder:
    # /home/.../h122/{%Y/%m/%d/%H}/
    folder_parts = Path(folder).parts

    sf_pattern = {}
    sf_read_keys = {}

    for part in folder_parts:

        if "{%Y" in part:
            sf_pattern["year_folder"] = "{year}"
            sf_read_keys["year_folder"] = "year"

        elif "%m" in part:
            sf_pattern["month_folder"] = "{month}"
            sf_read_keys["month_folder"] = "month"

        elif "%d" in part:
            sf_pattern["day_folder"] = "{day}"
            sf_read_keys["day_folder"] = "day"

        elif "%H" in part:
            sf_pattern["hour_folder"] = "{hour}"
            sf_read_keys["hour_folder"] = "hour"

    # filename:
    # W_IT-HSAF-ROME,SAT,SSM-ASCAT-*_C_LIIB_*.nc
    fn_pattern = filename.replace(
        "SSM-ASCAT-*",
        "SSM-ASCAT-METOP{sat}-6.25km-H122"
    )

    fn_pattern = fn_pattern.replace(
        "_C_LIIB_*.nc",
        "_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    )

    return fn_pattern, sf_pattern, sf_read_keys

# helper to configure formats
def configure_format(sf_pattern):

    def sf_read_fmt(timestamp, sat="[abc]"):

        values = {
            "year": f"{timestamp.year}",
            "month": f"{timestamp.month}".zfill(2),
            "day": f"{timestamp.day}".zfill(2),
            "hour": f"{timestamp.hour}".zfill(2),
            "satellite": f"metop_{sat.lower()}",
        }

        sf_read = {}

        for folder_key, folder_pattern in sf_pattern.items():

            field_names = re.findall(r"{(.*?)}", folder_pattern)

            sf_read[folder_key] = {
                field_name: values[field_name]
                for field_name in field_names
            }

        return sf_read

    return staticmethod(sf_read_fmt)


# ----------------------------------------------------------------------------------------------------------------------
# class to initialize product
class AscatH122SwathNRT(AscatSwathProduct):

    def __init__(self, config, *args, **kwargs):

        self.config = config

        self.source = config["source"]
        self.product = config.get("product", {})

        self.fn_pattern = self._set_fn_pattern()
        self.sf_pattern = self._set_sf_pattern()

        self.date_field_fmt = self.product.get(
            "date_field_fmt", "%Y%m%d%H%M%S"
        )
        self.grid_name = self.product.get(
            "grid_name", "fibgrid_6.25"
        )
        self.cell_fn_format = self.product.get(
            "cell_fn_format", "{:04d}.nc"
        )

        super().__init__(*args, **kwargs)

    def _set_fn_pattern(self):

        return self.product.get(
            "filename_pattern",
            self.source["filename"]
        )

    def _set_sf_pattern(self):

        folder = self.source["folder"]

        sf_pattern = {}

        if "%Y" in folder:
            sf_pattern["year_folder"] = "{year}"

        if "%m" in folder:
            sf_pattern["month_folder"] = "{month}"

        if "%d" in folder:
            sf_pattern["day_folder"] = "{day}"

        if "%H" in folder:
            sf_pattern["hour_folder"] = "{hour}"

        return sf_pattern

    def sf_read_fmt(self, timestamp, sat="[abc]"):

        timestamp = pd.Timestamp(timestamp)

        values = {
            "year": f"{timestamp.year}",
            "month": f"{timestamp.month:02d}",
            "day": f"{timestamp.day:02d}",
            "hour": f"{timestamp.hour:02d}",
            "sat": sat.lower(),
            "satellite": f"metop_{sat.lower()}",
        }

        folder_values = {}

        for folder_name, folder_pattern in self.sf_pattern.items():

            keys = re.findall(r"{(.*?)}", folder_pattern)

            folder_values[folder_name] = {
                key: values[key]
                for key in keys
            }

        return folder_values
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# product helpers and catalogs
swath_io_catalog = {
    "H122_DR": AscatH122SwathNRT,
    "H122_NRT": AscatH122SwathNRT,
}

swath_fname_regex_lookup = {
    "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP[ABC]-6.25km-H122_C_LIIB_.*_.*_.*____.nc":
        "H122",
}

def get_swath_product_id(filename):
    for pattern, swath_product_id in swath_fname_regex_lookup.items():
        if re.match(pattern, filename):
            return swath_product_id
    raise ValueError(
        f"Filename {filename} does not match any known swath product ID pattern."
    )
# ----------------------------------------------------------------------------------------------------------------------
