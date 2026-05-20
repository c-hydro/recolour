# Copyright (c) 2025, TU Wien
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of TU Wien, Department of Geodesy and Geoinformation
#      nor the names of its contributors may be used to endorse or promote
#      products derived from this software without specific prior written
#      permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL TU WIEN DEPARTMENT OF GEODESY AND
# GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import logging
import re

import numpy as np

from ascat.cell import RaggedArrayTs
from ascat.cell import OrthoMultiTimeseriesCell

from lib_utils_obj import split_path, create_sf_pattern
from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)

class BaseCellProduct:
    fn_format = "{:04d}.nc"

    @classmethod
    def preprocessor(cls, ds):
        return ds


class RaggedArrayCellProduct(BaseCellProduct):
    file_class = RaggedArrayTs
    sample_dim = "obs"
    instance_dim = "locations"

    @classmethod
    def preprocessor(cls, ds):
        if "row_size" in ds.variables:
            ds["row_size"].attrs["sample_dimension"] = cls.sample_dim
        if "locationIndex" in ds.variables:
            ds["locationIndex"].attrs["instance_dimension"] = cls.instance_dim
        if "location_id" in ds.variables:
            ds["location_id"].attrs["cf_role"] = "timeseries_id"
        if ds.attrs.get("featureType") is None:
            ds = ds.assign_attrs({"featureType": "timeSeries"})
        if ds.attrs.get("grid_mapping_name") is None:
            ds.attrs["grid_mapping_name"] = cls.grid_name
        return ds


class ErsCell(RaggedArrayCellProduct):

    @classmethod
    def preprocessor(cls, ds):
        if "obs" in ds.dims:
            chunk_dim = "obs"
        else:
            chunk_dim = "time"
        ds = super().preprocessor(ds).chunk({chunk_dim: -1})
        for var in ds.variables:
            if ds[var].dtype == np.float32:
                ds[var] = ds[var].where(ds[var] > -2147483600)
            if var == "alt":
                ds[var] = ds[var].where(ds[var] < 999999)

            parts = var.split("_")
            if parts[0] in ["fore", "mid", "aft"]:
                if parts[0] == "fore":
                    parts[0] = "for"
                ds = ds.rename({var: "_".join(parts[1:] + [parts[0]])})
        return ds


class ErsHCell(ErsCell):
    grid_name = "fibgrid_12.5"


class ErsNCell(ErsCell):
    grid_name = "fibgrid_25"


class AscatH129Cell(RaggedArrayCellProduct):
    grid_name = "fibgrid_6.25"


class AscatH130Cell(RaggedArrayCellProduct):
    grid_name = "fibgrid_6.25"


class AscatH122Cell(RaggedArrayCellProduct):
    grid_name = "fibgrid_6.25"


class AscatH121Cell(RaggedArrayCellProduct):
    grid_name = "fibgrid_12.5"


class AscatH139Cell(RaggedArrayCellProduct):
    grid_name = "fibgrid_12.5"


class AscatH29Cell(RaggedArrayCellProduct):
    grid_name = "fibgrid_12.5"


class AscatSIG0Cell6250m(RaggedArrayCellProduct):
    grid_name = "fibgrid_6.25"


class AscatSIG0Cell12500m(RaggedArrayCellProduct):
    grid_name = "fibgrid_12.5"


class OrthoMultiArrayCellProduct(BaseCellProduct):
    file_class = OrthoMultiTimeseriesCell
    sample_dim = "obs"
    instance_dim = "locations"

    @classmethod
    def preprocessor(cls, ds):
        if "location_id" in ds.variables:
            ds["location_id"].attrs["cf_role"] = "timeseries_id"
        if ds.attrs.get("featureType") is None:
            ds = ds.assign_attrs({"featureType": "timeSeries"})
        return ds


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


class AscatH129Swath(AscatSwathProduct):
    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-6.25km-H129_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {
        "satellite_folder": "metop_[abc]",
        "year_folder": "{year}",
        "month_folder": "{month}"
    }
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_6.25"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp, sat="[abc]"):
        sat = sat.lower()
        return {
            "satellite_folder": {
                "satellite": f"metop_{sat}"
            },
            "year_folder": {
                "year": f"{timestamp.year}"
            },
            "month_folder": {
                "month": f"{timestamp.month}".zfill(2)
            },
        }


class AscatH121Swath(AscatSwathProduct):
    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-12.5km-H121_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {
        "satellite_folder": "metop_[abc]",
        "year_folder": "{year}",
        "month_folder": "{month}"
    }
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_12.5"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp, sat="[abc]"):
        sat = sat.lower()
        return {
            "satellite_folder": {
                "satellite": f"metop_{sat}"
            },
            "year_folder": {
                "year": f"{timestamp.year}"
            },
            "month_folder": {
                "month": f"{timestamp.month}".zfill(2)
            },
        }


class AscatH122Swath(AscatSwathProduct):
    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-6.25km-H122_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {},
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_6.25"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp, sat="[abc]"):
        sat = sat.lower()
        return {
            "satellite_folder": {
                "satellite": f"metop_{sat}"
            },
        }

class AscatH122SwathNRT(AscatSwathProduct):

    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-6.25km-H122_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"

    sf_pattern_in = {
        "year_folder": "{year}",
        "month_folder": "{month}",
        "day_folder": "{day}",
        "hour_folder": "{hour}"
    }

    date_field_fmt_in = "%Y%m%d%H%M%S"
    date_field_fmt_out = "%Y%m%d_%H%M"
    grid_name = "fibgrid_6.25"
    cell_fn_format = "{:04d}"

    def __init__(self, path_in, filename_in, path_out, filename_out='{cell}.nc',
                 date_fmt_field_in='%Y%m%d%H%M%S', date_field_fmt_out='%Y%m%d',
                 cell_fn_format=cell_fn_format, grid_name=grid_name):

        self.path_in = path_in
        self.fn_pattern_in = filename_in

        self.path_out = path_out
        #filename_out = filename_out.replace('{cell}', cell_fn_format)
        self.fn_pattern_out = filename_out

        self.root_path_in, self.dynamic_path_in = split_path(self.path_in)
        self.root_path_out, self.dynamic_path_out = split_path(self.path_out)

        self.date_fmt_field_in = date_fmt_field_in
        self.date_field_fmt_out = date_field_fmt_out

        self.cell_fn_format = cell_fn_format
        self.grid_name = grid_name

        self.sf_pattern_in = create_sf_pattern(self.dynamic_path_in)
        self.sf_pattern_out = create_sf_pattern(self.dynamic_path_out)

    @staticmethod
    def fn_read_fmt(timestamp, date_fmt='%Y%m%d*', sat="[ABC]"):
        sat = sat.upper()

        return {
            "date": timestamp.strftime(date_fmt),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp):

        timestamp_parts = {}
        if hasattr(timestamp, "year"):
            timestamp_parts["year_folder"] = {
                "year": f"{timestamp.year}"
            }
        if hasattr(timestamp, "month"):
            timestamp_parts["month_folder"] = {
                "month": f"{timestamp.month:02d}"
            }
        if hasattr(timestamp, "day"):
            timestamp_parts["day_folder"] = {
                "day": f"{timestamp.day:02d}"
            }
        if hasattr(timestamp, "hour"):
            timestamp_parts["hour_folder"] = {
                "hour": f"{timestamp.hour:02d}"
            }
        if hasattr(timestamp, "minute"):
            timestamp_parts["minute_folder"] = {
                "minute": f"{timestamp.minute:02d}"
            }

        return timestamp_parts

    @staticmethod
    def fn_write_fmt(timestamp, date_fmt='%Y%m%d_%H00'):
        return {
            "date": timestamp.strftime(date_fmt),
            "cell": '{:4d}'
        }

    @staticmethod
    def sf_write_fmt(timestamp):

        timestamp_parts = {}
        if hasattr(timestamp, "year"):
            timestamp_parts["year_folder"] = {
                "year": f"{timestamp.year}"
            }
        if hasattr(timestamp, "month"):
            timestamp_parts["month_folder"] = {
                "month": f"{timestamp.month:02d}"
            }
        if hasattr(timestamp, "day"):
            timestamp_parts["day_folder"] = {
                "day": f"{timestamp.day:02d}"
            }
        if hasattr(timestamp, "hour"):
            timestamp_parts["hour_folder"] = {
                "hour": f"{timestamp.hour:02d}"
            }
        if hasattr(timestamp, "minute"):
            timestamp_parts["minute_folder"] = {
                "minute": f"{timestamp.minute:02d}"
            }

        return timestamp_parts

class AscatH122CellNRT(AscatSwathProduct):
    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-6.25km-H122_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {
        "year_folder": "{year}",
        "month_folder": "{month}",
        "day_folder": "{day}",
        "hour_folder": "{hour}"
    },
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_6.25"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def sf_read_fmt(timestamp, sat="[abc]"):
        sat = sat.lower()
        return {
            "year_folder": {
                "year": f"{timestamp.year}"
            },
            "month_folder": {
                "month": f"{timestamp.month}".zfill(2)
            },
            "day_folder": {
                "day": f"{timestamp.day}".zfill(2)
            },
            "hour_folder": {
                "hour": f"{timestamp.hour}".zfill(2)
            },
        }

    @staticmethod
    def sf_read_fmt(timestamp, sat="[abc]"):
        sat = sat.lower()
        return {
            "satellite_folder": {
                "satellite": f"metop_{sat}"
            },
            "year_folder": {
                "year": f"{timestamp.year}"
            },
            "month_folder": {
                "month": f"{timestamp.month}".zfill(2)
            },
            "day_folder": {
                "day": f"{timestamp.day}".zfill(2)
            },
            "hour_folder": {
                "hour": f"{timestamp.hour}".zfill(2)
            },
        }

class AscatH29Swath(AscatSwathProduct):
    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-12.5km-H29_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {}
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_12.5"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp):
        return {}


class AscatH130Swath(AscatSwathProduct):
    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-6.25km-H130_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {}
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_6.25"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp):
        return {}


class AscatH139Swath(AscatSwathProduct):
    fn_pattern = "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{sat}-12.5km-H139_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {}
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_12.5"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp):
        return {}


class AscatSIG0Swath6250m(AscatSwathProduct):
    """
    Class for reading ASCAT sigma0 swath data and writing it to cells.
    """
    fn_pattern = "W_IT-HSAF-ROME,SAT,SIG0-ASCAT-METOP{sat}-6.25_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {"satellite_folder": "metop_[abc]", "year_folder": "{year}"}
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_6.25"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        """
        Format a timestamp to search as YYYYMMDD*, for use in a regex
        that will match all files covering a single given date.

        Parameters
        ----------
        timestamp: datetime.datetime
            Timestamp to format

        Returns
        -------
        dict
            Dictionary of formatted strings
        """
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp, sat="[abc]"):
        sat = sat.lower()
        return {
            "satellite_folder": {
                "satellite": f"metop_{sat}"
            },
            "year_folder": {
                "year": f"{timestamp.year}"
            },
        }


class AscatSIG0Swath12500m(AscatSwathProduct):
    """
    Class for reading and writing ASCAT sigma0 swath data.
    """
    fn_pattern = "W_IT-HSAF-ROME,SAT,SIG0-ASCAT-METOP{sat}-12.5_C_LIIB_{placeholder}_{placeholder1}_{date}____.nc"
    sf_pattern = {"satellite_folder": "metop_[abc]", "year_folder": "{year}"}
    date_field_fmt = "%Y%m%d%H%M%S"
    grid_name = "fibgrid_12.5"
    cell_fn_format = "{:04d}.nc"

    @staticmethod
    def fn_read_fmt(timestamp, sat="[ABC]"):
        """
        Format a timestamp to search as YYYYMMDD*, for use in a regex
        that will match all files covering a single given date.

        Parameters
        ----------
        timestamp: datetime.datetime
            Timestamp to format

        Returns
        -------
        dict
            Dictionary of formatted strings
        """
        sat = sat.upper()
        return {
            "date": timestamp.strftime("%Y%m%d*"),
            "sat": sat,
            "placeholder": "*",
            "placeholder1": "*"
        }

    @staticmethod
    def sf_read_fmt(timestamp, sat="[abc]"):
        sat = sat.lower()
        return {
            "satellite_folder": {
                "satellite": f"metop_{sat}"
            },
            "year_folder": {
                "year": f"{timestamp.year}"
            },
        }


cell_io_catalog = {
    "H129": AscatH129Cell,
    "H130": AscatH130Cell,
    "H122": AscatH122Cell,
    "H122_NRT": AscatH122CellNRT,
    "H121": AscatH121Cell,
    "H139": AscatH139Cell,
    "H29": AscatH29Cell,
    "SIG0_6.25": AscatSIG0Cell6250m,
    "SIG0_12.5": AscatSIG0Cell12500m,
    "ERSH": ErsHCell,
    "ERSN": ErsNCell,
}

swath_io_catalog = {
    "H129": AscatH129Swath,
    "H130": AscatH130Swath,
    "H122": AscatH122Swath,
    "H122_NRT": AscatH122SwathNRT,
    "H121": AscatH121Swath,
    "H139": AscatH139Swath,
    "H29": AscatH29Swath,
    "SIG0_6.25": AscatSIG0Swath6250m,
    "SIG0_12.5": AscatSIG0Swath12500m,
}

swath_fname_regex_lookup = {
    "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP[ABC]-6.25km-H129_C_LIIB_.*_.*_.*____.nc":
        "H129",
    "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP[ABC]-6.25km-H130_C_LIIB_.*_.*_.*____.nc":
        "H130",
    "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP[ABC]-6.25km-H122_C_LIIB_.*_.*_.*____.nc":
        "H122",
    "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP[ABC]-12.5km-H121_C_LIIB_.*_.*_.*____.nc":
        "H121",
    "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP[ABC]-12.5km-H139_C_LIIB_.*_.*_.*____.nc":
        "H139",
    "W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP[ABC]-12.5km-H29_C_LIIB_.*_.*_.*____.nc":
        "H29",
    "W_IT-HSAF-ROME,SAT,SIG0-ASCAT-METOP[ABC]-6.25_C_LIIB_.*_.*_.*____.nc":
        "SIG0_6.25",
    "W_IT-HSAF-ROME,SAT,SIG0-ASCAT-METOP[ABC]-12.5_C_LIIB_.*_.*_.*____.nc":
        "SIG0_12.5",
}


def get_swath_product_id(filename):
    for pattern, swath_product_id in swath_fname_regex_lookup.items():
        if re.match(pattern, filename):
            return swath_product_id
    raise ValueError(
        f"Filename {filename} does not match any known swath product ID pattern."
    )
