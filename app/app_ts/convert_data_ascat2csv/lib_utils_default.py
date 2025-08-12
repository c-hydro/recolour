# ----------------------------------------------------------------------------------------------------------------------
# libraries
import argparse
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to parse command-line arguments
def parse_args() -> argparse.Namespace:

    p = argparse.ArgumentParser(description="ASCAT daily index builder (JSON-driven).")
    p.add_argument("-settings_file", required=True, help="Path to JSON config.")
    p.add_argument("-time-now", help="YYYY-MM-DD; used only if config lacks time_start/time_end.")
    p.add_argument("-require-exist", action="store_true",
                   help="Emit rows only when TIFF, CSV, and registry exist.")
    p.add_argument("-stdout", action="store_true",
                   help="Write a single combined CSV to stdout, ignoring 'out' pattern.")
    return p.parse_args()
# ----------------------------------------------------------------------------------------------------------------------
