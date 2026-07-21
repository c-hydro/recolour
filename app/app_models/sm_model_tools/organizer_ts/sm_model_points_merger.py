#!/usr/bin/env python3
"""
Merge daily point time-series CSV files over a backward day window.

For each configured variable, the application:

1. Defines a reference time, which may be "now" or a time in the past.
2. Searches one input CSV for each day from N days before through the
   reference day.
3. Merges all rows using the configured time column as the unique key.
4. When the same time step occurs in more than one file, keeps the row
   from the last file processed.
5. Writes one merged CSV per variable with rows sorted by time.

Daily files are processed from oldest to newest, so a newer daily file
has priority over an older daily file for duplicate time steps.

Example:
    python sm_model_points_merger.py \
        -settings_file sm_model_points_merger.json \
        -time "2026-07-10 23:00" \
        -days 10
"""

from __future__ import annotations

import argparse
import csv
import re
import json
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


def log_info(message: str) -> None:
    print(f"INFO  --> {message}")


def log_warning(message: str) -> None:
    print(f"WARN  --> {message}")


def parse_time(value: str, tz: timezone = timezone.utc) -> datetime:
    value = str(value).strip()

    if value.lower() == "now":
        return datetime.now(tz).replace(minute=0, second=0, microsecond=0)

    formats = (
        "%Y%m%d",
        "%Y%m%d%H",
        "%Y%m%d%H%M",
        "%Y-%m-%d",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M",
        "%Y-%m-%dT%H:%M:%S",
    )

    for fmt in formats:
        try:
            return datetime.strptime(value, fmt).replace(tzinfo=tz)
        except ValueError:
            pass

    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=tz)

    return parsed.astimezone(tz)


def build_daily_times(
    time_now: datetime,
    days_before: int,
    current_hour: Optional[int] = None,
    previous_hour: Optional[int] = None,
    include_reference_day: bool = True,
) -> List[datetime]:
    """
    Build input-file reference times.

    Historical files use previous_hour, including the current calendar day.

    The current execution file uses current_hour and is appended last, so it
    has priority when duplicate time steps are found.

    Example:
        time_now = 2026-06-18 15:00
        days_before = 3
        current_hour = 15
        previous_hour = 0

    Result:
        2026-06-15 00:00
        2026-06-16 00:00
        2026-06-17 00:00
        2026-06-18 00:00
        2026-06-18 15:00
    """

    if days_before < 0:
        raise ValueError("days_before must be >= 0")

    if current_hour is None:
        current_hour = time_now.hour

    if previous_hour is None:
        previous_hour = current_hour

    if not 0 <= current_hour <= 23:
        raise ValueError("current_hour must be between 0 and 23")

    if not 0 <= previous_hour <= 23:
        raise ValueError("previous_hour must be between 0 and 23")

    reference_day = time_now.replace(
        hour=0,
        minute=0,
        second=0,
        microsecond=0,
    )

    daily_times: List[datetime] = []

    # Historical files: oldest to newest.
    #
    # include_reference_day=True:
    #     include current calendar day at previous_hour.
    #
    # include_reference_day=False:
    #     stop at the day before the current calendar day.
    last_offset = 0 if include_reference_day else 1

    if last_offset <= days_before:
        for day_offset in range(days_before, last_offset - 1, -1):
            historical_time = (
                reference_day - timedelta(days=day_offset)
            ).replace(
                hour=previous_hour,
                minute=0,
                second=0,
                microsecond=0,
            )

            daily_times.append(historical_time)

    # Current execution file.
    current_time = reference_day.replace(
        hour=current_hour,
        minute=0,
        second=0,
        microsecond=0,
    )

    # Avoid adding the same file twice when both hours are equal and the
    # current reference day is already included.
    if current_time not in daily_times:
        daily_times.append(current_time)

    return daily_times


def format_path(
    template: str,
    time_step: datetime,
    time_now: datetime,
    var_name: str,
) -> str:
    return template.format(
        time=time_step,
        time_step=time_step,
        time_now=time_now,
        var_name=var_name,
        yyyy=time_step.strftime("%Y"),
        mm=time_step.strftime("%m"),
        dd=time_step.strftime("%d"),
        yyyymmdd=time_step.strftime("%Y%m%d"),
        yyyymmddhh=time_step.strftime("%Y%m%d%H"),
    )


def normalize_variables(cfg: Dict[str, Any]) -> List[str]:
    variables = cfg.get("variables")

    if isinstance(variables, list):
        normalized = [str(name).strip() for name in variables if str(name).strip()]
    elif isinstance(variables, dict):
        normalized = [
            str(name).strip()
            for name, var_cfg in variables.items()
            if not isinstance(var_cfg, dict) or var_cfg.get("enabled", True)
        ]
    elif variables is None:
        normalized = []
    else:
        raise TypeError("'variables' must be a list or dictionary")

    if not normalized:
        raise RuntimeError("No variables configured")

    return normalized


def read_csv_rows(
    file_path: str,
    delimiter: str,
    encoding: str,
    time_column: str,
) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(file_path, "r", newline="", encoding=encoding) as file_handle:
        reader = csv.DictReader(file_handle, delimiter=delimiter)

        if reader.fieldnames is None:
            raise RuntimeError(f"CSV file has no header: {file_path}")

        header = [str(name).strip() for name in reader.fieldnames]

        if time_column not in header:
            raise KeyError(
                f"Time column '{time_column}' not found in {file_path}. "
                f"Available columns: {header}"
            )

        rows: List[Dict[str, str]] = []
        for row_number, row in enumerate(reader, start=2):
            clean_row = {
                str(key).strip(): "" if value is None else str(value).strip()
                for key, value in row.items()
                if key is not None
            }

            time_value = clean_row.get(time_column, "")
            if not time_value:
                log_warning(
                    f"FILE {file_path} | ROW {row_number}: empty "
                    f"'{time_column}', row skipped"
                )
                continue

            rows.append(clean_row)

    # define natural keys
    def natural_key(value: str):
        return [
            int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", value)
        ]

    # Sort header keeping time_column first
    header = [time_column] + sorted(
        (column for column in header if column != time_column),
        key=natural_key,
    )

    # Reorder every row using the sorted header
    rows = [
        {column: row.get(column, "") for column in header}
        for row in rows
    ]

    return header, rows


def merge_headers(
    current_header: Optional[List[str]],
    new_header: Sequence[str],
    time_column: str,
    strict_headers: bool,
    file_path: str,
) -> List[str]:
    new_header = list(new_header)

    if current_header is None:
        return new_header

    if current_header == new_header:
        return current_header

    if strict_headers:
        raise RuntimeError(
            f"Header mismatch in {file_path}.\n"
            f"Expected: {current_header}\n"
            f"Found   : {new_header}"
        )

    merged = [time_column]

    for field_name in list(current_header) + new_header:
        if field_name != time_column and field_name not in merged:
            merged.append(field_name)

    log_warning(
        f"FILE {file_path}: header differs; using the union of all columns"
    )
    return merged


def parse_row_time(value: str, configured_format: Optional[str]) -> Tuple[int, Any]:
    value = str(value).strip()

    if configured_format:
        try:
            return 0, datetime.strptime(value, configured_format)
        except ValueError:
            pass

    formats = (
        "%Y%m%d%H",
        "%Y%m%d%H%M",
        "%Y%m%d",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d",
    )

    for fmt in formats:
        try:
            return 0, datetime.strptime(value, fmt)
        except ValueError:
            pass

    # Unknown formats are still deterministic and are sorted lexically
    # after recognized datetime formats.
    return 1, value


def write_merged_csv(
    file_path: str,
    header: Sequence[str],
    rows_by_time: Dict[str, Dict[str, str]],
    time_column: str,
    delimiter: str,
    encoding: str,
    time_format: Optional[str],
) -> None:
    Path(file_path).parent.mkdir(parents=True, exist_ok=True)

    sorted_time_keys = sorted(
        rows_by_time,
        key=lambda value: parse_row_time(value, time_format),
    )

    with open(file_path, "w", newline="", encoding=encoding) as file_handle:
        writer = csv.DictWriter(
            file_handle,
            fieldnames=list(header),
            delimiter=delimiter,
            extrasaction="ignore",
        )
        writer.writeheader()

        for time_key in sorted_time_keys:
            source_row = rows_by_time[time_key]
            output_row = {
                field_name: source_row.get(field_name, "")
                for field_name in header
            }
            writer.writerow(output_row)

def get_daily_time_template(
    rows: Sequence[Dict[str, str]],
    time_column: str,
    time_format: str,
) -> List[Tuple[int, int]]:
    """
    Extract the unique hour/minute combinations available in a daily file.

    Example:
        2026070100
        2026070101
        ...
        2026070123

    returns:
        [(0, 0), (1, 0), ..., (23, 0)]
    """
    daily_steps: List[Tuple[int, int]] = []

    for row in rows:
        time_value = row.get(time_column, "")

        if not time_value:
            continue

        try:
            time_step = datetime.strptime(time_value, time_format)
        except ValueError as exc:
            raise ValueError(
                f"Cannot parse time value '{time_value}' "
                f"using format '{time_format}'"
            ) from exc

        hour_minute = (time_step.hour, time_step.minute)

        if hour_minute not in daily_steps:
            daily_steps.append(hour_minute)

    daily_steps.sort()

    if not daily_steps:
        raise RuntimeError(
            "No valid time steps found in the reference daily file"
        )

    return daily_steps

from datetime import datetime
from typing import Dict, List, Sequence, Tuple

import pandas as pd


def create_missing_daily_rows(
    daily_time: datetime,
    daily_steps: Sequence[Tuple[int, int]],
    header: Sequence[str],
    time_column: str,
    time_format: str,
    missing_value: float,
    direction: str = "backward",
) -> List[Dict[str, str]]:
    """
    Create missing rows for the date defined by daily_time.

    daily_time is used only to select the year, month and day.
    Hours and minutes are taken exclusively from daily_steps.

    direction="forward":
        daily_steps are ordered from first to last.

    direction="backward":
        daily_steps are ordered from last to first.
    """

    direction = direction.strip().lower()

    if direction not in {"forward", "backward"}:
        raise ValueError(
            "direction must be either 'forward' or 'backward'"
        )

    if not daily_steps:
        raise ValueError("daily_steps cannot be empty")

    # Validate, remove duplicates and sort by hour/minute
    daily_steps_sorted = sorted(
        {
            (int(hour), int(minute))
            for hour, minute in daily_steps
        }
    )

    for hour, minute in daily_steps_sorted:
        if not 0 <= hour <= 23:
            raise ValueError(
                f"Invalid hour in daily_steps: {hour}"
            )

        if not 0 <= minute <= 59:
            raise ValueError(
                f"Invalid minute in daily_steps: {minute}"
            )

    if direction == "backward":
        selected_steps = reversed(daily_steps_sorted)
    else:
        selected_steps = iter(daily_steps_sorted)

    missing_rows: List[Dict[str, str]] = []

    for hour, minute in selected_steps:
        time_step = daily_time.replace(
            hour=hour,
            minute=minute,
            second=0,
            microsecond=0,
        )

        row = {
            field_name: str(missing_value)
            for field_name in header
        }

        row[time_column] = time_step.strftime(time_format)
        missing_rows.append(row)

    return missing_rows

def filter_merged_data_by_time(
    rows_by_time: Dict[str, Dict[str, str]],
    source_by_time: Dict[str, str],
    time_reference: datetime,
    time_format: str,
    keep_until_reference: bool = True,
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str]]:
    """
    Keep only merged rows with time <= time_reference.

    The same filtering is applied to rows_by_time and source_by_time.
    """

    if not keep_until_reference:
        return rows_by_time, source_by_time

    filtered_rows: Dict[str, Dict[str, str]] = {}
    filtered_sources: Dict[str, str] = {}

    removed_count = 0

    for time_key, row in rows_by_time.items():

        try:
            row_time = datetime.strptime(time_key, time_format)
        except ValueError:
            log_warning(
                f"TIME {time_key}: cannot parse using format "
                f"'{time_format}', row removed"
            )
            removed_count += 1
            continue

        # Align timezone awareness if needed
        if time_reference.tzinfo is not None:
            row_time = row_time.replace(tzinfo=time_reference.tzinfo)

        if row_time <= time_reference:
            filtered_rows[time_key] = row

            if time_key in source_by_time:
                filtered_sources[time_key] = source_by_time[time_key]
        else:
            removed_count += 1

    log_info(
        f"TIME FILTER <= {time_reference:%Y-%m-%d %H:%M}: "
        f"{len(filtered_rows)}/{len(rows_by_time)} rows kept, "
        f"{removed_count} removed"
    )

    return filtered_rows, filtered_sources

def sort_merged_data(
    rows_by_time: Dict[str, Dict[str, str]],
    source_by_time: Dict[str, str],
    time_format: str,
    sort_order: str = "ascending",
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str]]:
    """
    Sort rows_by_time and source_by_time using the same time-key order.

    sort_order:
        ascending
        descending
    """

    sort_order = str(sort_order).strip().lower()

    if sort_order not in {"ascending", "descending"}:
        raise ValueError(
            "sort_order must be either 'ascending' or 'descending'"
        )

    reverse = sort_order == "descending"

    def get_sort_key(time_key: str) -> datetime:
        try:
            return datetime.strptime(time_key, time_format)
        except ValueError as exc:
            raise ValueError(
                f"Cannot parse time key '{time_key}' "
                f"using format '{time_format}'"
            ) from exc

    sorted_time_keys = sorted(
        rows_by_time.keys(),
        key=get_sort_key,
        reverse=reverse,
    )

    sorted_rows = {
        time_key: rows_by_time[time_key]
        for time_key in sorted_time_keys
    }

    sorted_sources = {
        time_key: source_by_time[time_key]
        for time_key in sorted_time_keys
        if time_key in source_by_time
    }

    return sorted_rows, sorted_sources

def merge_variable(
    var_name: str,
    daily_times: Sequence[datetime],
    time_now: datetime,
    input_cfg: Dict[str, Any],
    output_cfg: Dict[str, Any],
) -> None:
    input_folder_template = input_cfg.get("folder_name", "")
    input_file_template = input_cfg["file_template"]
    input_delimiter = input_cfg.get("delimiter", ",")
    input_encoding = input_cfg.get("encoding", "utf-8-sig")
    missing_file_value = float(input_cfg.get("missing_file_value", -9998.0))
    missing_file_method = input_cfg.get("missing_file_method", "backward")

    output_folder_template = output_cfg.get("folder_name", "./output")
    output_file_template = output_cfg.get("file_template","obs_db_{var_name}_{time_now:%Y%m%d%H}.csv",)
    output_delimiter = output_cfg.get("delimiter", input_delimiter)
    output_encoding = output_cfg.get("encoding", "utf-8")
    output_sort_order = output_cfg.get("sort_order", "ascending").lower()
    output_until_reference = bool( output_cfg.get("keep_until_reference", True))

    time_column = input_cfg.get("time_column", "time")
    time_format = input_cfg.get("time_format")
    strict_headers = bool(input_cfg.get("strict_headers", True))
    fail_if_no_files = bool(input_cfg.get("fail_if_no_files", True))

    rows_by_time: Dict[str, Dict[str, str]] = {}
    source_by_time: Dict[str, str] = {}
    merged_header: Optional[List[str]] = None

    files_found = 0
    rows_read = 0
    duplicate_count = 0

    # get daily hour from daily times
    current_hour = daily_times[-1].hour

    log_info(f"VARIABLE {var_name}: start merge")

    daily_steps_template: Optional[List[Tuple[int, int]]] = None
    daily_header_template: Optional[List[str]] = None

    for file_time in daily_times:
        log_info(f"DAILY FILE TIME     : {file_time:%Y-%m-%d %H:%M}")

        folder_name = format_path(
            input_folder_template,
            file_time,
            time_now,
            var_name,
        )
        file_name = format_path(
            input_file_template,
            file_time,
            time_now,
            var_name,
        )
        file_path = os.path.join(folder_name, file_name)

        log_info(
            f"VARIABLE {var_name} | DAY {file_time:%Y-%m-%d} | "
            f"FILE {file_path}"
        )

        if os.path.exists(file_path):

            header, rows = read_csv_rows(
                file_path=file_path,
                delimiter=input_delimiter,
                encoding=input_encoding,
                time_column=time_column,
            )

            # Get the daily structure from the first available file.
            if daily_steps_template is None:
                daily_steps_template = get_daily_time_template(
                    rows=rows,
                    time_column=time_column,
                    time_format=time_format,
                )
                daily_header_template = list(header)

                log_info(
                    f"VARIABLE {var_name} | DAILY TEMPLATE: "
                    f"{len(daily_steps_template)} steps"
                )
                log_info(
                    f"VARIABLE {var_name} | AVAILABLE HOURS: "
                    + ", ".join(
                        f"{hour:02d}:{minute:02d}"
                        for hour, minute in daily_steps_template
                    )
                )

            else:
                current_daily_steps = get_daily_time_template(
                    rows=rows,
                    time_column=time_column,
                    time_format=time_format,
                )

                if current_daily_steps != daily_steps_template:
                    log_warning(
                        f"VARIABLE {var_name} | FILE {file_path}: "
                        f"daily time structure differs. "
                        f"Expected {len(daily_steps_template)} steps, "
                        f"found {len(current_daily_steps)}"
                    )

            merged_header = merge_headers(
                current_header=merged_header,
                new_header=header,
                time_column=time_column,
                strict_headers=strict_headers,
                file_path=file_path,
            )

            files_found += 1
            rows_read += len(rows)

        else:
            log_warning(
                f"VARIABLE {var_name} | FILE NOT FOUND: {file_path}"
            )

            if (
                    daily_steps_template is None
                    or daily_header_template is None
            ):
                log_warning(
                    f"VARIABLE {var_name} | Cannot fill missing day "
                    f"{daily_time:%Y-%m-%d}: daily template is not "
                    f"available yet"
                )
                continue

            header = list(daily_header_template)

            rows = create_missing_daily_rows(
                daily_time=daily_time,
                daily_steps=daily_steps_template,
                header=header,
                time_column=time_column,
                time_format=time_format,
                missing_value=missing_file_value,
                direction=missing_file_method
            )

            log_info(
                f"VARIABLE {var_name} | CREATED {len(rows)} "
                f"MISSING ROWS FOR {daily_time:%Y-%m-%d}"
            )

        for row in rows:
            time_key = row[time_column]

            if time_key in rows_by_time:
                duplicate_count += 1
                previous_file = source_by_time[time_key]

                log_info(
                    f"VARIABLE {var_name} | DUPLICATE {time_key}: "
                    f"replace row from {previous_file} "
                    f"with row from {file_path}"
                )

            rows_by_time[time_key] = row

            if os.path.exists(file_path):
                source_by_time[time_key] = file_path
            else:
                source_by_time[time_key] = "generated_missing_file"

    if files_found == 0:
        message = (
            f"VARIABLE {var_name}: no input files found in the requested period"
        )
        if fail_if_no_files:
            raise FileNotFoundError(message)
        log_warning(message)
        return

    if merged_header is None:
        raise RuntimeError(f"VARIABLE {var_name}: no valid CSV headers found")

    # filter data using time reference
    rows_by_time, source_by_time = filter_merged_data_by_time(
        rows_by_time=rows_by_time,
        source_by_time=source_by_time,
        time_reference=time_now,
        time_format=time_format,
        keep_until_reference=output_until_reference,
    )
    # sort data using method (ascending or descending)
    rows_by_time, source_by_time = sort_merged_data(
        rows_by_time=rows_by_time,
        source_by_time=source_by_time,
        time_format=time_format,
        sort_order=output_sort_order,
    )

    output_time = time_now.replace(
        hour=current_hour,
        minute=0,
        second=0,
        microsecond=0,
    )

    output_folder = format_path(
        output_folder_template,
        output_time,
        output_time,
        var_name,
    )
    output_file = format_path(
        output_file_template,
        output_time,
        output_time,
        var_name,
    )
    output_path = os.path.join(output_folder, output_file)

    write_merged_csv(
        file_path=output_path,
        header=merged_header,
        rows_by_time=rows_by_time,
        time_column=time_column,
        delimiter=output_delimiter,
        encoding=output_encoding,
        time_format=time_format,
    )

    log_info(f"VARIABLE {var_name}: files found       = {files_found}")
    log_info(f"VARIABLE {var_name}: rows read         = {rows_read}")
    log_info(f"VARIABLE {var_name}: duplicate rows    = {duplicate_count}")
    log_info(f"VARIABLE {var_name}: unique rows       = {len(rows_by_time)}")
    log_info(f"VARIABLE {var_name}: merged file       = {output_path}")

# main of application
def main() -> None:

    parser = argparse.ArgumentParser(description="Merge daily point time-series CSV files")
    parser.add_argument(
        "-settings_file",
        "--settings_file",
        required=True,
        help="JSON settings file",
    )
    parser.add_argument(
        "-time",
        "--time",
        dest="time_now",
        default=None,
        help='Reference time, for example "2026-07-10 23:00" or "now"',
    )
    parser.add_argument(
        "-days",
        "--days",
        dest="days_before",
        type=int,
        default=None,
        help="Number of days before the reference day",
    )

    parser.add_argument(
        "--hour-current-day",
        dest="hour_current_day",
        type=int,
        default=None,
        help="Input-file hour for the current execution day, from 0 to 23",
    )

    parser.add_argument(
        "--hour-previous-day",
        dest="hour_previous_day",
        type=int,
        default=None,
        help=(
            "Input-file hour used for the historical period, including "
            "the current calendar day, from 0 to 23"
        ),
    )


    args = parser.parse_args()

    cfg = json.loads(Path(args.settings_file).read_text(encoding="utf-8"))
    time_cfg = cfg.get("time", {})

    tz = timezone.utc
    configured_time = cfg.get(
        "time_now",
        time_cfg.get("time_now", "now"),
    )
    time_now = parse_time(
        args.time_now if args.time_now is not None else configured_time,
        tz=tz,
    )

    configured_days = time_cfg.get("days_before", 10)
    days_before = (
        args.days_before
        if args.days_before is not None
        else int(configured_days)
    )

    current_hour_cfg = time_cfg.get(
        "hour_current_day",
        time_cfg.get("daily_hour", time_now.hour),
    )

    previous_hour_cfg = time_cfg.get(
        "hour_previous_day",
        current_hour_cfg,
    )

    current_hour = (
        args.hour_current_day
        if args.hour_current_day is not None
        else int(current_hour_cfg)
    )

    previous_hour = (
        args.hour_previous_day
        if args.hour_previous_day is not None
        else int(previous_hour_cfg)
    )

    include_reference_day = bool(
        time_cfg.get("include_reference_day", True)
    )

    daily_times = build_daily_times(
        time_now=time_now,
        days_before=days_before,
        current_hour=current_hour,
        previous_hour=previous_hour,
        include_reference_day=include_reference_day,
    )

    log_info("POINTS MERGER - START")
    log_info(f"SETTINGS FILE       : {args.settings_file}")
    log_info(f"REFERENCE TIME      : {time_now:%Y-%m-%d %H:%M}")
    log_info(f"DAYS BEFORE         : {days_before}")
    log_info(f"CURRENT DAY HOUR    : {current_hour:02d}:00")
    log_info(f"PREVIOUS DAYS HOUR  : {previous_hour:02d}:00")
    log_info(f"INCLUDE REF DAY     : {include_reference_day}")
    log_info(f"NUMBER OF DAY FILES : {len(daily_times)}")

    if daily_times:
        log_info(
            f"FILE PERIOD         : "
            f"{daily_times[0]:%Y-%m-%d %H:%M} -> "
            f"{daily_times[-1]:%Y-%m-%d %H:%M}"
        )

    variables = normalize_variables(cfg)
    input_cfg = cfg.get("input", {})
    output_cfg = cfg.get("output", {})

    if "file_template" not in input_cfg:
        raise KeyError("Missing required configuration: input.file_template")

    for var_name in variables:
        merge_variable(
            var_name=var_name,
            daily_times=daily_times,
            time_now=time_now,
            input_cfg=input_cfg,
            output_cfg=output_cfg,
        )

    log_info("POINTS MERGER - END")


if __name__ == "__main__":
    main()
