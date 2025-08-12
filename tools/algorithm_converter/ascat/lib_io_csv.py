# libraries
import csv
import pandas as pd

# method to write a list of dictionaries to a CSV file
def write_csv(file_name, data, time_column='time', filter_column="ssm",
              filter_min=0, filter_max=100, filter_digits=5,
              delimiter=';'):

    if not data:
        raise ValueError("Data list is empty. Cannot create DataFrame.")

    df = pd.DataFrame(data)

    if time_column not in df.columns:
        raise KeyError(f"'{time_column}' column not found in data.")

    if filter_column not in df.columns:
        raise KeyError(f"'{filter_column}' column not found in data.")

    # Convert to datetime and drop invalid ones
    df[time_column] = pd.to_datetime(df[time_column], errors='coerce')
    df = df.dropna(subset=[time_column])

    # Filter by range
    df = df[(df[filter_column] >= filter_min) & (df[filter_column] <= filter_max)]

    # Ensure filter column is float and round to 3 decimals
    df[filter_column] = pd.to_numeric(df[filter_column], errors='coerce').round(filter_digits)

    # Set as index and sort
    df = df.set_index(time_column).sort_index()

    # Save to CSV
    df.to_csv(file_name, sep=delimiter, encoding="utf-8")
