# libraries
import logging
import json
from pathlib import Path

# method to setup configuration from JSON file
def setup_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

def setup_logging(log_cfg, date_str):

    log_dir_path = Path(log_cfg["base_path"])
    log_dir_path.mkdir(parents=True, exist_ok=True)

    log_file = log_dir_path / log_cfg["filename_pattern"].format(date=date_str)

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()  # Avoid duplicates if called multiple times

    # Formatter
    formatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console (prompt) handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


