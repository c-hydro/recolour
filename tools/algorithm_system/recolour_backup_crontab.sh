#!/bin/bash

# Define backup directory (you can change this path)
BACKUP_DIR="/home/cfmi.arpal.org/satsuolo/QNAPDEV_SaturazioneSuolo/info/crontab/"

# Create the directory if it doesn't exist
mkdir -p "$BACKUP_DIR"

# Get current date and time
DATE=$(date +%Y-%m-%d)

# Define the backup file name
BACKUP_FILE="$BACKUP_DIR/${DATE}_crontab_backup.txt"

# Save the crontab to the file
crontab -l > "$BACKUP_FILE"

# Confirm the backup
echo "Crontab saved to $BACKUP_FILE"

