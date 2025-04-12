#!/bin/bash
# cleanup_legacy_dirs.sh - Script to safely remove legacy directories after migration

# Set color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
RESET='\033[0m'

echo -e "${BOLD}Desmognathus_TE Legacy Directory Cleanup${RESET}"
echo -e "${YELLOW}WARNING: This script will permanently delete legacy directories.${RESET}"
echo -e "${YELLOW}Make sure you have completed migration and backed up any important data.${RESET}"
echo

# List legacy directories
legacy_dirs=(
    "Data"
    "Results"
    "Projects"
    "Output"
    "old_scripts"
)

echo -e "${BOLD}The following legacy directories will be removed:${RESET}"
for dir in "${legacy_dirs[@]}"; do
    if [ -d "$dir" ]; then
        size=$(du -sh "$dir" 2>/dev/null | cut -f1)
        echo -e " - ${YELLOW}$dir${RESET} (size: $size)"
    else
        echo -e " - ${YELLOW}$dir${RESET} (not found)"
    fi
done
echo

# Prompt for confirmation
read -p "Do you want to proceed with removal? (y/N): " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Operation cancelled.${RESET}"
    exit 0
fi

echo -e "${YELLOW}Creating backup archive of legacy directories...${RESET}"
backup_file="legacy_backup_$(date +%Y%m%d_%H%M%S).tar.gz"
tar -czf "$backup_file" "${legacy_dirs[@]}" 2>/dev/null
echo -e "${GREEN}Backup created: $backup_file${RESET}"

# Remove directories
echo -e "${YELLOW}Removing legacy directories...${RESET}"
for dir in "${legacy_dirs[@]}"; do
    if [ -d "$dir" ]; then
        rm -rf "$dir"
        if [ ! -d "$dir" ]; then
            echo -e " - ${GREEN}Removed: $dir${RESET}"
        else
            echo -e " - ${RED}Failed to remove: $dir${RESET}"
        fi
    else
        echo -e " - ${YELLOW}Skipped: $dir (not found)${RESET}"
    fi
done

echo
echo -e "${GREEN}Cleanup complete!${RESET}"
echo -e "${YELLOW}If you need to restore the backup, use:${RESET}"
echo -e "  tar -xzf $backup_file"
echo
echo -e "You can now run the following to create .gitkeep files for essential directories:"
echo -e "  python scripts/python/utils/create_gitkeep.py" 