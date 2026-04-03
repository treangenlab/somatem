#!/bin/bash
# No point copying the dirs; let's just use them from Source if Eddy gives permission


# Script to copy unified databases from source to destination
# Author: Generated for Somatem pipeline
# Date: $(date +%Y-%m-%d)

# Configuration
SOURCE_DB_PATH="/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025"
DEST_PATH="assets/databases"

# List of directories to copy (add more to this array as needed)
DIRS_TO_COPY=(
    "sylph_abf_030325"
    "ganon2_abvf_030325" 
    "k2_abfv_030325"
)

# Create destination directory if it doesn't exist
mkdir -p "$DEST_PATH"

echo "Starting database copy process..."
echo "Source: $SOURCE_DB_PATH"
echo "Destination: $DEST_PATH"
echo "----------------------------------------"

# Function to copy directory if it doesn't exist
copy_directory() {
    local dir_name="$1"
    local source_dir="$SOURCE_DB_PATH/$dir_name"
    local dest_dir="$DEST_PATH/$dir_name"
    
    echo "Processing: $dir_name"
    
    # Check if source directory exists
    if [ ! -d "$source_dir" ]; then
        echo "  ❌ ERROR: Source directory '$source_dir' does not exist!"
        return 1
    fi
    
    # Check if destination already exists
    if [ -d "$dest_dir" ]; then
        echo "  ⏭️  SKIP: Directory '$dir_name' already exists in destination"
        return 0
    fi
    
    # Copy the directory
    echo "  📁 COPYING: $source_dir -> $dest_dir"
    cp -r "$source_dir" "$dest_dir"
    
    if [ $? -eq 0 ]; then
        echo "  ✅ SUCCESS: '$dir_name' copied successfully"
        return 0
    else
        echo "  ❌ ERROR: Failed to copy '$dir_name'"
        return 1
    fi
}

# Main execution
success_count=0
error_count=0
skip_count=0

for dir in "${DIRS_TO_COPY[@]}"; do
    copy_directory "$dir"
    case $? in
        0) success_count=$((success_count + 1)) ;;
        1) error_count=$((error_count + 1)) ;;
        2) skip_count=$((skip_count + 1)) ;;
    esac
    echo ""
done

echo "----------------------------------------"
echo "Copy process completed!"
echo "✅ Successfully copied: $success_count"
echo "⏭️  Skipped (already exist): $skip_count" 
echo "❌ Errors: $error_count"

if [ $error_count -gt 0 ]; then
    echo "⚠️  Some directories failed to copy. Please check the errors above."
    exit 1
else
    echo "🎉 All operations completed successfully!"
    exit 0
fi
