#!/bin/bash

################################################################################
# OpenFOAM AD Library Renaming and Patching Script
################################################################################
#
# PURPOSE:
#   Renames all OpenFOAM AD libraries from lib*.so to lib*ADR.so or lib*ADF.so
#   and updates their internal NEEDED references using patchelf. Also patches
#   NEEDED references in executables located in the bin/ subdirectory.
#
# USAGE:
#   Dry-run mode (shows what will be changed):
#     ./rename_ad_libs.sh /path/to/platform/root --ADR
#     ./rename_ad_libs.sh /path/to/platform/root --ADF
#
#   Commit mode (applies changes):
#     ./rename_ad_libs.sh /path/to/platform/root --ADR --commit
#     ./rename_ad_libs.sh /path/to/platform/root --ADF --commit
#
# REQUIREMENTS:
#   - patchelf (must be installed separately)
#   - Bash 4.0+
#   - Standard Unix tools: find, sed, basename, dirname, mv
#   - Platform directory must contain lib/ and bin/ subdirectories
#
# EXAMPLES:
#   Test on your AD platform (ADR suffix):
#     ./rename_ad_libs.sh OpenFOAM-AD/platforms/linux64GccDPInt32OptADR --ADR
#
#   Apply the changes with ADR suffix:
#     ./rename_ad_libs.sh OpenFOAM-AD/platforms/linux64GccDPInt32OptADR --ADR --commit
#
#   Test with ADF suffix:
#     ./rename_ad_libs.sh OpenFOAM-AD/platforms/linux64GccDPInt32OptADF --ADF
#
#   Apply the changes with ADF suffix:
#     ./rename_ad_libs.sh OpenFOAM-AD/platforms/linux64GccDPInt32OptADF --ADF --commit
#
################################################################################

set -e

# ============================================================================
# Configuration and Input Validation
# ============================================================================

# Parse command-line arguments
LIB_DIR=""
SUFFIX=""
DRY_RUN=true

# Process arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ADR)
            SUFFIX="ADR"
            shift
            ;;
        --ADF)
            SUFFIX="ADF"
            shift
            ;;
        --commit)
            DRY_RUN=false
            shift
            ;;
        *)
            # First non-flag argument is the library directory
            if [ -z "$LIB_DIR" ]; then
                LIB_DIR="$1"
            fi
            shift
            ;;
    esac
done

# Use current directory if not specified
LIB_DIR="${LIB_DIR:-.}"

# Validate that suffix was specified
if [ -z "$SUFFIX" ]; then
    print_error() {
        echo -e "\033[0;31m$1\033[0m"
    }
    print_error "Error: You must specify either --ADR or --ADF suffix"
    echo
    echo "Usage: $0 /path/to/lib/directory --ADR|--ADF [--commit]"
    echo
    echo "Examples:"
    echo "  $0 /path/to/lib --ADR              (dry-run with ADR suffix)"
    echo "  $0 /path/to/lib --ADF              (dry-run with ADF suffix)"
    echo "  $0 /path/to/lib --ADR --commit     (apply changes with ADR suffix)"
    echo "  $0 /path/to/lib --ADF --commit     (apply changes with ADF suffix)"
    echo
    exit 1
fi

# ANSI color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'  # No Color

# ============================================================================
# Helper Functions
# ============================================================================

# Print info message in blue
print_info() {
    echo -e "${BLUE}$1${NC}"
}

# Print success message in green
print_success() {
    echo -e "${GREEN}$1${NC}"
}

# Print warning message in yellow
print_warning() {
    echo -e "${YELLOW}$1${NC}"
}

# Print error message in red
print_error() {
    echo -e "${RED}$1${NC}"
}

# ============================================================================
# Validation
# ============================================================================

print_info "=== OpenFOAM AD Library Renaming and Patching Script ==="
echo "Platform root directory: $LIB_DIR"
echo "Suffix mode: lib*${SUFFIX}.so"
echo

# Check if platform directory exists
if [ ! -d "$LIB_DIR" ]; then
    print_error "Error: Directory '$LIB_DIR' does not exist"
    exit 1
fi

# Construct paths to lib and bin directories
ACTUAL_LIB_DIR="$LIB_DIR/lib"
BIN_DIR="$LIB_DIR/bin"

# Check if lib directory exists
if [ ! -d "$ACTUAL_LIB_DIR" ]; then
    print_error "Error: Library directory '$ACTUAL_LIB_DIR' does not exist"
    exit 1
fi

# Check if bin directory exists
if [ ! -d "$BIN_DIR" ]; then
    print_error "Error: Binary directory '$BIN_DIR' does not exist"
    exit 1
fi

# Check if patchelf is available
if ! command -v patchelf &> /dev/null; then
    print_error "Error: patchelf is not installed"
    echo
    echo "Please install patchelf first:"
    echo "  Ubuntu/Debian:  sudo apt-get update && sudo apt-get install -y patchelf"
    echo "  CentOS/RHEL:    sudo yum install -y patchelf"
    echo "  macOS:          brew install patchelf"
    echo "  Other systems:  https://github.com/NixOS/patchelf"
    echo
    exit 1
fi

# Display mode (DRY_RUN is already set during argument parsing)
if [ "$DRY_RUN" = true ]; then
    print_info "Running in DRY-RUN mode - no changes will be applied"
else
    print_warning "Running in COMMIT mode - changes will be applied!"
fi
echo

# ============================================================================
# Step 1: Find all libraries and create mapping
# ============================================================================

print_info "Step 1: Finding all lib*.so files and creating library mapping..."

SO_FILES=()
declare -A LIB_MAPPING

# First, find libraries that need renaming (don't have the target suffix yet)
while IFS= read -r -d '' file; do
    # Skip if already has the target suffix
    if [[ ! "$file" =~ ${SUFFIX}\.so$ ]]; then
        SO_FILES+=("$file")
    fi
done < <(find "$ACTUAL_LIB_DIR" -name "lib*.so*" -type f -print0)

# Build mapping: for each library (renamed or not), map old name to new name
while IFS= read -r -d '' file; do
    basename_file=$(basename "$file")

    # Extract base name without .so suffix or version info
    # Example: libOpenFOAM.so.1.2.3 -> libOpenFOAM
    # Example: libOpenFOAMADR.so -> libOpenFOAM (remove existing suffix)
    base_name=$(echo "$basename_file" | sed "s/${SUFFIX}\.so.*//" | sed 's/\.so.*//')

    # Create mapping of old name -> new name with target suffix
    # Example: libOpenFOAM -> libOpenFOAMADR.so
    new_name="${base_name}${SUFFIX}.so"

    # Map the original (non-suffixed) library name to the new name
    old_name="${base_name}.so"
    LIB_MAPPING["$old_name"]="$new_name"

    # IMPORTANT: Also map the already-suffixed name to itself (for idempotency)
    # This ensures that running the script multiple times doesn't create double-suffixed names
    # Example: libOpenFOAMADR.so -> libOpenFOAMADR.so (no change needed)
    if [[ "$basename_file" =~ ${SUFFIX}\.so ]]; then
        LIB_MAPPING["$basename_file"]="$basename_file"
    fi
done < <(find "$ACTUAL_LIB_DIR" -name "lib*.so*" -type f -print0)

print_success "Found ${#LIB_MAPPING[@]} unique libraries"
print_success "Found ${#SO_FILES[@]} libraries that need renaming"
echo

# ============================================================================
# Step 2: Rename libraries that don't have the target suffix yet
# ============================================================================

if [ ${#SO_FILES[@]} -eq 0 ]; then
    print_info "Step 2: No libraries need renaming (all already have ${SUFFIX} suffix)"
else
    print_info "Step 2: Renaming libraries..."
fi
echo

# ============================================================================
# Step 3: Rename libraries that need renaming
# ============================================================================

RENAMED_COUNT=0

if [ ${#SO_FILES[@]} -gt 0 ]; then
    print_info "Step 3: Renaming libraries..."

    for file in "${SO_FILES[@]}"; do
        basename_file=$(basename "$file")

        # Extract base name without any .so suffix
        base_name=$(echo "$basename_file" | sed 's/\.so.*//')

        # Create new name with target suffix
        new_name="${base_name}${SUFFIX}.so"

        # Extract directory and create new path
        dir_path="${file%/*}"
        new_path="$dir_path/$new_name"

        if [ "$DRY_RUN" = true ]; then
            echo "  [DRY-RUN] mv $file -> $new_path"
        else
            # Show which file we're renaming and where
            echo "  Renaming: $basename_file in $dir_path"
            mv "$file" "$new_path"
            RENAMED_COUNT=$((RENAMED_COUNT + 1))
        fi
    done

    echo
else
    print_info "Step 3: All libraries already have ${SUFFIX} suffix, skipping rename"
    echo
fi

# ============================================================================
# Step 4: Patch NEEDED references using patchelf (libraries)
# ============================================================================

print_info "Step 4a: Fixing library dependencies with patchelf..."

PATCHED_COUNT=0

# Helper function to patch NEEDED entries in a file (library or executable)
patch_needed_entries() {
    local file_path="$1"
    local file_type="$2"  # "library" or "executable"

    print_warning "Processing $file_type: $(basename $file_path)"

    # Get current NEEDED entries
    current_needed=$(patchelf --print-needed "$file_path" 2>/dev/null || echo "")

    if [ -n "$current_needed" ]; then
        echo "  Current dependencies:"
        echo "$current_needed" | head -5 | while read -r lib; do
            echo "    - $lib"
        done
        echo

        # For each NEEDED entry, check if we need to patch it
        while IFS= read -r lib; do
            [ -z "$lib" ] && continue

            # Check if this library is in our mapping
            if [ -n "${LIB_MAPPING[$lib]}" ]; then
                new_lib="${LIB_MAPPING[$lib]}"

                if [ "$DRY_RUN" = true ]; then
                    echo "    [DRY-RUN] Replace: $lib -> $new_lib"
                else
                    echo "    Patching: $lib -> $new_lib"
                    # Use patchelf to replace the NEEDED entry
                    patchelf --replace-needed "$lib" "$new_lib" "$file_path"
                    PATCHED_COUNT=$((PATCHED_COUNT + 1))
                fi
            fi
        done <<< "$current_needed"
    fi

    echo
}

# Patch all renamed libraries
for old_name in "${!LIB_MAPPING[@]}"; do
    new_name="${LIB_MAPPING[$old_name]}"

    # Find all renamed libraries
    while IFS= read -r -d '' lib_path; do
        patch_needed_entries "$lib_path" "library"
    done < <(find "$ACTUAL_LIB_DIR" -name "${new_name}" -type f -print0)
done

# ============================================================================
# Step 4b: Patch NEEDED references in executables
# ============================================================================

print_info "Step 4b: Fixing executable dependencies with patchelf..."

# Find all ELF executables in bin directory (skip scripts, symlinks, etc.)
while IFS= read -r -d '' exe_path; do
    # Only patch if patchelf can read it (i.e., it's a valid ELF binary)
    # patchelf --print-needed returns nothing for non-ELF files or gives an error
    if patchelf --print-needed "$exe_path" &>/dev/null; then
        patch_needed_entries "$exe_path" "executable"
    fi
done < <(find "$BIN_DIR" -type f -print0)

# ============================================================================
# Step 5: Verification (only in commit mode)
# ============================================================================

if [ "$DRY_RUN" = false ]; then
    print_info "Step 5: Verification - Checking patched libraries and executables..."

    # Show sample of patched libraries
    count=0
    for new_name in $(echo "${LIB_MAPPING[@]}" | tr ' ' '\n' | sort -u | head -3); do
        lib_path=$(find "$ACTUAL_LIB_DIR" -name "$new_name" -type f | head -1)
        if [ -f "$lib_path" ]; then
            print_warning "Sample dependencies for library: $(basename $lib_path)"
            patchelf --print-needed "$lib_path" 2>/dev/null | head -3 | while read -r lib; do
                echo "  - $lib"
            done
            echo
            count=$((count + 1))
        fi
    done

    echo "  (showing first 3 libraries, total: ${#LIB_MAPPING[@]})"
    echo

    # Show sample of patched executables
    print_warning "Sample dependencies for executables:"
    count=0
    while IFS= read -r -d '' exe_path; do
        if [ $count -ge 2 ]; then
            break
        fi
        print_warning "  $(basename $exe_path):"
        patchelf --print-needed "$exe_path" 2>/dev/null | head -3 | while read -r lib; do
            echo "    - $lib"
        done
        count=$((count + 1))
    done < <(find "$BIN_DIR" -type f -print0)

    echo
fi

# ============================================================================
# Summary
# ============================================================================

print_success "=== Summary ==="

if [ "$DRY_RUN" = true ]; then
    print_warning "DRY-RUN completed. To apply changes, run:"
    echo "  $0 \"$LIB_DIR\" --${SUFFIX} --commit"
else
    print_success "Successfully renamed $RENAMED_COUNT libraries"
    print_success "Successfully patched $PATCHED_COUNT dependencies"
    echo
    echo "All libraries are now named lib*${SUFFIX}.so with updated internal references"
fi
