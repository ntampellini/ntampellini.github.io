#!/bin/bash



# Function to display file sizes in human-readable format and calculate the total size

display_sizes_and_sum() {

    local total_size=0

    echo "File Sizes:"    



    for file in "$@"; do

        if [[ -f "$file" ]]; then

            # Get the size of the file in bytes

            size=$(stat --printf="%s" "$file")

            total_size=$((total_size + size))



            # Convert the size to a human-readable format

            readable_size=$(ls -lh "$file" | awk '{print $5}')

            echo "$file: $readable_size"

        else

            echo "$file: [File does not exist or is not a regular file]"

        fi

    done



    # Convert the total size to a human-readable format

    readable_total=$(numfmt --to=iec $total_size)

    echo -e "\nTotal Size: $readable_total"

}



# Function to prompt for file deletion

prompt_delete_files() {

    echo -e "\nDo you want to delete these files? (yes/no)"

    read -r response



    if [[ "$response" == "yes" ]]; then

        for file in "$@"; do

            if [[ -f "$file" ]]; then

                rm "$file"

                echo "$file has been deleted."

            fi

        done

    else

        echo "No files were deleted."

    fi

}



# Check if any arguments were passed

if [[ $# -eq 0 ]]; then

    echo "Usage: $0 [file1 file2 ... | -m match_string1 [match_string2 ...]]"

    exit 1

fi



# Handle match string input

if [[ "$1" == "-m" ]]; then

    if [[ -z "$2" ]]; then

        echo "Error: At least one match string must be provided after -m."

        exit 1

    fi

    shift



    # Use find to gather files matching multiple patterns

    files=""

    for match_string in "$@"; do

        found_files=$(find . -type f -name "$match_string")

        files+="$found_files\n"

    done



    # Remove duplicates and empty lines

    files=$(echo -e "$files" | sort -u | sed '/^$/d')



    if [[ -z "$files" ]]; then

        echo "No files found matching the given patterns."

        exit 0

    fi



    # Convert files to an array for processing

    set -- $files

fi



# Display file sizes and total size

display_sizes_and_sum "$@"



# Prompt the user for file deletion

prompt_delete_files "$@"


