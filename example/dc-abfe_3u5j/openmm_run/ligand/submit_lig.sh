#!/bin/bash
# Find directories matching the patterns and store the results in an array
directories=(`find . -maxdepth 1 -type d \( -name 'sterics*' \)`)

# Function to extract the numeric part for sorting
extract_number() {
    echo "$1" | awk 'match($0, /([0-9]+)$/, a) {print a[1] + 0}'
}

sorted_dirs=()
for dir in "${directories[@]}";do
    basename=$(basename "$dir")
    number=$(extract_number "$basename")
    if [[ $basename == restraints* ]]; then
        sorted_dirs+=(`printf "1~%s~%012d\n" "$basename" "$number"`)
    elif [[ $basename == electrostatics* ]]; then
        sorted_dirs+=(`printf "2~%s~%012d\n" "$basename" "$number"`)
    elif [[ $basename == sterics* ]]; then
        sorted_dirs+=(`printf "3~%s~%012d\n" "$basename" "$number"`)
    fi
done

IFS=$'\n' sorted_dirs=($(sort -t"~" -k1,1n -k3,3n <<< "${sorted_dirs[*]}"))
unset IFS


# Output the sorted results
for dir in "${sorted_dirs[@]}"; do
    sorted_single_dir=`echo "$dir" | cut -d~ -f2`
    cd $sorted_single_dir
    sh md_run.sh
    cd ..
done

