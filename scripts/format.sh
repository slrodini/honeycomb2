#!/usr/bin/env bash

if [ -z "$(which clang-format)" ]; then
   echo -e "clang-format is not install... Aborting."
   exit -1
fi

major_version=$(clang-format --version | grep -oE '[0-9]+' | head -n1)
if [ "$major_version" -le 15 ]; then
   echo "Clang-format version $major_version is older than 16. Aborting..."
   exit -2
fi

DIRECTORY="../src/"
EXTENSION="cc"

files=$(find "$DIRECTORY" -type f -name "*.$EXTENSION")

if [ -z "$files" ]; then
   echo "No files with the extension .$EXTENSION found in $DIRECTORY"
else
   for file in $files; do
      echo "Processing file: $file"
      clang-format -i -style=file $file
   done
fi

DIRECTORY="../inc/honeycomb2/"
EXTENSIONS=("h" "hpp" "hh")

# Dynamically build the find command
find_cmd="find \"$DIRECTORY\" -type f \\( "
for ext in "${EXTENSIONS[@]}"; do
   find_cmd+=" -name \"*.$ext\" -o"
done
find_cmd=${find_cmd% -o} # Remove trailing '-o'
find_cmd+=" \\)"

# Execute and store results
files=$(eval "$find_cmd")

if [ -z "$files" ]; then
   echo "No files with the extension .$EXTENSION found in $DIRECTORY"
else
   for file in $files; do
      echo "Processing file: $file"
      clang-format -i -style=file $file
   done
fi
