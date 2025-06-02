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

for file in ./src/*.cc; do
   echo "Processing file: $file"
   clang-format -i -style=file $file
done

for file in ./inc/honeycomb2/*.hpp; do
   echo "Processing file: $file"
   clang-format -i -style=file $file
done

for file in ./tests/*.hpp; do
   echo "Processing file: $file"
   clang-format -i -style=file $file
done

for file in ./tests/*.cc; do
   echo "Processing file: $file"
   clang-format -i -style=file $file
done
