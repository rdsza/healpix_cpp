#!/bin/bash

# Remove all files and subdirectories in folders matching n_side*
for dir in nside_*/; do
    if [ -d "$dir" ]; then
        echo "Cleaning up $dir"
        rm -rf "$dir"/*
    fi
done

echo "Cleanup complete"