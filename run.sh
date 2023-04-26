#!/usr/bin/env bash

# Usage: 
# ./run.sh input_file_path NUM_THREADS

export OMP_NUM_THREADS=$2;

./ksd.out < "$1";

