#!/bin/bash

# Check if both input and output files are provided as arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

# Input file
input_file="$1"

# Output file
output_file="$2"

# Filter the rows based on the specified conditions using awk
awk -F '|' '{
    # Define column names
    name=$1; ra=$2; dec=$3; sigra=$4; sigdec=$5;
    w1=$6; w1sigpro=$7; w2=$8; w2sigpro=$9; w3=$10;
    w1flux=$11; w1sigflux=$12; w2flux=$13; w2sigflux=$14; w3flux=$15; w3sigflux=$16;

    # Check conditions
    if (w1 < 14.3 && w2 < 13.8 && w3 < 12.2 && (w1 - w2) > 0.05 && (w1 - w2) < 0.9 && (w2 - w3) > 1.25 && (w2 - w3) < 2.4) {
        print $0 > "'"$output_file"'";
    }
}' "$input_file"


