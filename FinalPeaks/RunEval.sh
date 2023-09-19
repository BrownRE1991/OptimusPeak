#!/bin/bash

if [ -f "Combined_Results.csv" ]; then
	rm Combined_Results.csv
fi

python3 Combine.py 

