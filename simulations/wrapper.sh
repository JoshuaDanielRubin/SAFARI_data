#!/bin/bash

# Define the filename of the Python script
python_script="main.py"

# Run the Python script with varying parameters
for W in 20
do
    for L in 100
    do
        for delta in 0.01 0.5
        do
            for k in 10 15
            do
                for N in 100 500 1000 1500 3000
                do
                    for unique_flag in "" "--unique"
                    do
                        if [ $W -ge $k ]; then
                            echo "Running with N=$N, L=$L, delta=$delta, k=$k, W=$W $unique_flag"
                            python3 $python_script --N $N --L $L --delta $delta --k $k --W $W $unique_flag
                        fi
                    done
                done
            done
        done
    done
done

