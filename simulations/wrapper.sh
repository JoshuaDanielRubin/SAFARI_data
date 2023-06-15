#!/bin/bash

# Define the filename of the Python script
python_script="simulate.py"

# Run the Python script with varying parameters
for N in 500 1000 1500
do
    for L in 80 100 120
    do
        for delta in 0.01 0.02 0.03
        do
            for k in 10 15 20
            do
                for W in 20 30 40
                do
                    if [ $W -ge $k ]; then
                        echo "Running with N=$N, L=$L, delta=$delta, k=$k, W=$W"
                        python3 $python_script --N $N --L $L --delta $delta --k $k --W $W
                    fi
                done
            done
        done
    done
done

