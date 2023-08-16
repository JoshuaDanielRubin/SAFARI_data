#!/bin/bash

# Define the filename of the Python script
python_script="main.py"

run_task() {
    W=$1
    L=$2
    delta=$3
    k=$4
    N=$5
    unique_flag=$6

    if [ $W -ge $k ]; then
        echo "Running with N=$N, L=$L, delta=$delta, k=$k, W=$W $unique_flag"
        python3 main.py --random-reference --N $N --L $L --delta $delta --k $k --W $W $unique_flag
    fi
}

export -f run_task

nice parallel -j60 run_task ::: 20 50 100 ::: 50 100 200 ::: 0.01 0.3 0.5 0.7 0.9 ::: 10 15 30 ::: 1000 2000 ::: "" "--unique"

