#!/bin/bash

# Add path to binary
export PATH=$PATH:"$PWD/../bin"

# Generate 2d grid
gridgenlaval2d \
    -a 0.25 \
    -d 0.5 \
    -x1 5.0 \
    -y1 8.205 \
    -h0 0.0 \
    -h1 -2.705 \
    -h2 0.0 \
    -h3 25.0 \
    -r0 5.5 \
    -r1 5.0 \
    -w1 10.0 \
    -w2 0.0 \
    -w3 50 \
    -a1 90.0 \
    -a2 0.674785758658841 \
    -l  50.0 \
    -f nozzle1

showme nozzle1.1
