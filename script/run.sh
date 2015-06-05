#!/bin/sh
mkdir -p result
dataset="kdd"
#dataset="synthetic"
T=11
P=10
Q=5
R1=100
R2=10000
K=20
theta=0.01
./blesimrank $dataset $T $P $Q $R1 $R2 $K $theta | tee info.log |ccze 
