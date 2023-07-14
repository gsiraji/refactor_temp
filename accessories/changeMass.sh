#!/bin/bash

for (( k=1; k<=10; k++ ))
do

    sed "s/mass,[^,]*/mass,1e-$k/" inertia_sim_parameters.csv > "massIs$k/inertia_sim_parameters.csv"

done
