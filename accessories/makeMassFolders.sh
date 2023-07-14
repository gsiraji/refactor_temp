#!/bin/bash

for (( k=1; k<=10; k++ ))
do
    # create folders named massIs$k, where k=1 to 10
    foldername="massIs$k"
    mkdir "$foldername"
    
    # copy the parameters.csv file into each folder and modify the mass value
    cp inertia_sim_parameters.csv "$foldername/inertia_sim_parameters.csv"
done



sed "s/mass,[^,]*/mass,1e-1/" inertia_sim_parameters.csv > "massIs1/inertia_sim_parameters.csv"
