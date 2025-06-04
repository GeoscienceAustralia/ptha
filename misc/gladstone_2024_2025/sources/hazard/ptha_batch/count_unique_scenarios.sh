#!/bin/bash

dirs=random_*/scenario_initial_conditions

for dir in $dirs; do
    ls $dir | wc -l
done
 