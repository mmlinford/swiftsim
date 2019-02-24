#!/bin/bash

i=0
for file in $( ls stellar_evolution_0*.hdf5 )
do
  echo "plotting data for snapshot $i from $file"
  python plot_sedov.py $file
  cp sedov.png "sedov_$i.png"
  i=$((i + 1))
done
