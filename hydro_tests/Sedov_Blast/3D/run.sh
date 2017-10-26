#! /bin/bash

if [ ! -e glassCube_64.hdf5 ]
then
  echo "Fetching initial glass files for the 3D Sedov blast tests..."
  ./getGlass.sh
fi

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Sedov_Blast/3D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Sedov_Blast/3D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 16 sedov.yml
    mv sedov_0005.hdf5 $folder/
    mv timesteps_16.txt $folder/
    python analyze.py $folder
  done
done