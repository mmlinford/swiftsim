#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e gravity_glassCube_32.hdf5 ]
then
    echo "Fetching initial gravity glass file for the constant cosmological box example..."
    ./getGlass.sh
fi

# Fetch the cooling tables
if [ ! -e ics_no_z.hdf5 ]
then
    echo "Generating initial conditions for the uniform cosmo box example..."
    python3 makeIC.py
fi

swift_location="../../../build/anarchy-const-lambda/examples/swift"

# Run SWIFT
$swift_location --hydro --cosmology --cooling --threads=4 cooling_redshift_dependence_low_z.yml 2>&1 | tee output_low.log
$swift_location --hydro --cooling --threads=4 cooling_redshift_dependence_no_z.yml 2>&1 | tee output_no.log

