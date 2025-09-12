#!/bin/bash

CMAKELISTS="CMakeLists.txt"
COMPILE=1

if [ $COMPILE -eq 1 ]; then
    echo "Compiling"
# Loop over different NS, NT, and LEVELS vales
  for NS in 256; do
  for LEVELS in 2 3 4; do
    NT=$NS
    # Update NS and NT (lines 22 and 23)
    sed -i "22s/set(NS \".*\")/set(NS \"${NS}\")/" "$CMAKELISTS"
    sed -i "23s/set(NT \".*\")/set(NT \"${NT}\")/" "$CMAKELISTS"
    # Update LEVELS (line 27)
    sed -i "27s/set(LEVELS \".*\")/set(LEVELS \"${LEVELS}\")/" "$CMAKELISTS"
    echo "--------Compiling for NS=${NS}, NT=${NT}, LEVELS=${LEVELS}--------"
    mkdir -p build
    cd build
    cmake ..
    make -j 8
    cd ..
  done
  done
  echo "All builds completed."
else
  echo "Compilation skipped."
fi

#run program 

for NS in 256; do
    NT=$NS
    rm parameters.dat
    printf "0 4 4 20 4 4" >> parameters.dat
    echo "--------Running for NS=${NS}, NT=${NT}, LEVELS=2--------"
    cd build
    mpirun -np 8 AMG_${NS}x${NT}_l2
    cd ..

    rm parameters.dat
    printf "0 8 8 20 4 4\n1 4 4 20 4 4" >> parameters.dat
    echo "--------Running for NS=${NS}, NT=${NT}, LEVELS=3--------"
    cd build
    mpirun -np 8 AMG_${NS}x${NT}_l3
    cd ..


    rm parameters.dat
    printf "0 8 8 20 4 4\n1 4 4 20 4 4\n2 2 2 20 2 2" >> parameters.dat
    echo "--------Running for NS=${NS}, NT=${NT}, LEVELS=4--------"
    cd build
    mpirun -np 8 AMG_${NS}x${NT}_l4
    cd ..


done
