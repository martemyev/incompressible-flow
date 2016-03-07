#!/bin/bash

# This script extrudes 2D media properties copying a binary file. The
# properties are the same in XY-plane which is not very suitable for the
# geophysical application (we typically want them to be copyied in XZ-plane).
# See extrude2D.m script for that.

for i in `seq 1 50`; do
  #echo "Extruded ${i} times"
  cat ../build/spe10_bottom.bin >> ../build/spe10_bottom_extruded_50.bin
done
