#!/bin/sh

#setFields

#for scaling down the coordinate system by 1000
#transformPoints "scale=(0.001 0.001 0.001)"

decomposePar

mpirun -np 4 rhoCentralFoam -parallel > log.rhoCentralFoam

reconstructPar
