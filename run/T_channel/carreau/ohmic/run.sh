#!/bin/bash


cd /home/skotapa/openfoam-tests/C_slip_50_PPM_70V_1
. /opt/openfoam9/etc/bashrc
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH

./Allclean

blockMesh


cp 0/C.org 0/C
cp  0/sigma.org 0/sigma

setFields
checkMesh -allGeometry
decomposePar

mpirun rheoEFoam -parallel