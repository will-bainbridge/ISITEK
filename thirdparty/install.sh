#!/bin/bash

# source

if [ ! -d GotoBLAS2 ]
then
	echo "download GotoBLAS from <http://www.tacc.utexas.edu/tacc-projects/gotoblas2> and extract to this directory"
	exit
fi

wget http://www.cise.ufl.edu/research/sparse/UFconfig/current/UFconfig.tar.gz
tar -xf UFconfig.tar.gz
rm -f UFconfig.tar.gz

wget http://www.cise.ufl.edu/research/sparse/amd/current/AMD.tar.gz
tar -xf AMD.tar.gz
rm -f AMD.tar.gz

wget http://www.cise.ufl.edu/research/sparse/umfpack/current/UMFPACK.tar.gz
tar -xf UMFPACK.tar.gz
rm -f UMFPACK.tar.gz

# configure

sed "s|^BLAS.*|BLAS = `pwd`/GotoBLAS2/libgoto2.a|g;s|^LAPACK.*||g;s|^UMFPACK_CONFIG.*|UMFPACK_CONFIG = -DNCHOLMOD|g" UFconfig/UFconfig.mk  > temp
mv -f temp UFconfig/UFconfig.mk

# build

cd GotoBLAS2
make USE_THREAD=0

cd ../UFconfig/
make

cd ../AMD
make

cd ../UMFPACK
make

cd ..
