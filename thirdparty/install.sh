#!/bin/bash

# OpenBLAS
wget https://github.com/xianyi/OpenBLAS/archive/v0.3.6.tar.gz
mv v0.3.6.tar.gz OpenBLAS-0.3.6.tar.gz
tar -xf OpenBLAS-0.3.6.tar.gz
mv OpenBLAS-0.3.6 OpenBLAS
(cd OpenBLAS && make)

# SuiteSparse
wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.4.0.tar.gz
tar -xf SuiteSparse-5.4.0.tar.gz
(cd SuiteSparse && make metis)
(cd SuiteSparse/SuiteSparse_config && make)
(cd SuiteSparse/AMD && make)
(cd SuiteSparse/CAMD && make)
(cd SuiteSparse/COLAMD && make)
(cd SuiteSparse/CCOLAMD && make)
(cd SuiteSparse/CHOLMOD && make)
(cd SuiteSparse/UMFPACK && make)
