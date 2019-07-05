#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

source /opt/intel/Compiler/11.1/064/bin/iccvars.sh intel64
icc <code.c> -o IandF_johann_g4_eletric_bash.c
 ./IandF_johann_g4_eletric_bash.c
 
