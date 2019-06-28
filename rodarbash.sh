#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

icc <code.c> -o  	IandF_johann_g4_eletric_bash1.c
 ./ 	IandF_johann_g4_eletric_bash1.c
 
