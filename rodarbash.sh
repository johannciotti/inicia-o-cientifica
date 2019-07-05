#!/bin/bash

for i in {1..8}; do 
  gcc -O3 -o x$i.x IandF_johann_g4_eletric_bash.c -lm
  echo 
done

for i in {1..8}; do 
   echo $i | srun -n1 ./x$i.x & 
   echo 
   sleep 1
   rm x$i.x 
