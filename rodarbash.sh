#!/bin/bash
for i in {1..2}; do 
  gcc -O3 -o x$i.x IandF_johann_g4_eletric_bash.c -lm
  echo 
done
#!/bin/bash

#$ -cwd
#$ -j y
#$ -S /bin/bash
#

source /opt/intel/Compiler/11.1/064/bin/iccvars.sh intel64
icc <code.c> -o ./x$i.x &
 ./x$i.x &
 
done

