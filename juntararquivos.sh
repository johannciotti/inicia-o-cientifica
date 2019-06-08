for i in {1..8}; do
   let "j=$i-1" 
   cat POK_CV_F_N100_p_gexc_g1_gel0.1_-$j.dat POK_CV_F_N100_p_gexc_g1_gel0.1_$i.dat > POK_CV_F_N100_p_gexc_g1_gel0.1_-$i.dat
   echo 
   rm POK_CV_F_N100_p_gexc_g1_gel0.1_$i.dat
   rm POK_CV_F_N100_p_gexc_g1_gel0.1_-$j.dat
done
