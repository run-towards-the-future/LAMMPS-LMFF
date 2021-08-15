##Figure 9: Speedup on one CG or one core of different platforms. 
##The running of Baseline and LMFF on SW processors just are on MPE of one CG, while others use an Intel core. 

##SW39000(Sunway SW39000)
 #Baseline
  bsub -cache_size 0  -o out.SW9.1CG.Baseline	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff
 #LMFF	 	
  bsub -cache_size 0  -o out.SW9.1CG.LMFF	    -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1 -b -cgsp 64  ../../LMFF/lmp_sunway_big -in in.replicate	

##SW26010(Sunway SW26010)
 #Baseline
  bsub -o out.SW5.1CG.Baseline	-q q_sw_expr -host_stack 256 -share_size 6144 -priv_size 4 -n 1 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff
 #LMFF	 	
  bsub -o out.SW5.1CG.LMFF	    -q q_sw_expr -host_stack 256 -share_size 6144 -priv_size 4 -n 1 -b -cgsp 64  ../../LMFF/lmp_sunway_big -in in.replicate
 	
##G6278C(Intel Xeon Gold 6278C)
 #Baseline
  mpirun -n 1 ../../Baseline/lmp_intel_cpu_intelmpi -in in.replicate.ilp_tersoff 2>1|tee out.G6278C.1Core.Baseline
 #LMFF
  mpirun -n 1 ../../LMFF/lmp_intel_cpu_intelmpi -in in.replicate 2>1|tee out.G6278C.1Core.LMFF

##G6148(Intel Xeon Gold 6148)
#Baseline
  mpirun -n 1 ../../Baseline/lmp_intel_cpu_intelmpi -in in.replicate.ilp_tersoff 2>1|tee out.G6148.1Core.Baseline
 #LMFF
  mpirun -n 1 ../../LMFF/lmp_intel_cpu_intelmpi -in in.replicate 2>1|tee out.G6148.1Core.LMFF

##E5(Intel Xeon E5 2680-v3 on Sunway Supercomputer)
#Baseline
  bsub -o out.E5.1Core.Baseline -q q_x86_expr -n 1 -b ../../Baseline/lmp_intel_cpu_intelmpi -in in.replicate.ilp_tersoff
 #LMFF	 	
  bsub -o out.E5.1Core.LMFF	    -q q_x86_expr -n 1 -b ../../LMFF/lmp_intel_cpu_intelmpi -in in.replicate
