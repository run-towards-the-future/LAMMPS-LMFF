##Table5: The Time Ratio of Hotspots on One CG of SW39000 processor

##SW39000(Sunway SW39000)
 #Baseline
  bsub -cache_size 0  -o out.SW9.1CG.Baseline	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff
 #Rft
  bsub -cache_size 0  -o out.SW9.1CG.Rft		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1 -b -cgsp 64  ../../Rft/lmp_sunway_big -in in.replicate.ilp_tersoff
 #LMFF	 	
  bsub -cache_size 0  -o out.SW9.1CG.LMFF	    -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1 -b -cgsp 64  ../../LMFF/lmp_sunway_big -in in.replicate	