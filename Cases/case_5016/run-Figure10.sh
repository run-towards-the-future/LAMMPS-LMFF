##Figure 10: Speedup on one node. Baseline and LMFF just use 6 MPEs of one node, while SWLMFF can run on 6 MPEs and 6 × 64 CPEs.

##SW39000(Sunway SW39000)
 #Baseline
  bsub -cache_size 0  -o out.SW9.1Node.Baseline	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff
 #LMFF	 	
  bsub -cache_size 0  -o out.SW9.1Node.LMFF	    -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6 -b -cgsp 64  ../../LMFF/lmp_sunway_big -in in.replicate
 #SWLMFF	 	
  bsub -cache_size 0  -o out.SW9.1Node.SWLMFF	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate	