##Figure 11: Speedup on one plugin board. Baseline and LMFF just use 8×6 MPEs, while SWLMFF can run on 8×6 MPEs and 8 × 6 × 64 CPEs. 
##The 16M, 32M, and 64M in y-axis are three atomic systems: 16,720,000 atoms (case-1672w), 33,440,000 atoms (case-3344w), and 66,880,000 atoms (case-6688w), respectively.
 #Figure 11(1): 64M

##SW39000(Sunway SW39000)
 #Baseline
  bsub -cache_size 0  -o out.SW9.1Plugin.Baseline	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 48 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff
 #LMFF	 	
  bsub -cache_size 0  -o out.SW9.1Plugin.LMFF		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 48 -b -cgsp 64  ../../LMFF/lmp_sunway_big -in in.replicate
 #SWLMFF	 	
  bsub -cache_size 0  -o out.SW9.1Plugin.SWLMFF		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 48 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate	
