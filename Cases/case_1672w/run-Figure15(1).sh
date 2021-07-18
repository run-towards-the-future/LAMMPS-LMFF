##Figure 15: Weak Scaling with large processes
 #Figure 15(1): Nlocal = 5542

 bsub -cache_size 0  -o out.1672w.n3072SX1		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 3072    	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 1	 	
 bsub -cache_size 0  -o out.1672w.n6144SX2		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6144	  	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 2	 	
 bsub -cache_size 0  -o out.1672w.n12288SX4		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 12288	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 4	 	
 bsub -cache_size 0  -o out.1672w.n24576SX8	 	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 24576	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8	 	
 bsub -cache_size 0  -o out.1672w.n49152SX16	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 49152	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16	
 bsub -cache_size 0  -o out.1672w.n98304SX32   	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 98304	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 32	
 bsub -cache_size 0  -o out.1672w.n196608SX64  	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 196608	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 64	

