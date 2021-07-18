##Figure 17: Strong Scaling with large processes.

 #Figure 17(1):1020M expr 
 #1,070,080,000 atoms
 bsub -cache_size 0  -o out.13376wn3072SX8		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 3072   	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8
 bsub -cache_size 0  -o out.13376wn6144SX8		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6144	 	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8
 bsub -cache_size 0  -o out.13376wn12288SX8 	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 12288	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8
 bsub -cache_size 0  -o out.13376wn24576SX8	  	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 24576	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8
 bsub -cache_size 0  -o out.13376wn49152SX8  	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 49152	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8
 bsub -cache_size 0  -o out.13376wn98304SX8 	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 98304	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8
 bsub -cache_size 0  -o out.13376wn196608SX8    -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 196608	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8


 #Figure 17(2):2041M expr 
 #2,140,160,000 atoms
 bsub -cache_size 0  -o out.13376wn3072SX16		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 3072   	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16
 bsub -cache_size 0  -o out.13376wn6144SX16		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6144	  	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16
 bsub -cache_size 0  -o out.13376wn12288SX16 	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 12288	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16
 bsub -cache_size 0  -o out.13376wn24576SX16	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 24576	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16
 bsub -cache_size 0  -o out.13376wn49152SX16  	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 49152	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16
 bsub -cache_size 0  -o out.13376wn98304SX16 	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 98304	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16
 bsub -cache_size 0  -o out.13376wn196608SX16 	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 196608	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16

 
 

