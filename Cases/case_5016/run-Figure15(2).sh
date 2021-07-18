##Figure 15: Weak Scaling with large processes
 #Figure 15(2): Nlocal = 3344
 bsub -cache_size 0  -o out.5016k.n3072SX2048		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 3072   -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 2048	 	
 bsub -cache_size 0  -o out.5016k.n6144SX4096		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6144	  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 4096	 	
 bsub -cache_size 0  -o out.5016k.n12288SX8192		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 12288  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8192	 	
 bsub -cache_size 0  -o out.5016k.n24576SX16384	 	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 24576  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16384	 	
 bsub -cache_size 0  -o out.5016k.n49152SX32768		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 49152  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 32768	
 bsub -cache_size 0  -o out.5016k.n98304SX65536		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 98304  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 65536	
 bsub -cache_size 0  -o out.5016k.n196608SX131072	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 196608 -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 131072	

