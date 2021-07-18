##Figure 18: Strong Scaling with small processes
 #64M expr
 #66,880,000 atoms
 bsub -cache_size 0  -o out.6688w.n48		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 48	  	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.6688w.n96		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 96	  	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.6688w.n192		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 192	  	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.6688w.n384		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 384	  	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.6688w.n768		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 768	  	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.6688w.n1536		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1536		-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.6688w.n3072		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 3072		-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.6688w.n6144		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6144		-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 


