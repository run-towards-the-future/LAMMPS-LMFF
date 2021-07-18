##Figure 17: Strong Scaling with large processes.
 #1254M expr 
 
 #Total: 1,314,914,304atoms(=5016 atoms * 262144; also is 1254M = 1314914304/1024/1024)
 bsub -cache_size 0  -o out.5016k.n3072SX262144	  -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 3072	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 262144
 bsub -cache_size 0  -o out.5016k.n6144SX262144	  -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6144	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 262144
 bsub -cache_size 0  -o out.5016k.n12288SX262144  -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 12288	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 262144
 bsub -cache_size 0  -o out.5016k.n24576SX262144  -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 24576	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 262144
 bsub -cache_size 0  -o out.5016k.n49152SX262144  -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 49152	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 262144
 bsub -cache_size 0  -o out.5016k.n98304SX262144  -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 98304	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 262144
 bsub -cache_size 0  -o out.5016k.n196608SX26214  -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 196608	-np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 262144
