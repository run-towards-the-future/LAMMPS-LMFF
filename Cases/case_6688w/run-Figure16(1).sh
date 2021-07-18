##Figure 16: Weak Scaling with small processes
 #Figure 16(1): Nlocal = 1360M 
 
 bsub -cache_size 0  -o out.6688wn48SX1			-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 48   -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 1
 bsub -cache_size 0  -o out.6688wn96SX2			-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 96   -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 2
 bsub -cache_size 0  -o out.6688wn192SX4		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 192  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 4
 bsub -cache_size 0  -o out.6688wn384SX8		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 384  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 8
 bsub -cache_size 0  -o out.6688wn768SX16	  	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 768  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 16
 bsub -cache_size 0  -o out.6688wn1536SX32		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1536 -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate -var SX 32


