# Layered Materials Force Field (LMFF)


## A.Abstract
* LAMMPS is one of the most popular Molecular Dynamic (MD) packages and is widely used in the field of physics, chemistry and materials simulation. Layered Materials Force Field (LMFF) is our expansion of the LAMMPS potential function based on the Tersoff potential and  inter-layer potential (ILP) in LAMMPS. LMFF is designed to study layered materials such as graphene and boron hexanitride. It is universal and does not depend on any platform. We have also carried out a series of optimizations on LMFF and the optimization work is carried out on the new generation of Sunway supercomputer, called SWLMFF. Experiments show that our implementation is efficient, scalable and portable. When generic LMFF is ported to Intel Xeon Gold 6278C, 2× performance improvement is achieved. For the optimized SWLMFF, the overall performance improvement is nearly 200-330× compared to the original ILP and Tersoff potentials. And SWLMFF has good parallel efficiency of 95%-100% under weak scaling with 2.7 million atoms on a single process. The maximum atomic system simulated by SWLMFF is close to 2^31 atoms. And nanosecond simulations in one day can be realized.

## B. Description
* Artifact name: LMFF
* [![DOI](https://zenodo.org/badge/386895896.svg)](https://zenodo.org/badge/latestdoi/386895896)
* Version: LAMMPS with the stable version of 3Mar2020 from https://download.lammps.org/tars/lammps-3Mar2020.tar.gz
* Key algorithms: Molecular Dynamics
* Compilation: swgcc
* Operating systems and versions: Sunway Raise OS with Linux kernel 3.8.0
* Libraries and versions: MVAPICH2 version 2.2a/SWCH MPI at 20210101
* Execution: By the job system.
* Output: Reported performance in LAMMPS.
* Publicly available?: Yes.

## C.SUMMARY OF THE EXPERIMENTS REPORTED

* In the paper, we evaluate the performance of LMFF and SWLMFF on Sunway (SW) processors and report the speedup on one core group (CG), one node, and one plug-in board. We analyze the strong and weak scalable efficiency and present the performance of the simulation with the time length of one nanosecond. The running of experiments in the paper is as follows:

### 1 Software

* There are four kinds of software in our experiments:

    (1) Baseline: which is the running of Tersoff and ILP potentials in LAMMPS with the stable version of “_3Mar2020_”. URL:
        https://download.lammps.org/tars/lammps-3Mar2020.tar.gz
        
    (2) Rft: which is the refactor of ILP potential compared to the Baseline. URL:
        https://github.com/run-towards-the-future/lammps-lmfftree/master/Rft
        
    (3) LMFF: which is a generic force field for layered materials. It is designed based on Tersoff and ILP potentials in LAMMPS with the stable version of  “_3Mar2020_”. It can run directly on Intel and SW processors. URL:
        https://github.com/run-towards-the-future/lammps-lmfftree/master/LMFF
   
    (4) SWLMFF: which is the optimization of LMFF on Sunway processors, so it can only run on SW processors. URL:
        https://github.com/run-towards-the-future/lammps-lmfftree/master/SWLMFF


### 2 Hardware
* Our experiments are mainly tested on Sunway supercomputer. The parameters of Sunway processors are presented in Table 1:

    **Table 1: The Parameters of Sunway Processors.** 

    | Processor    | SW26010           | SW39000                            |
    |--------------|-------------------|------------------------------------|
    | MPE Freq.    | 1.45GHz           | 2.25GHz                            |
    | Memory       | DDR3 8GB          | DDR4 16GB                          |
    | Cache        | 256KB             | 512KB                              |
    | CPE Freq.    | 1.45GHz           | 2.25GHz                            |
    | Vector width | 256 bit           | 512 bit                            |
    | Rergisters   | Vector × 32  | Scalar × 32 + Vector × 32 |
    | LDM size     | 64KB              | 256KB                              |
    | Cores        | (1+64) × 4   | (1+64) × 6                    |
    | Bandwidth    | 33.6GB/s × 4 | 50GB/s × 6                    |


* Computation resources for TaihuLight can be applied from: http://www.nsccwx.cn/guide/5d301c0724364f0351459268 .
* Computation resources for next-generation Sunway supercomputer are not currently publically available, but it is expected to be available after release. The information about the Sunway supercomputer is as follows:

    * The new generation of Sunway supercomputer system inherits and extends the Sunway TaihuLight architecture based on the new generation of Sunway (SW) high-performance heterogeneous many-core processors and interconnected network chips. Each node on the new generation of Sunway supercomputer system is equipped with a single SW processor (SW39000) that is subdivided into six core groups (CGs). Each CG contains a management processing element (MPE) and 64 computing processing elements (CPEs), and 16 GB attached DDR4 main memory. The main memory can be accessed by both the MPE and the 64 CPEs.
    * The MPE has 32 KB L1 instruction cache and 32 KB L1 data cache, with 512 KB L2 cache for both instructions and data. Each CPE has its own 32 KB L1 instruction cache and 256 KB data storage space. The data storage space can be configured as a local data memory (LDM) that is completely controlled by the user. Part of the data storage space can also be configured as a local data cache (LDCache) automatically managed by the hardware with limited coherence support. Data is transferred between LDM and main memory through direct memory access (DMA). And conventional load/store instructions can also be used to realize data transfer between LDCache and main memory.

### 3 Datasets

* All the datasets are produced by researchers in the field of superlubricity.
Please see the files in  _Cases/case*/_  with the prefix "_in_" and we provide the data download address : https://www.dropbox.com/sh/urg0yf19agzbagq/AADMoFsZY8dqUKvXgeKCvf3Ua?dl=0 .
We provide two kinds of cases for testing:

    (1) Six-layer graphene: case-5016 (5,016 atoms)
    
    (2) Double-layer graphene: case-1672w (16,720,000 atoms), case-3344w (33,440,000 atoms), case-6688W (66,880,000 atoms), and case-13376w (133,760,000 atoms)

### 4 Installation

* Baseline: which is the referenced software. Download the original code of LAMMPS from the official website with the stable version of “_3Mar2020_”.

* Rft, LMFF, and SWLMFF: Pull the source code of Rft, LMFF, and SWLMFF from our repository.


### 5 Build
* You can navigate to the following directories corresponding to each software in turn and compile them:
    ```
     cd Baseline/ 
     cd Rft/    
     cd LMFF/   
     cd SWLMFF/
    ```
    Before building, you can use:  
    ```
    make clean-all  
    ```
    to clean up the previous compilation and recompile the source code of four software in turn:

    (1) Baseline: which is the original code, please refer to the official document of LAMMPS.

    (2) Rft and LMFF: take the Intel Xeon Gold 6148 platform as an example, select the corresponding file from:`LMFF/MAKE/OPTIONS/Makefile.intel_*` provided by LAMMPS depending on the platform and compiler you are running on. For an example:
    ```
    make intel_cpu_intelmpi
    ```
    to compile the source code.
    
    (3) SWLMFF: taking the Sunway processor SW39000 platform as an example, use:
    ```
    make sunway_big
    ```
    to compile the source code.

**Caution** : 
    For SWLMFF, users can turn on or off  **–DSUNWAY**  option by modifying the Makefile file in: `SWLMFF/MAKE/MACHINES/Makefile.sunway_big` and choose to use or not to use optimization techniques.

### 6 Run
When the compilation is done, there will be a binary file named "_lmp*_".

* Baseline, Rft, and LMFF: they can run on many platforms, such as Intel, Sunway and so on. Taking the Intel Xeon Gold 6148 as an example, run using: 

    ```
    mpirun -n <nprocs> ./lmp_intel_cpu_intelmpi -in <InputScript>
    ```
     to start an experiment.
    
* SWLMFF: it just can run on the Sunway platform, using:

    ```
    bsub -cache_size 0 -I -J JobName -b -cgsp 64 -share_size 15000 -host_stack 512 -priv_size 4  -n <nprocs> ./lmp_sunway_big -in <InputScript>
    ```
    to start an experiment.
    You can submit more than one job at a time through the shell script ("_*.sh_") we provide, replacing the option **-I**  with the **-o** and specifying the name of the output file such as "_out.*_". The scripts are also in the directory "_Cases/case$*$/_".
    
**Caution** : 
    Due to the functional limitations of the LSF job management system used in this platform, when the **-o** command is used to save the output content of the run if the specified output file already exists, it will be appended to the original content of the file. Therefore, you can take the last 4 performance values in the output file for analysis. This is the result of this experiment.
    Of course, you can use the **rm** command to clear the history of results before the job is submitted to avoid misunderstandings.

### 7 Performance and Results
* When the job finishes, some files may have been generated. Among, the files "_log.lammps_" or "_out.*_" describe the performance of the software. For example:
    ```
        Performance: 1.140 ns/day, 21.061 hours/ns, 13.189 timesteps/s
    ```
In our experiments, we chose time steps per second ("_timesteps/s_") as the performance measure, because it is more accurate than the other two, "_ns/day_" and "_hours/ns_".

* For all experiments, energy minimization is firstly carried. Then run 200 steps of simulation. Thermodynamic results are output once every 10 steps. Each run consists of four different simulations and each simulation runs 50 steps and evaluates one performance. Thus, there are four different performance values that present different simulations. The performance will be written by LAMMPS in the console and log file named "_log.lammps_". We recommend conducting at least five experiments and taking the average of four simulations for each experiment as an observation.

* In addition to the time-break-down with the MPI timer provided by LAMMPS, we also call the General Purpose Timing Library (GPTL) to sample the time-consuming kernels during our optimization. 
The GPTL results are stored in the file \textit{LammpsTime.0}.

### 8 Experiments for Reproducibility
* We evaluate the performance of LMFF and SWLMFF on SW processors and report the speedup on one CG, one node, and one plug-in board. Please refer to the shell scripts in each directory of "_Cases/case$*$/_". These results are used to draw Figure 9, Figure 10, and Figure 11 in the paper.
* For the strong and weak scalability, we append **-var SX <n>** to the command line arguments to expand the atomic system scale to n times the original size in the x dimension. Please refer to the shell scripts in each directory of "_Cases/case$*$/_". These results are used to draw Figure 15, Figure 16, Figure 17, and Figure 18 in the paper.
\par For the simulation with the time length of one nanosecond, we evaluate the performance of SWLMFF and the referenced performance of LMFF for two cases:

    (1) Six-layer graphene: For the evaluation of 5,016 atoms (case-5016) with one CG (one process), the 200 steps of evaluations, LMFF should have a Timestep/s of ~ 0.48 and 0.042 ns/day, while SWLMFF should have a Timestep/s of ~ 11.5 and 1 ns/day.

    (2) Double-layer graphene: For the evaluation of 16,720,000 atoms (case-1672w) with 1,536 processes (256 nodes), the 200 steps of evaluations, LMFF should have a Timestep/s of ~ 0.7 and an 0.06 ns/day, while SWLMFF should have a Timestep/s of ~ 15 and 1.2 - 1.3 ns/day. 


* After all jobs are finished, results can be found in files, named ``out*". The tables and figures are made in our main paper according to the following steps: 

    **Step 1:** Filter Performance Values.
    ```
        grep -rin "Performance"  ./case*/out*
    ```
    Then, you can get information related to performance. And performance values can be filtered from that.
    
     **Step 2:** Record and Process Performance Values in Excel.
        Excel can be used to record these values and average them. 
    
     **Step 3:** Draw Tables and Figures in Paper.
        In our paper, Figures 9, 10, and 11 are made by Excel. Figures 15, 16, 17, and 18 are made by Latex.

* The scripts corresponding to the tables and pictures are as follows
    * Table 5
        ```
        cd ../Cases/case_5016/
        sh run-Table5.sh 
        ```

    * Table 6
        ```
        cd ../Cases/case_5016/
        sh run-Table6.sh 
        ```
    * Figure 9
        ```
        cd ../Cases/case_5016/
        sh run-Figure9.sh 
        ```
    * Figure 10 
        ```
        cd ../Cases/case_5016/
        sh run-Figure10.sh 
        ```
    * Figure 11

        ```
        cd ../Cases/case_6688w/
        sh run-Figure11(1).sh 
        cd ../Cases/case_3344w/
        sh run-Figure11(2).sh 
        cd ../Cases/case_1672w/
        sh run-Figure11(3).sh 
        ```
    * Figure 12
        ```
        cd ../Cases/case_5016/
        sh run-Figure12.sh 
        ```
    * Figure 13 and Figure 14
        ```
        cd ../Cases/case_5016/
        ```
        The visualizations of results are obtained by importing file, "_fric.lammpstrj_" into the Ovito software.
        The file, fric.lammpstrj is produced by the command dump in each input file,  "_case_5016/in*_".

    * Figure 15
        ``` 
        cd ../Cases/case_1672w/
        sh run-Figure15(1).sh 
        cd ../Cases/case_5016/
        sh run-Figure15(2).sh 
        ```
    * Figure 16
         ```
        cd ../Cases/case_6688w/
        sh run-Figure16(1).sh 
        cd ../Cases/case_13376w/
        sh run-Figure16(2).sh 
        ```
    * Figure 17
        ``` 
        cd ../Cases/case_13376w/
        sh run-Figure17(1)(2).sh 
        cd ../Cases/case_5016/
        sh run-Figure17(3).sh 
        ```
    * Figure 18
         ```
        cd ../Cases/case_6688w/
        sh run-Figure18(1).sh 
        cd ../Cases/case_13376w/
        sh run-Figure18(2).sh 
        ```
