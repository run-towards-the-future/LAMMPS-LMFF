# A Package for LMFF


## A.Abstract
* LAMMPS is one of the most popular Molecular Dynamic (MD) packages and is widely used in the field of physics, chemistry and materials simulation. Layered Materials Force Field (LMFF) is our expansion of the LAMMPS potential function based on the Tersoff potential and  inter-layer potential (ILP) in LAMMPS. LMFF is designed to study layered materials such as graphene and boron hexanitride. It is universal and does not depend on any platform. 

## B. Description
* Artifact name: USER-LMFF
* [![DOI](https://zenodo.org/badge/386895896.svg)](https://zenodo.org/badge/latestdoi/386895896)
* Version: LAMMPS with the stable version of 3Mar2020 from https://download.lammps.org/tars/lammps-3Mar2020.tar.gz
* Key algorithms: Molecular Dynamics
* Compilation: swgcc
* Operating systems and versions: Sunway Raise OS with Linux kernel 3.8.0
* Libraries and versions: MVAPICH2 version 2.2a/SWCH MPI at 20210101
* Execution: By the job system.
* Output: Reported performance in LAMMPS.
* Publicly available?: Yes.

## C.SUMMARY OF THE USER-LMFF PACKAGE
* Recently, we have made some improvements and optimizations on the virial stress part of the calculation of interlayer forces. And, based on the stable version of LAMMPS: 3Mar2020, we installed, compiled and tested it on Sunway TaihuLight.
* In order to facilitate the installation and use, we follow the conventional practice in LAMMPS and use the latest LMFF code as an installation package.
The installation process is as follows:
```
 make yes-USER-MISC
 make yes-MOLECULE
 make yes-MANYBODY 
 make yes-USER-LMFF
```
* **Note:** *make yes-USER-LMFF* needs to be executed at the end.

