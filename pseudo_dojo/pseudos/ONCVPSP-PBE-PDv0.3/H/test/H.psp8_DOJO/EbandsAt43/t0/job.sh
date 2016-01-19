#!/bin/bash

#PBS -q main
#PBS -N w0_t0
#PBS -l select=20:ncpus=1:vmem=2048mb:mpiprocs=1
#PBS -l pvmem=2048mb
#PBS -l walltime=4:0:0
#PBS -W group_list=tcos
#PBS -o /home/acad/ucl-naps/setten/software/pseudo_dojo/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/H/test/H.psp8_DOJO/EbandsAt43/t0/queue.qout
#PBS -e /home/acad/ucl-naps/setten/software/pseudo_dojo/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/H/test/H.psp8_DOJO/EbandsAt43/t0/queue.qerr
#PBS -r y

# Load Modules
module purge
module load compiler/intel/composerxe/2013_sp1.1.106
module load intelmpi
module load python/2.7

# Shell Environment
export PATH=/projects/acad/napsimu/software/abinit/7.11.5/dev/bin/:$PATH

cd /home/acad/ucl-naps/setten/software/pseudo_dojo/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/H/test/H.psp8_DOJO/EbandsAt43/t0
# Commands before execution
source ~/to_dojo
ulimit  # pre_run is a string in verbatim mode (note |)


mpirun -n 20 abinit < /home/acad/ucl-naps/setten/software/pseudo_dojo/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/H/test/H.psp8_DOJO/EbandsAt43/t0/run.files > /home/acad/ucl-naps/setten/software/pseudo_dojo/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/H/test/H.psp8_DOJO/EbandsAt43/t0/run.log 2> /home/acad/ucl-naps/setten/software/pseudo_dojo/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/H/test/H.psp8_DOJO/EbandsAt43/t0/run.err
