#!/bin/bash
#SBATCH --job-name="hyperp_DD"
#SBATCH -D .
#SBATCH --output=hyperp_DD.out
#SBATCH --error=hyperp_DD.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=160
#SBATCH --time=08:00:00
##SBATCH --qos=debug
#SBATCH --gres gpu:4


module purge
#module load cuda/9.1   gcc/6.4.0   openmpi/3.0.0 ffmpeg/4.0.2 opencv/3.4.1 cudnn/7.1.3 atlas/3.10.3 scalapack/2.0.2  fftw/3.3.7 szip/2.1.1 rdkit/2018.09.2 python/3.6.5_ML 
module load cuda/9.1  gcc/8.3.0 cuda/10.1 cudnn/7.6.4 nccl/2.4.8 tensorrt/6.0.1 openmpi/4.0.1 fftw/3.3.8 ffmpeg/4.2.1 opencv/4.1.1 atlas/3.10.3 scalapack/2.0.2 szip/2.1.1 python/3.7.4_ML 
python Deep_docking_CYP.py
