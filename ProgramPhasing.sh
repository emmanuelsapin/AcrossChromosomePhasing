#!/bin/sh
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=blanca-ibg  
#SBATCH --mem=1000M


echo Welocome to the script 

NbIndiv=$1
PathInput=$2 
PathOutput=$3
PathMAF=$4

/pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing -NbIndiv "$NbIndiv" -PathInput "$PathInput" -PathOutput "$PathOutput" -PathMAF "$PathMAF" 
