1. Introduction
This repository is to make accesible the code to do across chromosome phasing, comile it, and run it. 

2. Compilation
To compile the program a simple make command should do it although a executable file is also provided

3. Commands
The command to run the program is:
/pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing -NbIndiv 16942 -PathInput /pathofthehapfile/filesCHR -PathOutput /pathoftheouput/filesCHR -PathMAF /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt

The command to run the script
sbatch --ntasks=12 --mem=200000M  /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing.sh 16942 /pl/active/KellerLab/Emmanuel/Britishandirish/pedfile123chr /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt

3. Inputs
For the ease of usage a minimum number of parameters are needed. Those parameter are:
-NbIndiv 16942 is the number of individuals in the dataset
The number of individuals has to be for memory requirement lower than 1,000,000.

-PathInput /pathofthehapfile/filesCHR indicates the program where to find the hap files of the 22 chromosomes. /pathofthehapfile/filesCHR1.hap should be the name of the file corresponding to chromosome 1. The hap files 

-PathOutput /pathoftheouput/filesCHR indictae where to write the ped files after across chromosome phasing

-PathMAF /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt is the path of the allele frequency
