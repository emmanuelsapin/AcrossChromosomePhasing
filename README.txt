This repository is to make accesible the code to do across chromosome phasing. For the ease of usage a minimum number of parameters are available. Those parameter are:
-NbIndiv 16942 is the number of individuals in the dataset
-PathInput /pathofthehapfile/filesCHR indicates the program where to find the hap of of the 22 chromosomes. /pathofthehapfile/filesCHR1.hap should be the name of the file corresponding to chromosome 1
-PathOutput /pathoftheouput/filesCHR indictae where to write the hap file after across chromosome phasing
-PathMAF /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt is the path of the allele frequency

The command top run the program is:

/pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing -NbIndiv 16942 -PathInput /pathofthehapfile/filesCHR -PathOutput /pathoftheouput/filesCHR -PathMAF /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt

The command to run the script
sbatch --ntasks=12 --mem=120M  /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing.sh 16942 /pl/active/KellerLab/Emmanuel/Britishandirish/pedfile123chr /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt


					

