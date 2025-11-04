1. Introduction
This repository is to make accesible the code to do across chromosome phasing, comile it, and run it. 

2. Compilation
To compile the program a simple make command should do it although a executable file is also provided

3. Commands
The command to run the program is:
/pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing -NbIndiv 16942 -PathInput /pathofthehapfile/filesCHR -PathOutput /pathoftheouput/filesCHR -PathMAF /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt

The command to run the script
sbatch --ntasks=12 --mem=200000M  /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing.sh 16942 /pl/active/KellerLab/Emmanuel/Britishandirish/pedfile123chr /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt

--ntasks=12 is the number of cpus the program will use. For all loops on the number of individuals are paralelized using the command "#pragma omp parallel for" which use the number of cpus specified in "--ntasks=12". --mem=200000M is the memory need to run the program.
3. Input parameter
For the ease of usage a minimum number of parameters are needed. Those parameter are:
-NbIndiv 16942 is the number of individuals in the dataset
The number of individuals has to be for memory requirement lower than 1,000,000.

-PathInput /pathofthehapfile/filesCHR indicates the program where to find the hap files of the 22 chromosomes. /pathofthehapfile/filesCHR1.hap should be the name of the file corresponding to chromosome 1. Each line of the hap files corresponds to a SNP and should be composed of: chromosome_number RS_number position Allele1 Allele2 then the phenotype as shown in the expemple below:
16 rs79691782 111247 C T 0.0 0.0 0.0 1.0 1.0 1.0
16 rs79562482 111965 A C 1.0 0.0 0.0 0.0 1.0 1.0
For a dataset of two SNPs and three individuals being 
Individual 1: CC CA, 
Individual 2: CT AA, 
and Individual 3: TT AA

-PathOutput /pathoftheouput/filesCHR indicates where to write the output ped files after across chromosome phasing. In those ped file the first haplotype will be detected as coming from the same parent accross the 22 chromosomes.

-PathMAF /pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt is the path of the allele frequency. These allele frequency have to be the frequency of the dataset in the input files.
