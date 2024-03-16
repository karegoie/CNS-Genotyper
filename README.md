# CNS-Genotyper

Cas9 Mutant Genotyper with Next Generation Sequencing!

## Dependency:

- Bio = From reading fastq files to making alignments

- click = Supporting Command Line Interface

- xlsxwriter = Writing xlsx-shaped log file

## TODO Stack

- Build Pseudocode for using 'Hashmap' to shorten the calc time(100%).

- Rename the classes.

- gRNA-PAM match: Set the Position at Reference class!!!

- Make the code work...!

- Benchmarking: CRISPResso2 vs CNS-Genotyper 1.5

- If memory issue exists: Try using SQLite for Memoization!

## Data folder of CNS-Genotyper:

To run the program, you must put the NGS data in here:

only file with name '.fastq' and '.fastq.gz' will be read by the program,

and other files like this will not be used as a data.

### CNS-Genotyper.exe -x R2 -x undetermined 

To ignore some files after getting all data from the NGS raw result,

you can delete some of it manually, or using the config in the main.py code.

the config [-x R2] will ignore files with name 'R2', which mainly means Read 2.

