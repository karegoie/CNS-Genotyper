
# CNS-Genotyper

Cas9 Mutant Genotyper with Next Generation Sequencing!

## Dependency:

- Bio = From reading fastq files to making alignments

- click = Supporting Command Line Interface

- xlsxwriter = Writing xlsx-shaped log file

## how to use:

Windows: 

1. Download the most recent release(exe file).
2. Move the exe file to the workspace folder you will use.
3. Double-click the program.
4. (Please ignore the warning. As you can see, there is no virus in the software!)
5. The program will automatically generate the needing-files.
6. Edit the src/guide_RNA_seq.txt and src/reference_seq_set.txt files with your sequences.
7. Put all the fastq or fastq.gz files in the data/ folder.
8. Double-click the program again.

If you want to change the variables of program,

use it on a command-line interface.

.\CNS-Genotyper.exe --help will work well.



## TODO Stack

- Build Pseudocode for using 'Hashmap' to shorten the calc time(100%).

- Rename the classes(100%).

- gRNA-PAM match: Set the Position at Reference class!!! (100%)

- Make the code work...! (100%)

- Benchmarking: CRISPResso2 vs CNS-Genotyper 1.5: 1.6 times faster, memory works well.

## Data folder of CNS-Genotyper:

To run the program, you must put the NGS data in here:

only file with name '.fastq' and '.fastq.gz' will be read by the program,

and other files like this will not be used as a data.

### CNS-Genotyper.exe -x R2 -x undetermined 

To ignore some files after getting all data from the NGS raw result,

you can delete some of it manually, or using the config in the main.py code.

the config [-x R2] will ignore files with name 'R2', which mainly means Read 2.

