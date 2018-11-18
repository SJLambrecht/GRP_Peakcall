# GRP_Peakcall
Python script to define peaks from RNA-seq data using a .grp file as input

This script defines peaks/transcriptional units from RNA-seq data that were processed to the .grp file format. Furthermore, it calculates the fold change between a sample and a control for each detected peak. The output is stored as a new file in the .gff file format with the fold change in the 'score' column.

Executing the script:
Before the script is executed there are some changes that have to made to ensure proper functionality.

Most important:
- enter the file names or paths for the grp files for forward reads (#5) and reverse reads (#6)
- enter the genome size of your organism/dataset (#18)

The script will run now, but
- to calculate fold changes you need to specify in wich column are the data for the control (#4) and in which column are the data for the sample (#3). Note, the first column/data set has the index 0.
- you might want to specify based on wich data the peaks should be detected (#2). This can be the same as (#3) or (#4).
- if you run multiple analysis you might want to name your data sets (#1), this will be included in the output file name and prevents that previous data are overwritten and makes result files more easy to identify.
- if you plan to visualize your results using a suitable genome browser you can change the color of the peaks (#17) for each time you run the script to better disciminate between versions of your peak calling.
- if you need to you can assign a frame to your peaks (#16).

#7-#14 are multiple parameters that can be changed in case the results are not as you like them:

#7 is the base range over which is averaged to detect an increase in coverage. It might be good to increase the base range for low coverage data sets or if peaks are increasing only slowly.

#8, #9 is the threshold that has to be passed by the quotient between the two ranges defined by #7

#10, #11 if a start/end is found #10/#11 give the range for which the algorithm tests if there is a greater fold change between the two ranges defined by #7 close by.

#12, #13 defines the space between the two by #7 defined ranges while detecting starts (#12) or ends (#13). It might be good to increase the spacers if peaks are increasing only slowly.

#14 defines in percent how much coverage a peak has to have realtive to the average coverage per nucleotide
