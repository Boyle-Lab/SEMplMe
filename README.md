# SEMplMe
Implementation of the SEMplMe algorithm to integrate SEMpl with methylation data

Sierra S Nishizaki and Alan P Boyle, SEMplMe: A tool for integrating DNA methylation effects in transcription factor binding affinity predictions, BiorXiv 2020.


# System Requirements

## Hardware Requirements
Generation of a SEM requires variable RAM and disk storage based on the size of the initial PWM being considered. For minimal performance, we recommend a computer with the following specs:

RAM: 64+ GB  
CPU: 8+ cores, 3.4+ GHz/core

The runtime on this minimal system is approximately 38 CPU hours. Compile time is approximately 35 seconds.

## Software Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Ubuntu 18.04  
Packages: libcurl4-dev


# Installation
Clone a copy of the SEMplMe repository and submodules. This contains a copy of the SEMpl application that also has to be built:

```
git clone --recurse-submodules https://github.com/Boyle-Lab/SEMplMe.git
```

Build external libraries:
```
cd SEMplMe/lib/libBigWig
make
cd ..
make
mv */*.so .
cd ..
```

Symlink to bowtie index location (use your own index location):
```
ln -s /data/genomes/hg38/bowtie_index/ data
```

Build SEMpl
```
make
```
 
# Demo and Usage information
SEMplMe functions on SEMpl output, please run SEMpl before attempting to use SEMplMe outside of the demo. Example:

```
./iterativeSEM -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChrom
DnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V04161
01RawRep1.bigWig -TF_name HNF4A -genome data/hg19 -output results/HNF4A
```

For more information on SEMpl please go to  https://github.com/Boyle-Lab/SEMpl


## SEMplMe Demo 
SEMplMe requires whole genome bisulfite sequencing (WGBS) data, which can be downloaded from ENCODE. The following example will build the SEM with methylation for HNF4a in HepG2 cells given the example data including a precomputed example SEMpl output
```
cd examples
wget "https://www.encodeproject.org/files/ENCFF073DUG/@@download/ENCFF073DUG.bigWig"
cd ..
perl ./generateSignalMethylTable.pl --TF_name HNF4A --WGBS examples/ENCFF073DUG.wig
```

## Expected SEMplMe output

We include a small demo of SEMplMe for HNF4A in HepG2 cells. The expected output is:
```
Integrating methylation...Done
Creating SEM...Done
Creating R plot.............................................................................................................Done
```
