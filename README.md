# SEMplMe
perl implementation of the SEMplMe algorithm

# Installation of SEMpl
Clone a copy of the SEMpl repository and submodules:

```
git clone --recurse-submodules https://github.com/Boyle-Lab/SEMpl.git
```

Build external libraries:
```
cd SEMpl/lib/libBigWig
make
cd ..
make
mv */*.so .
cd ..
```

Symlink to bowtie index location (use your own index location):
```
ln -s /data/genomes/hg19/bowtie_index/ data
```

Build SEMpl
```
make
```

# Expected SEMplMe output

We include a small demo of SEMplMe for HNF4A in HepG2 cells. The expected output is:
```
Integrating methylation...Done
Creating SEM...Done
Creating R plot.............................................................................................................Done
```
 
# Usage information
SEMplMe functions on SEMpl output, please run SEMpl before attempting to use SEMplMe outside of the demo. Example:

```
./iterativeSEM -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChrom
DnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V04161
01RawRep1.bigWig -TF_name HNF4A -genome data/hg19 -output results/HNF4A
```

For more information on SEMpl please go to  https://github.com/Boyle-Lab/SEMpl


# SEMplMe Demo 
SEMplMe requires whole genome bisulfite sequencing (WGBS) data, which can be downloaded from ENCODE. The following example will build the SEM with methylation for HNF4a in HepG2 cells given the example data including a precomputed example SEMpl output
```
cd examples

wget "https://www.encodeproject.org/files/ENCFF073DUG/@@download/ENCFF073DUG.bigWig"

cd ..

perl ./generateSignalMethylTable.pl --TF_name HNF4A --WGBS examples/ENCFF073DUG.wig
```
