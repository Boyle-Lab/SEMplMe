# SEM_CPP
perl implementation of the SEMplMe algorithm

# Demo

We include a small demo of SEMplMe for HNF4A in HepG2 cells. The expected output is:
```
Integrating methylation...Done
Creating SEM...Done
Creating R plot.............................................................................................................Done
....
```
 
# Usage information
SEMplMe functions on SEMpl output, please run SEMpl before attempting SEMplMe:
https://github.com/Boyle-Lab/SEMpl

SEMplMe requires whole genome bisulfite sequencing (WGBS) data, which can be downloaded from ENCODE. The following example will build the SEM with methylation for HNF4a in HepG2 cells given the example data including an example SEMpl output
```
cd examples

wget "https://www.encodeproject.org/files/ENCFF073DUG/@@download/ENCFF073DUG.bigWig"

cd ..

perl ./generateSignalMethylTable.pl --TF_name HNF4A --WGBS examples/ENCFF073DUG.wig
```
