
#use strict;
#use warnings;

##READ TF FILE
open (TF_LIST, $ARGV[0]) || die "Missing List of TFs";
$finish = 1;

while (<TF_LIST>){
    my @line = split;
    $TF_file = $line[0];
    $TF = $line[1];
    $TF = uc $TF;
    $TF_name = $line[2];
    $cell = $line[3];
    $link = $line[4];
    $count = 0;
    $length = length($finish);

    unless ($cell eq "NA" || $link eq "NA"){
	if ($length > 0){
	    $length = 0;
	    print "PWM\t$TF_name\n";
	    
##DEFINE CELL TYPE FILE
	    print "cell\t$cell\n";
	    
	    if ($cell eq "HepG2") { $cell_type =  "ENCFF073DUG.wig";}
	    elsif ($cell eq "K562") { $cell_type =  "ENCFF872YSC.wig";}
	    elsif ($cell eq "GM12878") { $cell_type =  "ENCFF796NFQ.wig";}
	    elsif ($cell eq "IMR90") { $cell_type =  "ENCFF433WIE.wig";}
	    elsif ($cell eq "H1-hESC") { $cell_type =  "ENCFF770YJW.wig";}
	    elsif ($cell eq "GM23248") { $cell_type =  "ENCFF390POK.wig";}
	    else { next;}
	    

##RUN SEMplMe
	    $finish = system("perl ./generateSignalMethyTable.pl $TF_name $cell_type");
	}
    }
}
