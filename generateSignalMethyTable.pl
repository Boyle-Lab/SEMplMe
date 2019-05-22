test
# Useage: perl generateSignalMethyTable.pl name_of_file_in_results name_of_cell_line

#use strict;
#use warnings;
use File::Glob;
my $input = $ARGV[0];
my $cell_type = $ARGV[1];
my $OutputFolder = "." . "/" . "results/" . $input . "/final/";

##remove all past methylation files!!! (.sem, /METHYL, ect.)
#    system("rm","./results/$input/final/ALIGNMENT/*.me");

#unless (-d "../results/$input/METHYL"){ #alt   
#    glob(`mkdir ../results/$input/METHYL`); #alt        
#}

#1. Add methylation data to signal file
print STDERR "Integrating methylation";
my $done_me = &methyl_match();
print STDERR "Done\n";  

#2. Create a SEM
print STDERR "Creating SEM";
my @SEM = &generateSEM($input);
print STDERR "Done\n";

#3. Create a R plot
print STDERR "Creating R plot";
system ("../src/generateRmeplot.pl -TF_name $input -output $OutputFolder");
print STDERR "Done\n";

#All subs found below
sub methyl_match() {
    my $nuc_pos;
    my $nuc;
    my $pos;
    my $loci;
    
    my %methyl_locus_line;
    my $me_locus = "";
#    my $CpG = glob ("./examples/test2");
#    my $CpG = glob ("./examples/WGBS/ENCFF073DUG.wig"); # make general, add as argument
    my $CpG = glob ("./examples/WGBS/$cell_type");    
    open (methyl_file, $CpG) || die "$!\n";
    while (<methyl_file>){
	my @me_line = split(/\s+/, $_);
	unless ($me_line[3] == 0){
	    $me_locus = $me_line[0].':'.$me_line[2]; 
	    $methyl_locus_line{$me_locus} = $me_line[3];
	}
    }	
    
    my @sigfile = glob ("./results/$input/final/ALIGNMENT/*_pos*.signal");
    foreach my $m (0..$#sigfile){
	print ".";
	open (SIG_FILE, $sigfile[$m]) || die "$!\n";
	open (OUT, ">>$sigfile[$m].me") || die "$!\n";
	for ($sigfile[$m] =~ /([ACTG])_pos(\d+).signal/){
	    open (OUT_M, ">>./results/$input/final/ALIGNMENT/M_pos$2.signal.me") || die "$!\n";
	    open (OUT_W, ">>./results/$input/final/ALIGNMENT/W_pos$2.signal.me") || die "$!\n";
	}
	
	while(<SIG_FILE>){
	    chomp($_);
	    for ($sigfile[$m] =~ /([ACTG])_pos(\d+).signal/){
		$nuc_pos = $&;
		$nuc = $1;
		$pos = $2;
	    }
	    
	    my @line = split( /\s/, $_);
	    my $chr = $line[0];
	    my $kmer_start = $line[1];
	    my $kmer_end = $line[2];
	    my $kmer_seq = $line[3];
	    my $strand = $line[4];
	    my $signal = $line[5];
	    my $loci = "";
	    my $me_avg = $inv_me_avg = 0;
	    
	    if ($strand eq '+'){
		$loci = $kmer_start+$pos+1;
		$loci = $chr.':'.$loci;
	    }
	    elsif ($strand eq '-'){
		$loci = $kmer_end-$pos;
		$loci = $chr.':'.$loci;
	    }
	    
	    unless ($signal eq "NA" || $signal eq '-256.000000'){
		if ($nuc eq 'C'){
		    if (exists($methyl_locus_line{$loci})){
			$me_avg = ($methyl_locus_line{$loci}/100)*$signal;
			$inv_me_avg = $signal-$me_avg;
			unless ($methyl_locus_line{$loci} eq '0'){
			    print OUT_M "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\n";
			}
			unless($methyl_locus_line{$loci} eq '100'){
			    print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$inv_me_avg\n";
			}
		    }
		    else {
			print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\n";
		    }
		}
		elsif($nuc eq 'G'){
		    if (exists($methyl_locus_line{$loci})){
			$me_avg = ($methyl_locus_line{$loci}/100)*$signal;
			$inv_me_avg = $signal-$me_avg;
			unless ($methyl_locus_line{$loci} eq '0'){
			    print OUT_W "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\n";
			}
			unless($methyl_locus_line{$loci} eq '100'){
			    print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$inv_me_avg\n";
			}   
		    }
		    else {
			print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\n";
		    }	
		}
		else {
		    print OUT "$_\n";
		}
	    }
	}
	close SIG_FILE;	
	close methyl_file;
    }
}


sub generateSEM() {
    my $enum;
    my $scram;
    my %SNPEffect;
    my $enumerr;
    my $scramerr;
    my %STDErr;
    my $max = 0;
    my %stderr;
    my $nt;
    my $nt_m;
    my $nt_pos_nm;
    my $nt_pos_m;
    my $pos;
    my $total_nm;
    my $total_m;
    my $non_me_signal;
    my $me_signal;
    my $all_signal;
    my @stderr_calc;
#    my @stderr_me_calc;
    my $sqtotal;
    my %nt;
    my %nt_m;
    my %tot;
    my %tot_m;
    my $std;
    my $sterr;
    my %stderr_m;
   
    open (OUT, ">>./results/$input/final/BASELINE/$input.me.sem") || die "$!\n";
    open (TOT, ">>./results/$input/final/BASELINE/$input.me.totals") || die "$!\n";
    open (ERR, ">>./results/$input/final/BASELINE/$input.me.sterr") || die "$!\n";
    
    open (BASE, "./results/$input/final/BASELINE/baseline.maximums") || die "$!\n";
    while (<BASE>){
	my $line = $_;
        chomp ($line);
        my @fields = split ('\t', $line);
        if ($line =~ /Enumerated/) {
	    $enum = $fields[1];
	    $enumerr = $fields[4];
        }
        elsif ($line =~ /Scrambled/) {
	    $scram = $fields[1];
	    $scramerr = $fields[4];
	}
    }
    close (BASE);

#    my @me_files = glob ("./results/$input/final/ALIGNMENT/M_pos*.signal.me"); #test subset   
    my @me_files = glob ("./results/$input/final/ALIGNMENT/*_pos*.signal.me");
    my $length = (($#me_files+1)/6);
    foreach my $f (0..$#me_files){
        open (ME_FILE, $me_files[$f]) || die "$!\n";
	print ".";
        if ($me_files[$f] =~ /([ACTGMW])_pos(\d+).signal.me/){
	    $nt_m = $1;
            #$nt =$nt_m = $1;
            $pos = $2;
	    $nt_pos_m = $nt_m.'_'.$pos;
#	    print "$nt_pos_m\n";
#            $nt_pos_nm = $nt_pos_m = $nt.'_'.$pos;
	}
	$me_signal = $total_m = $all_signal = 0;
#	$non_me_signal = $me_signal = $total_nm = $total_m = 0;

	@stderr_me_calc = (); ####
        while (<ME_FILE>){
            my @line = split(/\s+/, $_);
	    chomp($line[5]);
	    
	    $total_m++;
	    $me_signal = $line[5];
	    $all_signal = $all_signal+$me_signal; ####   
	    push @stderr_me_calc, $me_signal;
	    
        }
	
	my $avg_me_signal;
	
	unless ($total_m == 0){
            $avg_me_signal = $all_signal/$total_m; ####
	    my $sqtotal = 0;
	    foreach my $val (@stderr_me_calc) {
		$sqtotal += ($avg_me_signal - $val) ** 2;
	    }
            my $std_m = ($sqtotal/(scalar(@stderr_me_calc) - 1)) ** 0.5;
	    my $sterr_m = $std_m/sqrt($total_m);
	    
	    unless( $avg_me_signal == 0){ ##added
		$nt_m{$nt_pos_m} = sprintf("%.4f", log2($avg_me_signal / $enum));
		$stderr_m{$nt_pos_m} = sprintf("%.4f", $sterr_m / $enum);
	    }	    
	}
	else { 
	    $nt_m{$nt_pos_m} = 0;
	    $stderr_m{$nt_pos_m} = 0; 
	}
	
        $tot_m{$nt_pos_m} = $total_m;
	close (ME_FILE);
    }
    print OUT "$input\tA\tC\tG\tT\tM\tW\n";
    print TOT "$input\tA\tC\tG\tT\tM\tW\n";
    print ERR "$input\tA\tC\tG\tT\tM\tW\n";
    foreach my $n (0..$length-1){
        my $loc = ($n+1);
        my $A = 'A'.'_'.$n;
        my $C = 'C'.'_'.$n;
        my $G = 'G'.'_'.$n;
        my $T = 'T'.'_'.$n;
	my $M = 'M'.'_'.$n;
	my $W =	'W'.'_'.$n;
        print OUT "$loc\t$nt_m{$A}\t$nt_m{$C}\t$nt_m{$G}\t$nt_m{$T}\t$nt_m{$M}\t$nt_m{$W}\n";
        print TOT "$loc\t$tot_m{$A}\t$tot_m{$C}\t$tot_m{$G}\t$tot_m{$T}\t$tot_m{$M}\t$tot_m{$W}\n";
        print ERR "$loc\t$stderr_m{$A}\t$stderr_m{$C}\t$stderr_m{$G}\t$stderr_m{$T}\t$stderr_m{$M}\t$stderr_m{$W}\n";
    }
    
    $scramerr = $scramerr / $enum;
    $enumerr = $enumerr / $enum;
    return($enumerr, $scramerr);
}


sub log2 {
    my $n = shift;
    return log ($n)/log(2);
}
