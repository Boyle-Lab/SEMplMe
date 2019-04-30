
# Useage: perl generateSignalMethyTable.pl name_of_file_in_results

#use strict;
#use warnings;
use File::Glob;
my $input = $ARGV[0];
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
#my @SEM = &generateSEM($input);
print STDERR "Done\n";

#3. Create a R plot
print STDERR "Creating R plot";
#system ("../src/generateRmeplot.pl -TF_name $input -output $OutputFolder");
print STDERR "Done\n";

#All subs found below

sub methyl_match() {
    my $nuc_pos;
    my $nuc;
    my $pos;
    my $loci;
    
    my %methyl_locus_line;
    my $me_locus = "";
    my $CpG = glob ("./examples/test2"); # make general
    open (methyl_file, $CpG) || die "$!\n";
    while (<methyl_file>){
	my @me_line = split(/\s+/, $_);
	if ($me_line[5]  eq "+"){
	    $me_locus = $me_line[0].':'.$me_line[2]; 
	}
	elsif ($me_line[5]  eq "-"){
            $me_locus = $me_line[0].':'.$me_line[2];                      
	}
	####if hash exists add to list
	$methyl_locus_line{$me_locus} = $_;
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
	    
	    unless ($signal eq "NA" ){ #|| $signal eq '-256.000000'){ ##added -256
		if ($nuc eq 'C'){
		    if (exists($methyl_locus_line{$loci})){
			print "$methyl_locus_line{$loci}\n";
			my @me_split = split('\s', $methyl_locus_line{$loci});
			chomp($me_split[9]);
			$me_avg = ($me_split[10]/100)*$signal;
			$inv_me_avg = $signal-$me_avg;
			unless ($me_split[10] eq '0'){
			    print OUT_M "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\n";
			}
			unless($me_split[10] eq '100'){
			    print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$inv_me_avg\n";
			}
		    }
		    else {
			print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\n";
		    }
		}
		elsif($nuc eq 'G'){
		    if (exists($methyl_locus_line{$loci})){
			my @me_split = split('\s', $methyl_locus_line{$loci});
			chomp($me_split[9]);
			$me_avg = ($me_split[10]/100)*$signal;
			$inv_me_avg = $signal-$me_avg;
			unless ($me_split[10] eq '0'){
			    print OUT_W "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\n";
			}
			unless($me_split[10] eq '100'){
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
    my @stderr_calc;
    my @stderr_me_calc;
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
    
    my @me_files = glob ("./results/$input/final/ALIGNMENT/*_pos*.signal.me");
    my $length = (($#me_files+1)/6);
    foreach my $f (0..$#me_files){
        open (ME_FILE, $me_files[$f]) || die "$!\n";
        print "$me_files[$f]\n";
	print ".";
        if ($me_files[$f] =~ /([ACTG])_pos(\d+).signal.me/){
            $nt =$nt_m = $1;
            $pos = $2;
            $nt_pos_nm = $nt_pos_m = $nt.'_'.$pos;
	}
	$non_me_signal = $me_signal = $total_nm = $total_m = 0;

	##### CHANGE FORMAT OF .ME FILES####
        while (<ME_FILE>){
            my @line = split(/\s+/, $_);
	    if ($line[6] < 1){
                $total_nm++;
                $non_me_signal = $non_me_signal+($line[5]*(1-$line[6]));
                push @stderr_calc, $non_me_signal;
            }
	    if ($line[6] > 0){
                $total_m++;
                $me_signal = $me_signal+($line[5]*$line[6]);
                push @stderr_me_calc, $me_signal;
            }
        }
	my $avg_non_me_signal;
	my $avg_me_signal;
	unless ($total_nm == 0){
            $avg_non_me_signal = $non_me_signal/$total_nm;
	    my $sqtotal = 0;
	    foreach my $val (@stderr_calc) {
		$sqtotal += ($avg_non_me_signal - $val) ** 2;
	    }
	    my $std = ($sqtotal/(scalar(@stderr_calc) - 1)) ** 0.5;
	    my $sterr = $std/sqrt($total_nm);
	    
	    $nt{$nt_pos_nm} = sprintf("%.4f", log2($avg_non_me_signal / $enum));
	    $stderr{$nt_pos_nm} = sprintf("%.4f", $sterr / $enum);
	}
	else {
	    $nt{$nt_pos_nm} = 0; 
	    $stderr{$nt_pos_nm} = 0;
	}
	
	unless ($total_m == 0){
            $avg_me_signal = $me_signal/$total_m;
	    my $sqtotal_m = 0;
	    foreach my $val (@stderr_me_calc) {
                $sqtotal += ($avg_me_signal - $val) ** 2;
	    }
            my $std_m = ($sqtotal/(scalar(@stderr_me_calc) - 1)) ** 0.5;
	    my $sterr_m = $std_m/sqrt($total_m);
	    
	    $nt_m{$nt_pos_m} = sprintf("%.4f", log2($avg_me_signal / $enum));
	    $stderr_m{$nt_pos_m} = sprintf("%.4f", $sterr_m/ $enum);
	}	    
	else { 
	    $nt_m{$nt_pos_m} = 0;
	    $stderr_m{$nt_pos_m} = 0; 
	}
	
        $tot{$nt_pos_nm} = $total_nm;
        $tot_m{$nt_pos_m} = $total_m;
	
	close (ME_FILE);
    }
    print OUT "$input\tA\tC\tG\tT\tCm\tGm\n";
    print TOT "$input\tA\tC\tG\tT\tCm\tGm\n";
    print ERR "$input\tA\tC\tG\tT\tCm\tGm\n";
    foreach my $n (0..$length){
        my $loc = ($n+1);
        my $A = 'A'.'_'.$n;
        my $C = 'C'.'_'.$n;
        my $G = 'G'.'_'.$n;
        my $T = 'T'.'_'.$n;
        print OUT "$loc\t$nt{$A}\t$nt{$C}\t$nt{$G}\t$nt{$T}\t$nt_m{$C}\t$nt_m{$G}\n";
        print TOT "$loc\t$tot{$A}\t$tot{$C}\t$tot{$G}\t$tot{$T}\t$tot_m{$C}\t$tot_m{$G}\n";
        print ERR "$loc\t$stderr{$A}\t$stderr{$C}\t$stderr{$G}\t$stderr{$T}\t$stderr_m{$C}\t$stderr_m{$G}\n";
    }
    
    $scramerr = $scramerr / $enum;
    $enumerr = $enumerr / $enum;
    return($enumerr, $scramerr);
}


sub log2 {
    my $n = shift;
#    return log ($n)/log(2);
}
