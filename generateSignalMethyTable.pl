
# Useage: perl generateSignalMethyTable.pl name_of_file_in_results name_of_cell_line

use strict;
use warnings;

use File::Glob;
my $input = $ARGV[0];
my $cell_type = $ARGV[1];
my $OutputFolder = "./results/" . $input . "/final/"; #7/30

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
system ("Rscript ./src/plotSEMme_Functions.R  -TF_name $input -output $OutputFolder");
print STDERR "Done\n";


#All subs found below
sub methyl_match() {

#### Calculate methylation for alignment files
    my $nuc_pos;
    my $nuc;
    my $pos;
    my $loci;
    my $me_locus;
    
    my %methyl_locus_line;
#=pod
    my $CpG = glob ("./examples/WGBS/$cell_type");    
    open (methyl_file, $CpG) || die "$!\n";
    while (<methyl_file>){
	my @me_line = split(/\t/, $_);
	unless($me_line[0]=~/\#/){
	    chomp $me_line[3];
	    unless ($me_line[3] == 0){
		$me_locus = $me_line[0].':'.$me_line[2]; 
		$methyl_locus_line{$me_locus} = $me_line[3];
	    }
	}
    }
    
    my @sigfile = glob ("./results/$input/final/ALIGNMENT/*_pos*.signal");
    foreach my $m (0..$#sigfile){
	open (SIG_FILE, $sigfile[$m]) || die "$!\n";
	open (OUT, ">>$sigfile[$m].me") || die "$!\n";
	if ($sigfile[$m] =~ /([CG])_pos(\d+).signal/){
	    print ".";
	    if ($1 eq 'C') {
		open (OUT_M, ">>./results/$input/final/ALIGNMENT/M_pos$2.signal.me") || die "$!\n";
	    } 
	    if ($1 eq 'G') {
		open (OUT_W, ">>./results/$input/final/ALIGNMENT/W_pos$2.signal.me") || die "$!\n";
	    }
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
	    my $me_avg = 0;
	    my $inv_me_avg = 0;
	    my $me_prop = 0;
	    my $inv_me_prop = 0;
	    
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
			$me_prop = $methyl_locus_line{$loci}/100; #9/3                                              
			$inv_me_prop = 1 - $methyl_locus_line{$loci}/100; #9/3  
			unless ($methyl_locus_line{$loci} eq '0'){
			    print OUT_M "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\t$me_prop\n"; #9/3
			}
			unless($methyl_locus_line{$loci} eq '100'){
			    print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$inv_me_avg\t$inv_me_prop\n"; #9/3
			}
		    }
		    else {
			print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\t1\n"; #9/3
		    }
		}
		elsif($nuc eq 'G'){
		    if (exists($methyl_locus_line{$loci})){
			$me_avg = ($methyl_locus_line{$loci}/100)*$signal;
			$inv_me_avg = $signal-$me_avg;
			$me_prop = $methyl_locus_line{$loci}/100; #9/3                                              
			$inv_me_prop = 1 - $methyl_locus_line{$loci}/100; #9/3  
			unless ($methyl_locus_line{$loci} eq '0'){
			    print OUT_W "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\t$me_prop\n"; #9/3
			}
			unless($methyl_locus_line{$loci} eq '100'){
			    print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$inv_me_avg\t$inv_me_prop\n"; #9/3
			}   
		    }
		    else {
			print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\t1\n"; #9/3
		    }	
		}
		else {
		    print OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\t1\n"; #9/3

		}
	    }
	}
	close SIG_FILE;	
	close methyl_file;
	close OUT;
	close OUT_M;
	close OUT_W;
    }
#=cut

    ### Calculate methylation for background files
    my @background = glob ("./results/$input/final/BASELINE/*_kmer.signal");
    my %background_avgs;
   
    foreach my $m (0..$#background){
	open (BACK_FILE, $background[$m]) || die "$!\n";
	open (BACK_OUT, ">>$background[$m].me") || die "$!\n";
	print ".";
	
	my $kount = 0;
	my $all_signal = 0;
	
	while(<BACK_FILE>){
	    chomp($_);
	    
	    my @line = split( /\t/, $_);
	    my $chr = $line[0];
	    my $kmer_start = $line[1];
	    my $kmer_end = $line[2];
	    my $kmer_seq = $line[3];
	    my $strand = $line[4];
	    my $signal = $line[5];
	    my $loci = "";
	    my $me_avg = 0;
	    my $inv_me_avg = 0;
	    my $me_prop = 0;
	    my $inv_me_prop = 0;
	    my $length = length($kmer_seq);
	    my @seq_split = split('',$kmer_seq);
	    
	    foreach my $l (0..$length-1){
		my $nuc = $seq_split[$l];
		if ($strand eq '+'){
		    $loci = $kmer_start+$l+1;
		    $loci = $chr.':'.$loci;
		}
		
		elsif ($strand eq '-'){
		    $loci = $kmer_end-$l;
		    $loci = $chr.':'.$loci;
		}
		
		unless ($signal eq "NA" || $signal eq '-256.000000'){
		    $kount++;
		    $all_signal += $signal;
		    
		    if ($nuc eq 'C'){
			if (exists($methyl_locus_line{$loci})){
			    $me_avg = ($methyl_locus_line{$loci}/100)*$signal;
			    $inv_me_avg = $signal-$me_avg;
			    $me_prop = $methyl_locus_line{$loci}/100; #9/2
			    $inv_me_prop = 1 - $methyl_locus_line{$loci}/100; #9/2  
			    unless ($methyl_locus_line{$loci} eq '0'){
				print BACK_OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\t$me_prop\n"; #9/2
			    }
			    unless($methyl_locus_line{$loci} eq '100'){
				print BACK_OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$inv_me_avg\t$inv_me_prop\n"; #9/2
			    }
			}
			else {
			    if ($l == $length-1){
				print BACK_OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\t1\n"; #9/2
			    }
			}
		    }
		    elsif($nuc eq 'G'){
			if (exists($methyl_locus_line{$loci})){
			    $me_avg = ($methyl_locus_line{$loci}/100)*$signal;
			    $inv_me_avg = $signal-$me_avg;
			    unless ($methyl_locus_line{$loci} eq '0'){
				print BACK_OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$me_avg\t$me_prop\n"; #9/2
			    }
			    unless($methyl_locus_line{$loci} eq '100'){
				print BACK_OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$inv_me_avg\t$inv_me_prop\n"; #9/2
			    }   
			}
			else {
			    if ($l == $length-1){
				print BACK_OUT "$chr\t$kmer_start\t$kmer_end\t$kmer_seq\t$strand\t$signal\t1\n"; #9/2
			    }
			}
		    }
		    else {
			if ($l == $length-1){
			    print BACK_OUT "$_\t1\n"; 
			}
		    }
		}
	    }
	}	    
	$background[$m] =~ /.+BASELINE\/(.+)/;
	my $file = $1;
	chomp($file);	
	$background_avgs{$file} = $all_signal/$kount;
	
	close BACK_FILE;
	close BACK_OUT;
    }
#=cut    
    
####Make baseline.me
    my @back_files = glob ("./results/$input/final/BASELINE/*signal.me");                         
    open (BASELINE, ">","./results/$input/final/BASELINE/baseline.maximums.me") || die "$!\n";

    foreach my $m (0..$#back_files){                                                                                              
	open (BACK_FILE_ME, $back_files[$m]) || die "$!\n";   
	my $kount = 0;
	my $square = 0;
	my $sum = 0;
	my $stddev = 0;
	my $stderr = 0;
	my $avg = 0;
	
	$back_files[$m] =~ /.+BASELINE\/(.+).me/;
	my $file = $1;
	$avg = $background_avgs{$file};
	while(<BACK_FILE_ME>){
	    my @line = split('\t',$_);
	    chomp($line[5]);
	    $kount++;
	    $square = ($line[5]-$avg)**2;
	    $sum += $square;
	}
	$stddev = sqrt((1/($kount-1))*$sum);
	$stderr = ($stddev/sqrt($kount));
	
	print BASELINE "$file\t$avg\t$kount\t$stddev\t$stderr\n";
	close BACK_FILE_ME;
    }
    close BASELINE;
}    

#### Generate SEM from methylation alignment files

sub generateSEM() {
    my $enum;
    my $scram;
    my $enumerr;
    my $scramerr;
    my $nt_m;
    my $nt_pos_m;
    my $pos;
    my $total_m;
    my $me_signal;
    my $all_signal;
    my %nt_m;
    my %tot_m;
    my %stderr_m;
   
    open (OUT, ">>./results/$input/final/$input.me.sem") || die "$!\n";
    open (TOT, ">>./results/$input/final/$input.me.totals") || die "$!\n";
    open (ERR, ">>./results/$input/final/$input.me.sterr") || die "$!\n";
    
    open (BASE, "./results/$input/final/BASELINE/baseline.maximums.me") || die "$!\n";
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
	print ".";
	if ($me_files[$f] =~ /([ACTGMW])_pos(\d+).signal.me/){
	    $nt_m = $1;
            $pos = $2;
	    $nt_pos_m = $nt_m.'_'.$pos;
	}

	$me_signal = $total_m = $all_signal = 0;
	my @stderr_me_calc = ();
	
        while (<ME_FILE>){
            my @line = split(/\s+/, $_);
	    chomp($line[5]);
	    
	    unless ($line[5] < 0){
		$total_m = $line[6] + $total_m; #9/2 $total_m++;
		$me_signal = $line[5];
		$all_signal = $all_signal+$me_signal;
		push @stderr_me_calc, $me_signal;
	    }
	}
	
	my  $avg_me_signal = 0;
	my $sqtotal = 0;
	
	unless ($total_m == 0){
            $avg_me_signal = $all_signal/$total_m;
	    foreach my $val (@stderr_me_calc) {
		$sqtotal += ($avg_me_signal - $val) ** 2;
	    }
            my $std_m = ($sqtotal/(scalar(@stderr_me_calc) - 1)) ** 0.5;
	    my $sterr_m = $std_m/sqrt($total_m);
	    
	    unless($avg_me_signal < 0){
		$nt_m{$nt_pos_m} = sprintf("%.4f", log2($avg_me_signal/$enum));
		$stderr_m{$nt_pos_m} = sprintf("%.4f", $sterr_m/$enum);
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
