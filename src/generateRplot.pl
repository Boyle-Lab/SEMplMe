#!/usr/bin/perl
# Natalie Ng, June 8, 2013
# This program is to generate an R input file that will plot out the score.
# The R input file will also be run in the background in R, and a pdf
# with the plot will be outputted.

use File::Basename;
use File::Copy;

use strict;
my $script = basename ($0);
my $TFname;
my $output;
my %StdErrHash;

if ($#ARGV < 0) {&usage;}

while ($_ = shift (@ARGV)) {
        if ($_ =~ /^-TF_name/) {
                $TFname = shift(@ARGV);
        }
        if ($_ =~ /^-output/) {
                $output = shift(@ARGV);
        }

}

# Step 1: Generate the R input file
&generate_input;

# Step 2: Run the R script in the background and export a pdf of the image.
&run_R;

exit 0;

#===============================submethods===================
sub usage {
        my $script = basename ($0);
        print <<"End";

        Usage:
                $script -TF_name <TF_name>

        Output:
                An R plot

End
        exit 0;

}

sub generate_input {
	my $Rout = $output . "/generateRinput.input"; 
	open (OUT_HANDLE, "> $Rout") || die;
	print OUT_HANDLE <<"End";
	pdf("$output/$TFname\_semplot.me.pdf")
	source("src/plotSEMme_Functions.R")
	plotSEM("$output", "$TFname", error=TRUE)
	dev.off ()
End
	close (OUT_HANDLE);
}

sub run_R {
	`R --vanilla < $output/generateRinput.input`;
	`rm $output/generateRinput.input`;
}
