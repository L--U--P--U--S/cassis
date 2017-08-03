#!/usr/bin/perl

###############################################################################
# Convert MEME output (in HTML format) to SiTaR input files (in FASTA format) #
# if the predicted motif is also present in the backbone gene                 #
###############################################################################

use strict;
use warnings;

# set organism (folder),
# plus-minus variations (folders)
# and names of backbone genes (folders)
my $organism   = "flavus";
my @plus_minus = qw( +0_-10 +0_-5 +10_-0 +4_-0 +5_-0 +8_-0 alle +0_-4 +0_-8  +10_-10 +4_-4 +5_-5 +8_-8 );
#my @plus_minus = qw(anr);
my @secondary_metabolites = qw(Aflatoxin);

print "These motifs also have occurrences in the backbone promoter:\n";

for my $sm (@secondary_metabolites)
{
	for my $plus_minus_dir (@plus_minus)
	{
		# open MEME file in particular directory
		# --> [organism]/[backbone gene]/[plus-minus]/meme.html
		open( my $meme_file, "<", $organism . "/meme/" . $sm . "/$plus_minus_dir" . "/meme.html" ) or die $!;
		my $motif = 0;
		my $print = 0;
		my @sequences;
		my $locus;

		while (<$meme_file>)
		{
			my $line = $_;

			# section with motif sequences starts
			if ( $line =~ m/<!--data for motif (\d)-->/ )
			{
				$motif     = $1;    # name of motif (mostly "motif 1", "motif 2" or "motif 3")
				@sequences = ();    # --> reset sequences
			}

			# line with locus name and sequence
			if ( $line =~ m/.+_bp__contig_.+  \( \d+\) ([ACGT]+).*/ )
			{
				push( @sequences, $1 );    # gather sequences
			}

			# line with BACKBONE locus name and sequence
			if ( $line =~ m/(.+)__\d+_bp__contig_\d+__BACKBONE  \( \d+\) [ACGT]+.*/ )
			{
				$print = 1;                # motif in backbone sequence --> print this motif data
				$locus = $1;               # name of backbone locus
			}

			# end of motif data section (\>)
			# found some motif data ($motif != 0)
			# motif in backbone gene ($print != 0)
			# gathered some sequences (@sequences != 0)
			# --> write motif sequences to FASTA file
			if ( index( $line, "\">" ) == 0 && $print != 0 && $motif != 0 && @sequences != 0 )
			{
				print "$plus_minus_dir $motif\n";
				
				open( my $tf_file, ">",
					$organism . "/meme/" . $sm . "/$plus_minus_dir" . "/$sm" . "_" . $locus . "_" . $plus_minus_dir . "_Motif" . $motif . ".fa" )
				  or die $!;
				print $tf_file ">$plus_minus_dir" . "_" . "$motif\n";
				print $tf_file join( "\n", @sequences );
				close $tf_file;

				$print = 0;
			}
		}
		close $meme_file;
	}
}
