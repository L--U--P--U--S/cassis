#!/usr/bin/perl

# TODO option to provide a set of genes as "guidance" for the cluster prediction

##########################################
#        ..::| C A S S I S |::..         #
#                                        #
# Cluster ASSignment by Islands of Sites #
##########################################

# Copyright (C) 2015 Leibniz Institute for Natural Product Research and
# Infection Biology -- Hans-Knoell-Institute (HKI)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Dieses Programm ist Freie Software: Sie können es unter den Bedingungen
# der GNU General Public License, wie von der Free Software Foundation,
# Version 3 der Lizenz oder (nach Ihrer Wahl) jeder neueren
# veröffentlichten Version, weiterverbreiten und/oder modifizieren.
#
# Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber
# OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
# Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
# Siehe die GNU General Public License für weitere Details.
#
# Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
# Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
#
# Contact: gianni.panagiotou@leibniz-hki.de, thomas.wolf@leibniz-hki.de
#
# Cite: If you use CASSIS, please cite https://doi.org/10.1093/bioinformatics/btv713

use 5.14.2;

use strict;
use warnings;

use Getopt::Long;                # handle command line arguments
use File::Path qw(make_path);    # create directories (and SUBdirectories!)
use File::Basename;              # parse file paths into directory, filename and suffix
use File::Copy::Recursive qw(dircopy);
use File::Copy;
use File::Compare;               # compare files --> check if same file already present
use Roman;                       # conversion between roman and arabic numbers
use Bio::SeqIO;                  # handle sequences in fasta format
use List::AllUtils qw(sum max min any none uniq);
use IPC::System::Simple qw(system);    # improved behavior of Perl's "system()" call
use File::Which;                       # emulation of "which" command line tool
use Sys::CpuAffinity;                  # parallelization
use Parallel::ForkManager;             # parallelization
use Text::CSV;                         # read and write csv files
use Text::Autoformat;                  # wrap text automatically (especially helpful for usage text)

# debug
#use Data::Dumper;
#$Data::Dumper::Terse = 1;
#use Devel::Size qw(size total_size);

#########
# usage #
#########

sub usage
{
	my $usage = autoformat(
		"Copyright (C) 2015 Leibniz Institute for Natural Product Research and Infection Biology -- Hans-Knoell-Institute (HKI).\n\n"
		  . "This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions. See the file COPYING, which you should have received along with this program, for details.\n\n"
		  . "      ###########################################################\n"
		  . "      # CASSIS - (C)luster (ASS)ignment by (I)slands of (S)ites #\n"
		  . "      ###########################################################\n\n"
		  . "CASSIS is a tool to precisely predict secondary metabolite (SM) gene clusters around a given anchor (or backbone) gene. Genes encoding SMs tend to be clustered. Gene clusters are small groups of normally up to 20 genes; tightly co-localized, co-regulated, and participating in the same metabolic pathway. The predictions are based on transcription factor binding sites shared by promoter sequences of the putative cluster genes.\n\n"
		  . "\n"
		  . "usage: cassis.pl <parameters>\n\n" . "\n"
		  . "required parameters:\n\n"
		  . "   --annotation, -a <file>\n\n"
		  . "      Genome annotation file. A simple text file in tabular format with at least five columns. These are:\n\n"
		  . "      gene <string> | contig <string> | start position <int> | stop position <int> | strand <+ or->.\n\n"
		  . "      The column separator must be tab (\\t). The annotation file must be sorted ascending by contig, start and stop. Start positions must be smaller than stop positions. Contig names have to coincide with the genome sequence file.\n\n"
		  . "   --genome, -g <file>\n\n"
		  . "      Genomic sequence file. A multiFASTA file containing the DNA sequences of all contigs of the species. Contig names have to coincide with the annotation file.\n\n"
		  . "\n"
		  . "optional parameters:\n\n"
		  . "   --anchor, -b <ID>\n\n"
		  . "      Feature ID of the clusters anchor gene. The ID has to coincide with the annotation file. This gene will be the starting point of the cluster prediction.\n\n"
		  . "   --cluster, -c <name>\n\n"
		  . "      Name of the gene cluster to predict. Creates a sub-directory with the given name to separate predictions for different clusters but same species (in the same working directory). Will be set to the cluster's anchor gene ID, if ommitted.\n\n"
		  . "   --dir, -d <directory>\n\n"
		  . "      Working  directory. All directories and files generated by CASSIS will go in here. Choosing different directories for e.g. different species might be a good idea.\n\n"
		  . "   --fimo, -f [<p-value cut-off>]\n\n"
		  . "      Enables the motif search with FIMO. Setting a maximum p-value, which will restrict the number of binding sites found by FIMO, is optional. If no p-value cut-off is specied, the default value of 0.00006 will be used. To run the motiv search, the MEME suite must be installed on your system. See http://meme-suite.org/doc/download.html\n\n"
		  . "   --frequency, -fq <frequency cut-off between 0 and 100>\n\n"
		  . "      By default, CASSIS will skip the cluster prediction for a certain motif, if this motif results in a number of binding sites found in more than 14 % of all promoters in the genome. By changing this paramater you may set a different frequency cut-off.\n\n"
		  . "   --gap-length, -gl <0-5>\n\n"
		  . "      Sets the maximum number of promoters in a row WITHOUT binding sites, which are allowed inside a cluster prediction. Allowed are 0, 1, 2, 3, 4, and 5. The default value is 2. Only applies if the cluster prediction has been enabled (see parameter --prediction).\n\n"
		  . "   --help, -h\n\n"
		  . "      Show the help page, which you are currently reading.\n\n"
		  . "   --meme, -m [<e-value cut-off>]\n\n"
		  . "      Enables the motif prediction around the anchor gene with MEME. Setting a maximum e-value for found motifs is optional. MEME will stop, if it reaches the e-value, regardless if any motif has been found or not. If no e-value cut-off is specified, the default value of 1.0e+005 will be used. The cut-off -1 has a special meaning: This means NO cut-off and MEME will not stop until it finds at least one motif, regardless of its e-value. To run the motiv prediction, the MEME suite must be installed on your system. See http://meme-suite.org/doc/download.html\n\n"
		  . "   --mismatches, -mm <allowed mismatches>\n\n"
		  . "      Setting a maximum number of allowed mismatches (0-3) per binding site sequence is optional. The default is 0. See --sitar.\n\n"
		  . "   --num-cpus, -n <number>\n\n"
		  . "      Set number of CPUs (or CPU cores) to use. The motif prediction (parameter --meme) and the motif search (parameter --fimo) steps will then run with multiple forks in parallel. CASSIS uses only one CPU by default.\n\n"
		  . "   --prediction, -p\n\n"
		  . "      Enables the cluster prediction, including the motif prediction via paramter --meme and the motif search via parameter --fimo. By default it's disabled.\n\n"
		  . "   --sitar, -s <file>\n\n"
		  . "      Alternative to MEME and FIMO: Enables the motif search with SiTaR. You have to provide a file in (multi) FASTA format with binding site sequences of at least one transcription factor.\n\n"
		  . "   --verbose, -v\n\n"
		  . "      By default, only the most important information is printed. Use this parameter if you prefer a more verbose output.\n\n" . "\n"
		  . "runtime (Intel Xeon, running at 2.7 GHz):\n\n"
		  . "   Using 2 CPUs (via --num-cpus), predicting the cluster to a given anchor gene takes about 40 min. Using 4 CPUs, it takes about 22 min. And using 60 CPUs, about 3 min.",
		{ justify => "full", all => 1 }
	  )
	  . "\n\n"
	  . "command-line example:" . "\n"
	  . "cassis.pl --dir fungi/fumigatus/ --annotation fungi/fumigatus/A_fumigatus_Af293_version_s03-m04-r22_features.csv --genome fungi/fumigatus/A_fumigatus_Af293_version_s03-m04-r22_chromosomes.fasta --anchor Afu6g09660 --cluster Gliotoxin --meme 1.0e+005 --fimo 0.00006 -frequency 14 --prediction --gap-length 2 --num-cpus 2 --verbose\n"
	  . "\n\n"
	  . "Contact: gianni.panagiotou\@leibniz-hki.de, thomas.wolf\@leibniz-hki.de"
	  . "\n\n"
	  . "Cite: If you use CASSIS, please cite https://doi.org/10.1093/bioinformatics/btv713";

	# print README to file
	open( my $readme, ">", "README" ) or die "Cannot write to file \"./README\".\n", $!;
	say $readme $usage;
	close $readme;

	# and show it on the screen
	if ( which("tty") )
	{
		if ( system( "tty", "-s" ) == 0 )    # check if stdin is a terminal
		{
			system( "less", "README" );
		}
		else
		{
			system( "cat", "README" );
		}
	}
	else                                     # no tty --> possibly running on windows
	{
		system("type README");
	}

	exit;
}

##########################
# Command line arguments #
##########################
sub command_line_arguments { }    # the only purpuse of such empty subs is to create an entry in Eclipse' outline view for faster navigation

my $verbose         = 0;          # print a lot of information (1) or only a minimum (0)
my $dir             = "";         # working directory --> organism
my $annotation_file = "";
my $genome_file     = "";

my $cluster_prediction = 0;       # enable cluster prediction yes/no? --> default: no
my $max_gap_length     = 2;       # maxium gap length allowed in cluster predictions --> default: 2
my $meme;                         # predict motifs around anchor gene with MEME (e-value cut-off)?
my @meme_parameters = qw(-dna -nostatus -mod anr -nmotifs 1 -minw 6 -maxw 12 -revcomp);
my $fimo;                         # compute frame scores based on FIMO results (p-value cut-off)?
my $sitar;                        # compute frame scores based on SiTaR results?
my $mismatches;                   # maximum number of allowed mismatches, if applying SiTaR instead of MEME and FIMO
my $frequency_cutoff = 14;        # binding site frequency cut-off --> default: 14 (determinded by parameter estimation of set of known gene clusters)
my $backbone_id;                  # gene ID of the clusters anchor/backbone protein (NRPS, PKS, …)
my $cluster_name;                 # name of the SM gene cluster (aka subfolder)
my $num_cpus = 1;                 # number of CPUs to use --> default: 1 (CPU-) core
my $help;                         # display help for command line arguments?

# no command line arguments specified at all?
usage() if !@ARGV;

# save whole command line
my $commandline = join( " ", ( fileparse($0) )[0], @ARGV );

# read arguments from command line
GetOptions(
	"dir|d=s"         => \$dir,
	"annotation|a=s"  => \$annotation_file,
	"genome|g=s"      => \$genome_file,
	"meme|m:s"        => \$meme,
	"anchor|b=s"      => \$backbone_id,
	"cluster|c=s"     => \$cluster_name,
	"fimo|f:s"        => \$fimo,
	"sitar|s=s"       => \$sitar,
	"mismatches|mm=i" => \$mismatches,
	"frequency|fq=i"  => \$frequency_cutoff,
	"prediction|p"    => \$cluster_prediction,
	"gap-length|gl=i" => \$max_gap_length,
	"num-cpus|n=i"    => \$num_cpus,
	"verbose|v"       => \$verbose,
	"help|h"          => \$help
) or die "[ERROR] Please check your command line arguments.\nStopped";

# or asked for help?
usage() if $help;

say "(0) " . localtime . " - SETTINGS …";

# check ANNOTATION argument
if ( $annotation_file =~ m/^([-\w\.\/\+]+)$/ and -s $1 )
{
	$annotation_file = $1;
	say "Genome annotation file: $annotation_file";
}
elsif ( !$annotation_file )
{
	die "[ERROR] Please specify an ANNOTATION file via the --annotation or -a argument.\nStopped";
}
else
{
	die "[ERROR] in ANNOTATION (--annotation, -a) argument \"$annotation_file\".\nStopped";
}

# check GENOME argument
if ( $genome_file =~ m/^([-\w\.\/\+]+)$/ and -s $1 )
{
	$genome_file = $1;
	say "Genomic sequence file: $genome_file";
}
elsif ( !$genome_file )
{
	die "[ERROR] Please specify a GENOME file via the --genome or -g argument.\nStopped";
}
else
{
	die "[ERROR] in GENOME (--genome, -g) argument \"$genome_file\".\nStopped";
}

# check DIR argument
if ($dir)
{
	if ( $dir =~ m/^([-\w\.\/\+]+)$/ and -d $1 )
	{
		$dir = $1;
		$dir .= "/" if ( substr( $dir, -1, 1 ) ne "/" );
		make_path($dir);
		say "Working directory: $dir";
	}
	else
	{
		die "[ERROR] in DIR (--dir, -d) argument \"$dir\".\nStopped";
	}
}
elsif ( ( fileparse($annotation_file) )[1] eq ( fileparse($genome_file) )[1] )
{
	$dir = ( fileparse($annotation_file) )[1];    # or $genome_file - doesn't matter, because they are the same
	say "Working directory: $dir";
}
else
{
	die
"[ERROR] Either make sure the annotation (--annotation, -a) and sequence file (--genome, -g) are located in the same directory or specify a working directory via --dir.\nStopped";
}

# check BACKBONE and CLUSTER argument
if ($backbone_id)
{
	if ( $backbone_id =~ m/^([-\w\.:]+)$/ )
	{
		$backbone_id = $1;
		say "Anchor gene: $backbone_id";

		if ($cluster_name)
		{
			if ( $cluster_name =~ m/^([-\w\.:]+)$/ )
			{
				$cluster_name = $1;
				say "Cluster name: $cluster_name";
			}
			else
			{
				die "[ERROR] in CLUSTER (--cluster, -c) argument \"$cluster_name\".\nStopped";
			}
		}
		else
		{
			$cluster_name = $backbone_id;
		}
	}
	else
	{
		die "[ERROR] in ANCHOR gene (--anchor, -b) argument \"$backbone_id\".\nStopped";
	}
}

# check PREDICTION argument
if ($cluster_prediction)
{
	$meme = "" if !defined $meme;
	$fimo = "" if !defined $fimo;
	say "Cluster prediction: yes";
}
else
{
	say "Cluster prediction: no";
}

# check MOTIF SEARCH (FIMO/SiTaR) argument
my $fimo_exec  = which("fimo");     # get location of fimo executable, may be undefined
my $sitar_exec = which("sitar");    # get location of sitar executable, may be undefined
if ( defined $fimo or defined $sitar )
{
	if ( defined $sitar )
	{
		die "[ERROR] Please get a copy of SiTaR and add it to your PATH environment.\nStopped" if !$sitar_exec;

		if ( $sitar =~ m/^([-\w\.\/\+]+)$/ and -s $1 )
		{
			# sitar overrides/deactivates meme and fimo settings
			$meme = undef;
			$fimo = undef;

			$sitar = $1;
			$mismatches = "" if !defined $mismatches;
			say "SiTaR binding site file: $sitar";
		}
		elsif ( !$sitar )
		{
			die "[ERROR] Please specify a BINDING SITE file via the --sitar or -s argument.\nStopped";
		}
		else
		{
			die "[ERROR] in BINDING SITE (--sitar, -s) argument \"$sitar\".\nStopped";
		}

		if ( $mismatches eq "" )
		{
			$mismatches = 0;
			say "Binding site search with: SiTaR, up to $mismatches mismatche(s) allowed (default)";
		}
		elsif ( $mismatches =~ m/^(\d)$/ and $1 >= 0 and $1 <= 3 )
		{
			$mismatches = $1;
			say "Binding site search with: SiTaR, up to $mismatches mismatches allowed";
		}
		else
		{
			die "[ERROR] in SiTaR mismatch (--mismatches, -mm) value \"$mismatches\", which should be 0-3.\nStopped";
		}
	}
	elsif ( defined $fimo )
	{
		die "[ERROR] Cannot find FIMO. Please install the MEME suite and add it to your PATH environment.\nStopped" if ( !$fimo_exec );

		$meme = "" if ( !defined $meme );

		if ( $fimo eq "" )
		{
			# default value determinded by parameter estimation of set of known gene clusters
			$fimo = 0.00006;
			say "Binding site search with: FIMO, p-value cut-off " . sprintf( "%e", $fimo ) . " (default)";
		}
		elsif ( $fimo =~ m/^(0\.\d+)\z/ )
		{
			$fimo = $1;
			say "Binding site search with: FIMO, p-value cut-off " . sprintf( "%e", $fimo );
		}
		else
		{
			die "[ERROR] in FIMO p-value cut-off (--fimo, -f) argument \"$fimo\", which should rather look like e.g. \"0.0001\".\nStopped";
		}
	}
}

# check MOTIF PREDICTION (MEME) argument
my $meme_exec = which("meme");    # get location of meme executable, may be undefined
if ( defined $meme )
{
	die "[ERROR] Cannot find MEME. Please install the MEME suite and add it to your PATH environment.\nStopped"           if !$meme_exec;
	die "[ERROR] Please specify an ANCHOR gene via --anchor, if you want to predict motifs around it with MEME.\nStopped" if !$backbone_id;

	if ( $meme eq "" )
	{
		# default value determinded by parameter estimation of set of known gene clusters
		$meme = "1.0e+005";
		say "Motif prediction with: MEME, e-value cut-off 1.0e+005 (default)";
	}
	elsif ( $meme =~ m/(\d\.\de[+|-]\d{3})\z/ )
	{
		$meme = $1;
		say "Motif prediction with: MEME, e-value cut-off $meme";
	}
	else
	{
		die "[ERROR] in MEME e-value cut-off (--meme, -m) argument \"$meme\", which should look like e.g. \"1.2e-003\".\nStopped";
	}

	push( @meme_parameters, "-evt", $meme );
}

# check BINDING SITE FREQUENCY argument
if ( defined $fimo )
{
	if ( $frequency_cutoff =~ m/^(\d{1,2})$/ and $frequency_cutoff > 0 and $frequency_cutoff < 100 )
	{
		if ( $frequency_cutoff == 14 )
		{
			say "Binding site frequency cut-off: $frequency_cutoff % (default)";
		}
		else
		{
			say "Binding site frequency cut-off: $frequency_cutoff %";
		}
		$frequency_cutoff = $frequency_cutoff / 100;
	}
	else
	{
		die
		  "[ERROR] in binding site FREQUENCY cut-off (--frequency, -fq) argument, which must be greater than 0 and lower than 100, f.e. 20.\nStopped";
	}
}

# check GAP LENGTH argument
if ($cluster_prediction)
{
	if ( $max_gap_length =~ m/^(\d)$/ and $max_gap_length >= 0 and $max_gap_length <= 5 )
	{
		if ( $max_gap_length == 2 )
		{
			say "Maximum gap length: $max_gap_length (default)";
		}
		else
		{
			say "Maximum gap length: $max_gap_length";
		}
	}
	else
	{
		die
"[ERROR] in GAP LENGTH (--gap-length, -gl) argument \"$max_gap_length\": Please specify a maximum gap length of either 0, 1, 2, 3, 4, or 5.\nStopped";
	}
}

# check CPU argument (if it's not the default value)
# --> only necessary if we are going to run MEME or FIMO
if ( $meme or $fimo )
{
	if ( $num_cpus != 1 )
	{
		die "[ERROR] in CPU (--num-cpus, -n) argument \"$num_cpus\": please specify a number greater than 0.\nStopped"
		  if ( $num_cpus < 1 );    # die if zero or negative
		$num_cpus = Sys::CpuAffinity::getNumCpus() if ( $num_cpus > Sys::CpuAffinity::getNumCpus() );    # set to installed CPUs if greater
		say "CPUs to use: $num_cpus";
	}
	else
	{
		say "CPUs to use: 1 (default, no parallelization)";
	}
}

# print whole command line
say "\n" . "Command line summary:\n" . $commandline if $verbose;

# additional scipts should be in the same directory like the main script itselfubu
my $self_dir = ( fileparse(__FILE__) )[1];
require $self_dir . "cassis_meme.pl";

# made all checks --> set new and UNTAINTED path
delete @ENV{ "IFS", "CDPATH", "ENV", "BASH_ENV" };    # make %ENV safer
my @minimal_path = qw(/bin /usr/bin);                 # "select" minimal path
my @untainted;

# path to executable may be undefined, if not in original (potentially tainted) $PATH; maybe don't need that one this time
for ( grep { defined } ( $meme_exec, $fimo_exec, $sitar_exec ) )    # $R_exec
{
	# remove executable (keep dir only) and trailing "/" to match usual style in $PATH variables
	if ( $_ =~ m/^(\/[-\w\/]+)\/[-\w\.]+$/ )
	{
		push( @untainted, $1 );
	}
	else
	{
		die "[ERROR] \"$_\" seems to be an insecure path.\nStopped";
	}

}

# prevent multiple occurrences of same path
$ENV{PATH} = join( ":", uniq( @minimal_path, @untainted ) );

#########################
# Additional parameters #
#########################

# assumed start and end of promotors related to TSS (+1)
my $upstream_tss        = 1000;    # nucleotides upstream TSS
my $downstream_tss      = 50;      # nucleotides downstream TSS
my $complete_intergenic = 0;       # use complete intergenic region --> ignore $upstream_tss (but keep $downstream_tss)

#####################
# Genome annotation #
#####################
sub genome_annotation { }

say "\n(1) Processing genome ANNOTATION …";

my @features;                      # {ID, first, last, strand, contig}

my $csv = Text::CSV->new(
	{
		allow_loose_quotes  => 1,      # bad practice in csv format (blame creator of the input file?)
		allow_loose_escapes => 1,      # bad practice in csv format (blame creator of the input file?)
		sep_char            => "\t",
		eol                 => $/,
		binary              => 1
	}
);

# open annotation file
open( my $annotation_io, "<", $annotation_file ) or die "Cannot read from genome annotation file \"$annotation_file\".\n", $!;

# analyse each line of the annotation file
while ( my $row = $csv->getline($annotation_io) )
{
	# jump to next line if current line is:
	# - empty
	# - comment (! or #)
	next if ( "@$row" eq "" or any { index( $row->[0], $_ ) == 0 } ( "!", "#" ) );

	# gene <string> | contig <string> | start position <string> | end position <int> | strand <+ or-> | … ### taken right from usage string

	my $ID            = $row->[0];
	my $contig_number = $row->[1];
	my $first         = $row->[2];
	my $last          = $row->[3];
	my $strand        = $row->[4];    # W = +; C = -

	$strand =~ s/\r$//;               # remove trailing carriage return (if present) --> necessary to correctly handle windows newline encoding "\r\n"

	die "[ERROR] Invalid whitespace in contig string \"$contig_number\".\nStopped" if ( index( $contig_number, " " ) != -1 );
	die "[ERROR] No gene name in line \"@$row\".\nStopped"                         if ( !$ID );
	die "[ERROR] No contig name in line \"@$row\".\nStopped"                       if ( !$contig_number );
	die "[ERROR] No gene start position in line \"@$row\".\nStopped" if ( !$first and $first != 0 );    # 0 is false, but still a valid start pos
	die "[ERROR] No gene stop position in line \"@$row\".\nStopped"  if ( !$last );
	die "[ERROR] No strand information in line \"@$row\".\nStopped"  if ( !$strand );
	die "[ERROR] Start position is not an positive integer in line \"@$row\".\nStopped" if ( $first !~ /^\d+$/ );
	die "[ERROR] Stop position is not an positive integer in line \"@$row\".\nStopped"  if ( $last !~ /^\d+$/ );
	die "[ERROR] Start position identical to stop position in line \"@$row\".\nStopped" if ( $first == $last );
	die "[ERROR] Start position AFTER stop position in line \"@$row\".\nStopped"        if ( $first > $last and $strand eq "+" );

	( $first, $last ) = ( $last, $first ) if ( $first > $last );                                        # switch start/stop position

	push( @features, { ID => $ID, first => $first, last => $last, strand => $strand, contig => $contig_number } );
}

if ( !$csv->eof )
{
	$csv->error_diag;
	die $csv->error_input;
}
close $annotation_io;

die "There is a problem with the annotation file.\nStopped" if ( @features == 0 );
die "Anchor $backbone_id doesn't appear in annotation file \"$annotation_file\".\nStopped"
  if ( $backbone_id and none { $_->{ID} eq $backbone_id } @features );
say "\"$annotation_file\"";
say "--> " . @features . " features";

###################
# Genome sequence #
###################
sub genome_sequence { }

say "\n(2) Processing genomic SEQUENCE …";

my %sequences;    # contig => sequence
my $sequence_io = Bio::SeqIO->new( -file => $genome_file );

while ( my $sequence = $sequence_io->next_seq() )
{
	my $contig_id = $sequence->display_id();
	die "Looks like there is more then one contig $contig_id. \""
	  . $sequence->display_id()
	  . "\" looks like another incarnation of the already processed contig number $contig_id, at least to me.\nStopped"
	  if exists $sequences{$contig_id};
	$sequences{$contig_id} = lc( $sequence->seq() );
}
say "\"$genome_file\"";
say "--> " . keys(%sequences) . " contigs";

#############
# Promoters #
#############
sub promoters { }

say "\n(3) Processing PROMOTERS …";

my @promoters;               # {ID, first, last, contig}
my $backbone_promoter_nr;    # number of promoter of backbone locus
my $backbone_contig;         # yet unkown, if genome-wide

my @faulty_promoters_short;  # to console
my @faulty_promoters_long;   # to file

my $contig_number;
my $contig_sequence;
my $contig_length;

my $promoter_sequences_file = $dir . "PROMOTERS/all_promoter_sequences.fasta";
my $promoter_positions_file = $dir . "PROMOTERS/all_promoter_positions.csv";

my $skip = 0;                # indicator for "special case 9"

# extract promoters from the genome sequence
make_path( ( fileparse($promoter_sequences_file) )[1] );
open( my $promoter_sequences_io, ">", $promoter_sequences_file )
  or die "Cannot write to promoter sequences file \"$promoter_sequences_file\".\n", $!;
open( my $promoter_positions_io, ">", $promoter_positions_file )
  or die "Cannot write to promoter positions file \"$promoter_positions_file\".\n", $!;
say "Writing sequences to \"$promoter_sequences_file\"";
say "Writing positions to \"$promoter_positions_file\"";
say $promoter_positions_io join( "\t", "#", "promoter", "start", "end", "length", "contig" );    # print head to positions file

# compute borders of the promoter for each locus
for ( 0 .. $#features )
{
	$contig_number   = $features[$_]{contig};
	$contig_sequence = \$sequences{$contig_number};

	if ( !$sequences{$contig_number} )
	{
		die "[ERROR] Contig \"$contig_number\" (gene $features[$_]{ID}) does not occur in the genome sequence file \"$genome_file\".\nStopped";
	}
	else
	{
		$contig_length = length( $sequences{$contig_number} );
	}

	# promoter region = complete intergenic region
	if ($complete_intergenic)
	{
		if ( $_ == 0 and $features[$_]{strand} eq "+" )    # first gene in the contig
		{
			$upstream_tss = $features[$_]{first} - 1;
		}
		elsif ( $_ == $#features and $features[$_]{strand} eq "-" )    # last gene in contig
		{
			$upstream_tss = $contig_length - $features[$_]{last} - 1;
		}
		elsif ( $features[$_]{strand} eq "-" )
		{
			$upstream_tss = $features[ $_ + 1 ]{first} - $features[$_]{last} - 1;
		}
		else
		{
			$upstream_tss = $features[$_]{first} - $features[ $_ - 1 ]{last} - 1;
		}
	}

	# two genes share the same promotor > special case 9
	# did computation with first gene > skip second gene
	if ($skip)
	{
		$skip = 0;
	}

	# only one locus within current contig
	elsif (
		( $_ == 0 and $contig_number ne $features[ $_ + 1 ]{contig} )    # first locus
		or (    $_ == $#features
			and $contig_number ne $features[ $_ - 1 ]{contig} )          # last locus
		or (    $contig_number ne $features[ $_ - 1 ]{contig}
			and $contig_number ne $features[ $_ + 1 ]{contig} )
	  )                                                                  # elsewhere
	{
		if ( $features[$_]{strand} eq "+" )
		{
			if (    ( ( $features[$_]{first} - $upstream_tss ) >= 0 )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{                                                            # 1
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( ( $features[$_]{first} - $upstream_tss ) < 0 )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{                                                            # 2
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => 1,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( ( $features[$_]{first} - $upstream_tss ) >= 0 )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{                                                            # 3
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{last}
					}
				);
			}
			elsif ( ( ( $features[$_]{first} - $upstream_tss ) < 0 )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{                                                            # 7
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => 1,
						last  => $features[$_]{last}
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		elsif ( $features[$_]{strand} eq "-" )
		{
			if (    ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) <= $contig_length ) )
			{    # 4
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) > $contig_length ) )
			{    # 5
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $contig_length
					}
				);
			}
			elsif ( ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) <= $contig_length ) )
			{    # 6
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ + 1 ]{first} >= $features[$_]{last} - $upstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) > $contig_length ) )
			{    # 8
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $contig_length
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		else
		{
			die "Strand-problem at $features[$_]{ID}.\nStopped";
		}
	}

	# first gene of the contig
	# AND NOT special case 9
	# [ $_ -1 ] not present (first gene of genome) or
	# current contig number differs from last contig number
	elsif (
		( $_ == 0 or $contig_number ne $features[ $_ - 1 ]{contig} and $_ != $#features )
		and !(
			    ( ( $features[$_]{strand} eq "-" ) and ( $features[ $_ + 1 ]{strand} eq "+" ) )
			and ( $features[$_]{last} + $upstream_tss >= $features[ $_ + 1 ]{first} - $upstream_tss )
		)
	  )
	{
		if ( $features[$_]{strand} eq "+" )
		{
			if (    ( ( $features[$_]{first} - $upstream_tss ) >= 0 )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{    # 1
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( ( $features[$_]{first} - $upstream_tss ) < 0 )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{    # 2
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => 1,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( ( $features[$_]{first} - $upstream_tss ) >= 0 )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{    # 3
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{last}
					}
				);
			}
			elsif ( ( ( $features[$_]{first} - $upstream_tss ) < 0 )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{    # 7
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => 1,
						last  => $features[$_]{last}
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		elsif ( $features[$_]{strand} eq "-" )
		{
			if (    ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( $features[ $_ + 1 ]{first} > $features[$_]{last} + $upstream_tss ) )
			{    # 4
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( $features[ $_ + 1 ]{first} <= $features[$_]{last} + $upstream_tss ) )
			{    # 5
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $features[ $_ + 1 ]{first} - 1
					}
				);
			}
			elsif ( ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss )
				and ( $features[ $_ + 1 ]{first} > $features[$_]{last} + $upstream_tss ) )
			{    # 6
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ + 1 ]{first} <= $features[$_]{last} + $upstream_tss )
				and ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss ) )
			{    # 8
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $features[ $_ + 1 ]{first} - 1
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		else
		{
			die "Strand-problem at $features[$_]{ID}.\nStopped";
		}
	}

	# last gene of contig
	# [ $_ +1 ] not present (last gene of genome) or
	# current contig number differs from next contig number
	elsif ( ( $_ == $#features or $contig_number ne $features[ $_ + 1 ]{contig} )
		and !$skip )
	{
		if ( $features[$_]{strand} eq "+" )
		{
			if (    ( $features[ $_ - 1 ]{last} < $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{    # 1
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ - 1 ]{last} >= $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{    # 2
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[ $_ - 1 ]{last} + 1,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ - 1 ]{last} < $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{    # 3
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{last}
					}
				);
			}
			elsif ( ( $features[ $_ - 1 ]{last} >= $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{    # 7
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[ $_ - 1 ]{last} + 1,
						last  => $features[$_]{last}
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		elsif ( $features[$_]{strand} eq "-" )
		{
			if (    ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) <= $contig_length ) )
			{    # 4
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) > $contig_length ) )
			{    # 5
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $contig_length
					}
				);
			}
			elsif ( ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) <= $contig_length ) )
			{    # 6
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ + 1 ]{first} <= $features[$_]{last} + $upstream_tss )
				and ( ( $features[$_]{last} + $upstream_tss ) > $contig_length ) )
			{    # 8
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $contig_length
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		else
		{
			die "Strand-problem at $features[$_]{ID}";
		}
	}

	# special-case 9
	elsif ( ( ( $features[$_]{strand} eq "-" ) and ( $features[ $_ + 1 ]{strand} eq "+" ) )
		and ( $features[$_]{last} + $upstream_tss >= $features[ $_ + 1 ]{first} - $upstream_tss ) )
	{
		if (    ( $features[$_]{last} > $features[$_]{first} + $downstream_tss )
			and ( $features[$_]{first} < $features[$_]{last} - $downstream_tss ) )
		{    # 9.1+4
			push(
				@promoters,
				{
					ID    => $features[$_]{ID} . "+" . $features[ $_ + 1 ]{ID},
					first => $features[$_]{last} - $downstream_tss,
					last  => $features[ $_ + 1 ]{first} + $downstream_tss
				}
			);
		}
		elsif ( ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
			and ( $features[ $_ + 1 ]{first} + $downstream_tss >= $features[ $_ + 1 ]{last} ) )
		{    # 9.3+4
			push(
				@promoters,
				{
					ID    => $features[$_]{ID} . "+" . $features[ $_ + 1 ]{ID},
					first => $features[$_]{last} - $downstream_tss,
					last  => $features[ $_ + 1 ]{last}
				}
			);
		}
		elsif ( ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss )
			and ( $features[ $_ + 1 ]{last} > $features[ $_ + 1 ]{first} + $downstream_tss ) )
		{    # 9.1+6
			push(
				@promoters,
				{
					ID    => $features[$_]{ID} . "+" . $features[ $_ + 1 ]{ID},
					first => $features[$_]{first},
					last  => $features[ $_ + 1 ]{first} + $downstream_tss
				}
			);
		}
		elsif ( ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss )
			and ( $features[ $_ + 1 ]{first} + $downstream_tss >= $features[ $_ + 1 ]{last} ) )
		{    # 9.3+6
			push(
				@promoters,
				{
					ID    => $features[$_]{ID} . "+" . $features[ $_ + 1 ]{ID},
					first => $features[$_]{first},
					last  => $features[ $_ + 1 ]{last}
				}
			);
		}
		else
		{
			die "Problem with promoter $features[$_]{ID}.\nStopped";
		}
		$skip = 1;
	}

	# normal-cases
	elsif ( !$skip )
	{
		if ( $features[$_]{strand} eq "+" )
		{
			if (    ( $features[ $_ - 1 ]{last} < $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{    # 1
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ - 1 ]{last} >= $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{last} > $features[$_]{first} + $downstream_tss ) )
			{    # 2
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[ $_ - 1 ]{last} + 1,
						last  => $features[$_]{first} + $downstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ - 1 ]{last} < $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{    # 3
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first} - $upstream_tss,
						last  => $features[$_]{last}
					}
				);
			}
			elsif ( ( $features[ $_ - 1 ]{last} >= $features[$_]{first} - $upstream_tss )
				and ( $features[$_]{first} + $downstream_tss >= $features[$_]{last} ) )
			{    # 7
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[ $_ - 1 ]{last} + 1,
						last  => $features[$_]{last}
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		elsif ( $features[$_]{strand} eq "-" )
		{
			if (    ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( $features[ $_ + 1 ]{first} > $features[$_]{last} + $upstream_tss ) )
			{    # 4
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[$_]{first} < $features[$_]{last} - $downstream_tss )
				and ( $features[ $_ + 1 ]{first} <= $features[$_]{last} + $upstream_tss ) )
			{    # 5
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{last} - $downstream_tss,
						last  => $features[ $_ + 1 ]{first} - 1
					}
				);
			}
			elsif ( ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss )
				and ( $features[ $_ + 1 ]{first} > $features[$_]{last} + $upstream_tss ) )
			{    # 6
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $features[$_]{last} + $upstream_tss
					}
				);
			}
			elsif ( ( $features[ $_ + 1 ]{first} <= $features[$_]{last} + $upstream_tss )
				and ( $features[$_]{first} >= $features[$_]{last} - $downstream_tss ) )
			{    # 8
				push(
					@promoters,
					{
						ID    => $features[$_]{ID},
						first => $features[$_]{first},
						last  => $features[ $_ + 1 ]{first} - 1
					}
				);
			}
			else
			{
				die "Problem with promoter $features[$_]{ID}.\nStopped";
			}
		}
		else
		{
			die "Strand-problem at $features[$_]{ID}.\nStopped";
		}
	}

	# negative start position
	# or stop position "beyond" contig
	# --> might happen in very small contigs
	$promoters[-1]{first} = 1 if ( $promoters[-1]{first} < 1 );
	$promoters[-1]{last} = $contig_length if ( $promoters[-1]{last} > $contig_length );

	# write positions and sequences to files
	if ( !$skip )
	{
		my $promoter_length = $promoters[-1]{last} - $promoters[-1]{first} + 1;
		my $promoter_sequence = substr( $$contig_sequence, $promoters[-1]{first} - 1, $promoter_length );

		my $min_promoter_length = 6;
		my $max_promoter_length = ( $upstream_tss + $downstream_tss ) * 2 + 1;
		my $invalid_promoter_sequence;

		# check for valid promoter length
		if ( $promoter_length < $min_promoter_length or $promoter_length > $max_promoter_length )
		{
			$invalid_promoter_sequence = "length";
		}
		# check for existence of a, c, g and t for at least one time in the promoter sequence
		elsif ( index( $promoter_sequence, "a" ) == -1 )
		{
			$invalid_promoter_sequence = "a";
		}
		elsif ( index( $promoter_sequence, "c" ) == -1 )
		{
			$invalid_promoter_sequence = "c";
		}
		elsif ( index( $promoter_sequence, "g" ) == -1 )
		{
			$invalid_promoter_sequence = "g";
		}
		elsif ( index( $promoter_sequence, "t" ) == -1 )
		{
			$invalid_promoter_sequence = "t";
		}

		if ($invalid_promoter_sequence)
		{
			my $error_long;
			if ( $invalid_promoter_sequence eq "length" )
			{
				$error_long = "strange length - promoter ignored:\n";
			}
			else
			{
				$error_long = "$invalid_promoter_sequence is missing - promoter ignored:\n";    # actually SiTaR doesn't like such missings
			}
			$error_long .= join(
				"\n",
				"promoter " . $promoters[-1]{ID},
				#				"promoter number " . $#promoters ,
				"start " . $promoters[-1]{first},
				"stop " . $promoters[-1]{last},
				"length $promoter_length",
				"contig $contig_number",
			) . "\n";

			die "The promoter sequence of the anchor gene is invalid." . "\n$error_long" . "\nStopped"
			  if ( $backbone_id and $promoters[-1]{ID} =~ m/(^\Q$backbone_id\E$)|(^\Q$backbone_id+\E)|(\Q+$backbone_id\E$)/ );

			push( @faulty_promoters_short, $promoters[-1]{ID} );
			push( @faulty_promoters_long,  $error_long );

			pop(@promoters);    # remove last promoter
		}
		else
		{
			# print promoter positions
			say $promoter_positions_io
			  join( "\t", ( $#promoters + 1 ), $promoters[-1]{ID}, $promoters[-1]{first}, $promoters[-1]{last}, $promoter_length, $contig_number );

			# print promoter sequences
			$promoter_sequence =~ s/(.{60})/$1\n/g;    # maximum of 60 characters per line
			say $promoter_sequences_io ">" . $promoters[-1]{ID} . "__" . $promoter_length . "_bp__contig_" . $contig_number;
			say $promoter_sequences_io $promoter_sequence . "\n";

			# save promoters' contig number
			$promoters[-1]{contig} = $contig_number;

			# if motif prediction (MEME) enabled:
			# backbone_nr not saved so far and current promoter ID includes backbone_id
			# --> save number and contig of backbone gene
			# warning: the promoter ID is not equal to the backbone ID, nor the promoter header! (it's somewhere "in between")
			if ( $backbone_id and !$backbone_promoter_nr and $promoters[-1]{ID} =~ m/(^\Q$backbone_id\E$)|(^\Q$backbone_id+\E)|(\Q+$backbone_id\E$)/ )
			{
				$backbone_promoter_nr = @promoters;
				$backbone_contig      = $contig_number;
			}
		}
	}

	# check if promoter names are unique
	if ( @promoters >= 2 and $promoters[-1]{ID} eq $promoters[-2]{ID} )
	{
		die "[ERROR] Promoter \""
		  . $promoters[-1]{ID}
		  . "\" occurs twice (promoter numbers "
		  . ( @promoters - 1 ) . " and "
		  . (@promoters)
		  . "). Please check your annotation file \"$annotation_file\". Most likely it contains overlapping gene annotations.\nStopped";
	}
}
close $promoter_sequences_io;
close $promoter_positions_io;

my $faulty_promoters_file = $dir . "faulty_promoters";
if (@faulty_promoters_short)
{
	say "\n[WARNING] Faulty promoters:\n" . join( "\n", @faulty_promoters_short );

	open( my $faulty_promoters_io, ">", $faulty_promoters_file )
	  or die "Cannot write to faulty promoters file \"$faulty_promoters_file\".\n", $!;
	say $faulty_promoters_io join( "\n", @faulty_promoters_long );
	close $faulty_promoters_io;

	say @faulty_promoters_short . " ignored promoter(s) in total, see file \"$faulty_promoters_file\"";
}

say "--> " . @promoters . " promoter regions";

# if cluster prediction enabled --> check number of promoter regions on the anchor gene contig and warn if small
if ($cluster_prediction)
{
	my $promoter_regions_on_backbone_contig = 0;
	map { $promoter_regions_on_backbone_contig++ if ( $_->{contig} eq $backbone_contig ) } @promoters;

	if ( $promoter_regions_on_backbone_contig < 3 )
	{
		say
"[ERROR] The contig \"$backbone_contig\" of the anchor gene \"$backbone_id\" has less than 3 promoter regions. Cluster prediction not possible on such a small contig.";
		exit;
	}
	elsif ( $promoter_regions_on_backbone_contig < 40 )
	{
		say
"[WARNING] The contig \"$backbone_contig\" of the anchor gene \"$backbone_id\" has only $promoter_regions_on_backbone_contig promoter regions. Predictions on small contigs may lead to partial clusters (cut by the contig border).";
	}
}

############################################
# Predict overrepresented motifs with MEME #
############################################
sub motifs { }

# MEME/FIMO:
my @predicted_motifs;    # { motif (plus_minus), motif-score, count [later] }
my @predicted_motifs_sort;

# SiTaR:
my @given_motifs;        # { FASTA header ID, count [later] }
my @given_motifs_sort;

if ($meme)
{
	say "\n(4) Predicting overrepresented MOTIFS around the anchor gene …";

	@predicted_motifs = meme( $dir, $promoter_sequences_file, $cluster_name, $backbone_id, $backbone_contig, $verbose, \@meme_parameters, $num_cpus );
	@predicted_motifs_sort = sort { $a->{score} <=> $b->{score} } @predicted_motifs;

	if ( !@predicted_motifs )
	{
		say "No motifs --> No cluster prediction for anchor gene \"$backbone_id\" (promoter #$backbone_promoter_nr).   :(";
		exit;
	}
	else
	{
		say "--> best motif score:      " . $predicted_motifs_sort[0]{motif} . "   " . $predicted_motifs_sort[0]{score};
		say "--> 2. best motif score:   " . $predicted_motifs_sort[1]{motif} . "   " . $predicted_motifs_sort[1]{score} if $predicted_motifs[1];
		say "--> 3. best motif score:   " . $predicted_motifs_sort[2]{motif} . "   " . $predicted_motifs_sort[2]{score} if $predicted_motifs[2];
	}
}
elsif ($sitar)
{
	say "\n(4) Predicting overrepresented MOTIFS around the anchor gene …";
	say "Skipping --> No need to predict motifs with MEME.";
	say "Instead, running SiTaR with binding sites given in file \"$sitar\"";

	my $sequence_io = Bio::SeqIO->new( -file => $sitar, -format => "fasta" );
	while ( my $sequence = $sequence_io->next_seq() )
	{
		die "[ERROR] The sequence IDs of \"$sitar\" are not unique. \"" . $sequence->display_id() . "\" occurs at least twice.\nStopped"
		  if ( any { $_->{motif} eq $sequence->display_id() } @given_motifs );
		push( @given_motifs, { motif => $sequence->display_id() } );
	}
}
else
{
	# motif prediction with MEME not activated and motif search with SiTaR not activated, too
	# --> skip any further steps and exit
	exit;
}

#############################
# Search motifs genome-wide #
# with FIMO or SiTaR        #
#############################

if ($fimo)
{
	say "\n(5) Running FIMO …";
	make_path("$dir$cluster_name/fimo/");
	fimo();
}
elsif ($sitar)
{
	say "\n(5) Running SiTaR …";
	make_path("$dir$cluster_name/sitar/");
	sitar();
}
else
{
	# binding site search with FIMO or SiTaR not activated
	# skip any further steps and exit
	exit;
}

sub fimo
{
	my $read_from_previous = 0;

	# multiple FIMO instances running in parallel, if at least 2 CPUs, no forks if less CPUs
	my $pm;    # parallel forks
	( $num_cpus >= 2 ) ? ( $pm = new Parallel::ForkManager($num_cpus) ) : ( $pm = new Parallel::ForkManager(0) );

	# progress bar (sort of), but updates -- thanks to \r -- the current line, instead of printing to the next one
	my $progress = 1;
	$pm->run_on_finish( sub { print "\rProgress: " . $progress++ . "/" . @predicted_motifs if $verbose; } );

	for ( 0 .. $#predicted_motifs )
	{
		my %motif = %{ $predicted_motifs[$_] };    # { motif, score }
		my ( $plus, $minus ) = split( "_", $motif{motif} );
		my $fimo_dir = "$dir$cluster_name/fimo/$motif{motif}/";

		# check if fimo output was already present, at least once
		$read_from_previous = 1 if ( !$read_from_previous and -d $fimo_dir and -s $fimo_dir . "fimo.txt" );
        # TODO fimo.tsv if meme suite >= 5.0.0

		$pm->start and next;                       # fork
		if ( !( -d $fimo_dir and -s $fimo_dir . "fimo.txt" ) )    # motif not searched previously? --> run FIMO
        # TODO fimo.tsv if meme suite >= 5.0.0
		{
			my @args = (
				# "-verbosity", 1, "-motif", 1, "-thresh", $fimo, "-oc", $fimo_dir, "$dir$cluster_name/meme/" . $plus . "_" . $minus . "/meme.html",
				# "-motif 1" is not recgonised anymore
				# have to provide "-motif <id>", like YTCGCTYYCCGA, in never FIMO versions
				# --> omit this arg, because our MEME output contains only one motif, anyway
				"-verbosity", 1, "-thresh", $fimo, "-oc", $fimo_dir, "$dir$cluster_name/meme/" . $plus . "_" . $minus . "/meme.html",
				$promoter_sequences_file
			);
			system( "fimo", @args );
		}
		$pm->finish;                                              # exit child process
	}

	$pm->wait_all_children;                                       # wait for all child processes to finish

	say "\n(No new search for at least one motif. Using existing FIMO files.)" if ( $read_from_previous and $verbose );
}

sub sitar
{
	# TODO SiTaR, also introduce score cut-off, like FIMO?

	my $sitar_dir = "$dir$cluster_name/sitar/";
	my @args = ( "mismatch=$mismatches", "tf=$sitar", "seq=$promoter_sequences_file", "output=$sitar_dir" . "sitar.csv" );
	system( "sitar", @args );                                     # progress bar provided by sitar, thanks
}

##############################
# Binding sites per promoter #
##############################
sub binding_sites_per_promoter { }

say "\n(6) Counting BINDING SITES per promoter …";

# @counts: number of binding sites per promoter
# each motif and corresponding cluster prediction will get his own @counts array ref
my @counts_refs;
my $bs_file;

# either ref to predicted_motifs(_sort) from FIMO or ref to given_motifs(_sort) from SiTaR
my $motifs;
my $motifs_sort;

$motifs = \@predicted_motifs if $fimo;
$motifs = \@given_motifs     if $sitar;

printf( "%7s   %s\n", "motif", "frequency" ) if $verbose;

for ( @{$motifs} )
{
	printf( "%7s   ", $_->{motif} ) if $verbose;

	$bs_file = "$dir$cluster_name/fimo/$_->{motif}/fimo.txt" if $fimo; # TODO fimo.tsv if meme suite >= 5.0.0
	$bs_file = "$dir$cluster_name/sitar/sitar.csv"           if $sitar;

	# open file with binding site positions
	open( my $bs_io, "<", $bs_file ) or die "Cannot read from binding sites file \"$bs_file\".\n", $!;

	# check version of binding sites file
	my $header = <$bs_io>;

	my $newer_fimo_version;
	if ( index( $header, "# motif_id" ) == 0 )
	{
		# new header (v4.11.4 and newer):
		# "# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence"
		$newer_fimo_version = 1;
	}
	elsif ( index( $header, "#pattern name" ) == 0 )
	{
		# old header:
		# "#pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence"
		$newer_fimo_version = 0;
	}
	else
	{
		die "[ERROR] Unrecognized format of binding sites file \"$bs_file\".\nStopped";
	}

	# extract locus names
	my @promoters_with_bs;
	# don't need to skip, because of header check above
	# <$bs_io>;    # skip first line (comment, header)
	while ( my $row = $csv->getline($bs_io) )
	{
		# FIMO bs file contains results for only one motif
		if ($fimo)
		{
			# save promoter ID ("sequence name") and FIMO p-value, also crop "__x_bp__contig_x" part
			if ($newer_fimo_version)
			{
				push( @promoters_with_bs, { ID => substr( $row->[2], 0, index( $row->[2], "__" ) ), score => $row->[7] } );
			}
			else
			{
				push( @promoters_with_bs, { ID => substr( $row->[1], 0, index( $row->[1], "__" ) ), score => $row->[6] } );
			}
		}
		# SiTaR bs file contains mixed results for all motifs
		elsif ( $sitar and $row->[1] eq $_->{motif} )
		{
			# save promoter ID ("sequence name") and SiTaR score, also crop "__x_bp__contig_x" part
			push( @promoters_with_bs, { ID => substr( $row->[0], 0, index( $row->[0], "__" ) ), score => $row->[5] } );
		}
	}

	if ( !$csv->eof )
	{
		$csv->error_diag;
		die $csv->error_input;
	}
	close $bs_io;

	# fast variant to fill @counts
	# (instead of looping through @promoters and @promoters_with_bs and checking IDs)
	# but still one of the most time-consuming parts of CASSIS, except MEME and FIMO
	my %promoters_temp;
	for ( 0 .. $#promoters )    # fill hash with promoter IDs as keys
	{
		$promoters_temp{ $promoters[$_]{ID} } = [ 0, $_, "-" ];
		#                |                        |  |    |
		#                |                        |  |    value 2: p-value (FIMO) / score (SiTaR)
		#                key: promoter ID         |  value 1: promoter number
		#                                         value 0: binding sites
	}

	for (@promoters_with_bs)
	{
		# count number of binding sites per promoter
		$promoters_temp{ $_->{ID} }[0]++;

		if ($fimo)
		{
			# save LOWEST p-value per promoter
			$promoters_temp{ $_->{ID} }[2] = $_->{score} if ( $promoters_temp{ $_->{ID} }[2] eq "-" or $_->{score} < $promoters_temp{ $_->{ID} }[2] );
		}
		elsif ($sitar)
		{
			# save HIGHEST score per promoter
			$promoters_temp{ $_->{ID} }[2] = $_->{score} if ( $promoters_temp{ $_->{ID} }[2] eq "-" or $_->{score} > $promoters_temp{ $_->{ID} }[2] );
		}
	}

	# read from right to left:
	# 1. get all values from hash --> array ref [binding sites, promoter numbers, p-values/score]
	# 2. sort by promoter numbers --> will be the array index
	# 3. push only number of binding sites and p-value/score into array
	#
	my @counts = map { { bs => $_->[0], score => $_->[2] } } sort { $a->[1] <=> $b->[1] } values %promoters_temp;    # somewhat "slow"

	# fraction of promoters with binding sites compared to all promoters
	my $num_promoters_with_bs = grep { $_->{bs} != 0 } @counts;
	my $frequency = $num_promoters_with_bs / @counts;
	printf(
		"%5s binding sites in   %14s   %23s\n",
		scalar(@promoters_with_bs),
		"$num_promoters_with_bs promoters",
		sprintf( "%.2f", $frequency * 100 ) . " % of all promoters"
	) if $verbose;

	# write binding sites per promoter to file
	my $bs_per_promoter_file;
	if ($fimo)
	{
		$bs_per_promoter_file = ( fileparse($bs_file) )[1] . "bs_per_promoter_fimo.csv";
	}
	elsif ($sitar)
	{
		my $dir = ( fileparse($bs_file) )[1] . $_->{motif};
		mkdir($dir) or die "Cannot create directory \"$dir\".\n", $! unless -e $dir;
		$bs_per_promoter_file = "$dir/bs_per_promoter_sitar.csv";
	}

	open( my $bs_per_promoter_io, ">", $bs_per_promoter_file )
	  or die "Cannot write to binding sites per promoter file \"$bs_per_promoter_file\".\n", $!;

	# table head
	say $bs_per_promoter_io join( "\t", "#", "promoter", "binding sites", "lowest p-value", "contig" ) if $fimo;
	say $bs_per_promoter_io join( "\t", "#", "promoter", "binding sites", "highest score",  "contig" ) if $sitar;

	map { say $bs_per_promoter_io join( "\t", $_ + 1, $promoters[$_]{ID}, $counts[$_]{bs}, $counts[$_]{score}, $promoters[$_]{contig} ) }
	  ( 0 .. $#counts );    # "slow", but Text::XS is even slower and we have to write the binding sites to a file anyway
	close $bs_per_promoter_io;

	@counts = ("too high") if ( $frequency > $frequency_cutoff );    # frequency too high --> mark invalid
	@counts = ("too low")  if ( $frequency == 0 );                   # frequency too low --> mark invalid

	push( @counts_refs, \@counts );
}

######################
# cluster prediction #
######################
sub cluster_prediction { }

exit if !$cluster_prediction;

say "\n(7) Computing CLUSTER PREDICTIONS …";

my @skipped;
my $too_many;
my $valid_motifs;

for ( 0 .. $#{$motifs} )
{
	if ( $counts_refs[$_][0] eq "too high" )
	{
		push( @skipped, sprintf( "%7s   Found too many binding sites. *", $$motifs[$_]{motif} ) );
		$$motifs[$_]     = "";    # mark motif, which resulted in too many binding sites, to get removed
		$counts_refs[$_] = "";
		$too_many        = 1;
	}
	elsif ( $counts_refs[$_][0] eq "too low" )
	{
		push( @skipped, sprintf( "%7s   Didn't find any binding sites.", $$motifs[$_]{motif} ) );
		$$motifs[$_]     = "";    # mark motif, which resulted in no binding sites, to get removed
		$counts_refs[$_] = "";
	}
	elsif ( $counts_refs[$_][ $backbone_promoter_nr - 1 ]{bs} == 0 )
	{
		push( @skipped, sprintf( "%7s   Promoter of anchor gene has no binding site.", $$motifs[$_]{motif} ) );
		$$motifs[$_]     = "";    # mark motif, which resulted in no binding site in the anchor promoter, to get removed
		$counts_refs[$_] = "";
	}
}

@counts_refs = grep { $_ ne "" } @counts_refs;    # cleanup @counts_refs

if ($fimo)
{
	@predicted_motifs = grep { $_ ne "" } @predicted_motifs;    # also cleanup @predicted_motifs (FIMO)
	@predicted_motifs =
	  map { { motif => $predicted_motifs[$_], counts => $counts_refs[$_] } } ( 0 .. $#predicted_motifs );    # connect motifs and their counts
	@predicted_motifs_sort = sort { $a->{motif}{score} <=> $b->{motif}{score} } @predicted_motifs;           # sort connected structure by motif score

	# refresh references (necessary?)
	$motifs      = \@predicted_motifs;
	$motifs_sort = \@predicted_motifs_sort;
}
elsif ($sitar)
{
	@given_motifs = grep { $_ ne "" } @given_motifs;                                                         # also cleanup @given_motifs (SiTaR)
	@given_motifs =
	  map { { motif => $given_motifs[$_], counts => $counts_refs[$_] } } ( 0 .. $#given_motifs );            # connect motifs and their counts
	      #	@given_motifs_sort = sort { $a->{motif}{score} <=> $b->{motif}{score} } @given_motifs;           # sort connected structure by motif score
	@given_motifs_sort = @given_motifs;    # SiTaR --> no motif scores --> nothing to sort by

	# refresh references (necessary?)
	$motifs      = \@given_motifs;
	$motifs_sort = \@given_motifs_sort;
}

make_path( $dir . $cluster_name . "/CLUSTER/" );
my $cluster_predictions_file = $dir . $cluster_name . "/CLUSTER/" . $backbone_id . "_all_predictions.csv";
open( my $cluster_predictions_io, ">", $cluster_predictions_file )
  or die "Cannot write to cluster predictions file \"$cluster_predictions_file\".\n", $!;
$csv->print( $cluster_predictions_io, [ "# motif", "motif score", "first promoter", "last promoter", "first gene", "last gene", "length" ] );

printf( "%7s   %8s   %11s\n", "motif", "score", "prediction" ) if $verbose;

my @cluster_predictions;
for ( 0 .. $#{$motifs_sort} )
{
	my $motif  = $$motifs_sort[$_]{motif};
	my @counts = @{ $$motifs_sort[$_]{counts} };

	my $first = $backbone_promoter_nr - 1;
	my $last  = $backbone_promoter_nr - 1;

	# upstream
	for ( my $i = $backbone_promoter_nr - 1 ; $i > 0 ; $i-- )
	{
		# promoter with binding site(s)
		# … 1 …
		if ( $counts[ $i - 1 ]{bs} >= 1 )
		{
			$first -= 1;
		}

		# no binding site, gap with length 1
		# … 1 0 1  …
		elsif ( $i - 2 >= 0 and $max_gap_length >= 1 and $counts[$i]{bs} >= 1 and $counts[ $i - 1 ]{bs} == 0 and $counts[ $i - 2 ]{bs} >= 1 )
		{
			$first -= 2;
			$i     -= 1;
		}

		# no binding site, gap with length 2
		# … 1 0 0 1  …
		elsif ( $i - 3 >= 0
			and $max_gap_length >= 2
			and $counts[$i]{bs} >= 1
			and $counts[ $i - 1 ]{bs} == 0
			and $counts[ $i - 2 ]{bs} == 0
			and $counts[ $i - 3 ]{bs} >= 1 )
		{
			$first -= 3;
			$i     -= 2;
		}

		# no binding site, gap with length 3
		# … 1 0 0 0 1  …
		elsif ( $i - 4 >= 0
			and $max_gap_length >= 3
			and $counts[$i]{bs} >= 1
			and $counts[ $i - 1 ]{bs} == 0
			and $counts[ $i - 2 ]{bs} == 0
			and $counts[ $i - 3 ]{bs} == 0
			and $counts[ $i - 4 ]{bs} >= 1 )
		{
			$first -= 4;
			$i     -= 3;
		}

		# no binding site, gap with length 4
		# … 1 0 0 0 0 1  …
		elsif ( $i - 5 >= 0
			and $max_gap_length >= 4
			and $counts[$i]{bs} >= 1
			and $counts[ $i - 1 ]{bs} == 0
			and $counts[ $i - 2 ]{bs} == 0
			and $counts[ $i - 3 ]{bs} == 0
			and $counts[ $i - 4 ]{bs} == 0
			and $counts[ $i - 5 ]{bs} >= 1 )
		{
			$first -= 5;
			$i     -= 4;
		}

		# no binding site, gap with length 5
		# … 1 0 0 0 0 0 1  …
		elsif ( $i - 6 >= 0
			and $max_gap_length >= 5
			and $counts[$i]{bs} >= 1
			and $counts[ $i - 1 ]{bs} == 0
			and $counts[ $i - 2 ]{bs} == 0
			and $counts[ $i - 3 ]{bs} == 0
			and $counts[ $i - 4 ]{bs} == 0
			and $counts[ $i - 5 ]{bs} == 0
			and $counts[ $i - 6 ]{bs} >= 1 )
		{
			$first -= 6;
			$i     -= 5;
		}

		# gap too long, stop upstream cluster extension
		else
		{
			last;
		}
	}

	# downstream
	for ( my $i = $backbone_promoter_nr - 1 ; $i < $#counts ; $i++ )
	{
		# promoter with binding site(s)
		# … 1 …
		if ( $counts[ $i + 1 ]{bs} >= 1 )
		{
			$last += 1;
		}

		# no binding site, gap with length 1
		# … 1 0 1  …
		elsif ( $i + 2 <= $#counts and $max_gap_length >= 1 and $counts[$i]{bs} >= 1 and $counts[ $i + 1 ]{bs} == 0 and $counts[ $i + 2 ]{bs} >= 1 )
		{
			$last += 2;
			$i    += 1;
		}

		# no binding site, gap with length 2
		# … 1 0 0 1  …
		elsif ( $i + 3 <= $#counts
			and $max_gap_length >= 2
			and $counts[$i]{bs} >= 1
			and $counts[ $i + 1 ]{bs} == 0
			and $counts[ $i + 2 ]{bs} == 0
			and $counts[ $i + 3 ]{bs} >= 1 )
		{
			$last += 3;
			$i    += 2;
		}

		# no binding site, gap with length 3
		# … 1 0 0 0 1  …
		elsif ( $i + 4 <= $#counts
			and $max_gap_length >= 3
			and $counts[$i]{bs} >= 1
			and $counts[ $i + 1 ]{bs} == 0
			and $counts[ $i + 2 ]{bs} == 0
			and $counts[ $i + 3 ]{bs} == 0
			and $counts[ $i + 4 ]{bs} >= 1 )
		{
			$last += 4;
			$i    += 3;
		}

		# no binding site, gap with length 4
		# … 1 0 0 0 0 1  …
		elsif ( $i + 5 <= $#counts
			and $max_gap_length >= 4
			and $counts[$i]{bs} >= 1
			and $counts[ $i + 1 ]{bs} == 0
			and $counts[ $i + 2 ]{bs} == 0
			and $counts[ $i + 3 ]{bs} == 0
			and $counts[ $i + 4 ]{bs} == 0
			and $counts[ $i + 5 ]{bs} >= 1 )
		{
			$last += 5;
			$i    += 4;
		}

		# no binding site, gap with length 5
		# … 1 0 0 0 0 0 1  …
		elsif ( $i + 6 <= $#counts
			and $max_gap_length >= 5
			and $counts[$i]{bs} >= 1
			and $counts[ $i + 1 ]{bs} == 0
			and $counts[ $i + 2 ]{bs} == 0
			and $counts[ $i + 3 ]{bs} == 0
			and $counts[ $i + 4 ]{bs} == 0
			and $counts[ $i + 5 ]{bs} == 0
			and $counts[ $i + 6 ]{bs} >= 1 )
		{
			$last += 6;
			$i    += 5;
		}

		# gap too long, stop downstream cluster extension
		else
		{
			last;
		}
	}

	# save details of the predicted cluster
	push(
		@cluster_predictions,
		{
			motif        => $motif->{motif},
			motif_score  => $motif->{score},
			first_number => $first + 1,
			last_number  => $last + 1,
			length       => $last - $first + 1,
			first_name   => $promoters[$first]{ID},
			last_name    => $promoters[$last]{ID}
		}
	);
	$cluster_predictions[-1]{motif_score} = "-" if $sitar;    # motifs/binding sites given by SiTaR input have no score assigned

	printf(
		"%7s   %8s   %11s\n",
		$cluster_predictions[-1]{motif},
		$cluster_predictions[-1]{motif_score},
		$cluster_predictions[-1]{first_number} . " - " . $cluster_predictions[-1]{last_number}
	) if ($verbose);
	$csv->print(
		$cluster_predictions_io,
		[
			$cluster_predictions[-1]{motif},       $cluster_predictions[-1]{motif_score}, $cluster_predictions[-1]{first_number},
			$cluster_predictions[-1]{last_number}, $cluster_predictions[-1]{first_name},  $cluster_predictions[-1]{last_name},
			$cluster_predictions[-1]{length}
		]
	);
}
close $cluster_predictions_io;

say "\n" . join( "\n", @skipped ) if ( $verbose and @skipped );
say "(* occurs in more than " . ( $frequency_cutoff * 100 ) . " % of all promoters in the genome)" if ( $verbose and $too_many );

###########
# Results #
###########
sub results { }

say "\n" . "(8) RESULTS …";

# save only most abundant cluster predictions for results table
my @cluster_predictions_most_abundant;

if ( !@{$motifs} or !@cluster_predictions )
{
	# @predicted_motifs or @cluster_predictions is emtpy --> not a single valid motif/prediction remains
	say "No cluster prediction for anchor gene \"$backbone_id\" (promoter #$backbone_promoter_nr).   :(";
	exit;
}
else
{
	if ($verbose)
	{
		print autoformat "See file \"$cluster_predictions_file\" for a detailed output.";
		print autoformat
"All predictions are given in the form \"upstream cluster border (first gene) - downstream cluster border (last gene)\". Motif scores and abundances are given in the form \"first/last\".";
	}

	print autoformat "Out of "
	  . @cluster_predictions
	  . " cluster prediction(s) in total, this is/are the best for the anchor gene \"$backbone_id\", promoter #$backbone_promoter_nr:";

	my %count_first_number;
	my %count_last_number;
	my %count_first_motif_score_number;
	my %count_last_motif_score_number;

	my %count_first_name;
	my %count_last_name;
	my %count_first_motif_score_name;
	my %count_last_motif_score_name;

	# count occurrences
	# and save first (== lowest) score for each prediction (@cluster_prediction is ordered by motif score)
	for (@cluster_predictions)
	{
		# count by promoter number
		$count_first_number{ $_->{first_number} }++;
		$count_last_number{ $_->{last_number} }++;
		$count_first_motif_score_number{ $_->{first_number} } = { score => $_->{motif_score}, motif => $_->{motif} }
		  if ( !defined $count_first_motif_score_number{ $_->{first_number} } );
		$count_last_motif_score_number{ $_->{last_number} } = { score => $_->{motif_score}, motif => $_->{motif} }
		  if ( !defined $count_last_motif_score_number{ $_->{last_number} } );

		# count by gene id/name
		$count_first_name{ $_->{first_name} }++;
		$count_last_name{ $_->{last_name} }++;
		$count_first_motif_score_name{ $_->{first_name} } = { score => $_->{motif_score}, motif => $_->{motif} }
		  if ( !defined $count_first_motif_score_name{ $_->{first_name} } );
		$count_last_motif_score_name{ $_->{last_name} } = { score => $_->{motif_score}, motif => $_->{motif} }
		  if ( !defined $count_last_motif_score_name{ $_->{last_name} } );
	}

	# find most abundant (highest occurrence)
	# by number or name - doesn't matter, should be the same
	my @abundances_first = sort { $a <=> $b } uniq values %count_first_number;    # only unique abundances, sorted
	my @abundances_last  = sort { $a <=> $b } uniq values %count_last_number;
	my $i = $#abundances_first;                                                   # index of last element
	my $j = $#abundances_last;

	my @warnings;                                                                 # to collect warning messages
	my $next_less_abundant = 1;                                                   # init value, to run while loop at least once

	while ($next_less_abundant)
	{
		# most abundant, second most abundant, third most … --> i/j will decrease by one each round
		my $most_first = $abundances_first[$i];
		my $most_last  = $abundances_last[$j];

		# collect most abundant promoter numbers and gene names
		my @most_abundant_first_numbers =
		  sort { $count_first_motif_score_number{$a}{score} <=> $count_first_motif_score_number{$b}{score} }
		  grep { $count_first_number{"$_"} == $most_first } keys %count_first_number;
		my @most_abundant_last_numbers =
		  sort { $count_last_motif_score_number{$a}{score} <=> $count_last_motif_score_number{$b}{score} }
		  grep { $count_last_number{"$_"} == $most_last } keys %count_last_number;
		my @most_abundant_first_names =
		  sort { $count_first_motif_score_name{$a}{score} <=> $count_first_motif_score_name{$b}{score} }
		  grep { $count_first_name{"$_"} == $most_first } keys %count_first_name;
		my @most_abundant_last_names =
		  sort { $count_last_motif_score_name{$a}{score} <=> $count_last_motif_score_name{$b}{score} }
		  grep { $count_last_name{"$_"} == $most_last } keys %count_last_name;

		# print most abundant (in all combinations, if more than one with same abundance)
		for my $first ( 0 .. $#most_abundant_first_names )
		{
			for my $last ( 0 .. $#most_abundant_last_names )
			{
				my $first_name = ( $most_abundant_first_names[$first] =~ m/^([^\+]+)\+?/ )[0];
				my $last_name  = ( $most_abundant_last_names[$last]   =~ m/\+?([^\+]+)$/ )[0];

				# save cluster prediction
				# --> pretty printing afterwards
				# --> checking inclusion of faulty promoters afterwards
				push(
					@cluster_predictions_most_abundant,
					[
						$first_name,                                                                       # first gene name
						$last_name,                                                                        # last gene name
						$most_abundant_first_numbers[$first],                                              # first promoter number
						$most_abundant_last_numbers[$last],                                                # last promoter number
						$most_abundant_last_numbers[$last] - $most_abundant_first_numbers[$first] + 1,     # cluster length
						$count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{score},    # left border motif score
						$count_last_motif_score_number{ $most_abundant_last_numbers[$last] }{score},       # right border motif score
						$most_first, $most_last, scalar(@cluster_predictions)                              # abundance (left, right, total)
					]
				);

				# copy MEME and FIMO or SiTaR files of predicted cluster to new folder
				# same motif for upstream and downstream cluster border (always true for SiTaR)
				if ( $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif} eq
					$count_last_motif_score_number{ $most_abundant_last_numbers[$last] }{motif} )
				{
					if ($fimo)                                                                             # implies MEME
					{
						my $source_dir_meme =
						  $dir . $cluster_name . "/meme/" . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};
						my $source_dir_fimo =
						  $dir . $cluster_name . "/fimo/" . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};
						my $dest_dir =
						    $dir
						  . $cluster_name
						  . "/CLUSTER/"
						  . $first_name . "_to_"
						  . $last_name . "/"
						  . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};

						dircopy( $source_dir_meme, "$dest_dir/meme" );
						dircopy( $source_dir_fimo, "$dest_dir/fimo" );
						move( "$dest_dir/meme/binding_sites.fasta",      $dest_dir );
						move( "$dest_dir/meme/promoters.fasta",          $dest_dir );
						move( "$dest_dir/fimo/bs_per_promoter_fimo.csv", $dest_dir );
					}
					elsif ($sitar)
					{
						my $source_dir =
						  $dir . $cluster_name . "/sitar/" . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};
						my $dest_dir =
						    $dir
						  . $cluster_name
						  . "/CLUSTER/"
						  . $first_name . "_to_"
						  . $last_name . "/"
						  . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};

						dircopy( $source_dir, $dest_dir );
						copy( $bs_file, $dest_dir );
					}
				}

				# different motifs (only applies to MEME/FIMO)
				else
				{
					# upstream
					my $source_dir_meme =
					  $dir . $cluster_name . "/meme/" . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};
					my $source_dir_fimo =
					  $dir . $cluster_name . "/fimo/" . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};
					my $dest_dir =
					    $dir
					  . $cluster_name
					  . "/CLUSTER/"
					  . $first_name . "_to_"
					  . $last_name
					  . "/upstream_border_"
					  . $count_first_motif_score_number{ $most_abundant_first_numbers[$first] }{motif};
					dircopy( $source_dir_meme, "$dest_dir/meme" );
					dircopy( $source_dir_fimo, "$dest_dir/fimo" );
					move( "$dest_dir/meme/binding_sites.fasta",      $dest_dir );
					move( "$dest_dir/meme/promoters.fasta",          $dest_dir );
					move( "$dest_dir/fimo/bs_per_promoter_fimo.csv", $dest_dir );

					# downstream
					$source_dir_meme = $dir . $cluster_name . "/meme/" . $count_last_motif_score_number{ $most_abundant_last_numbers[$last] }{motif};
					$source_dir_fimo = $dir . $cluster_name . "/fimo/" . $count_last_motif_score_number{ $most_abundant_last_numbers[$last] }{motif};
					$dest_dir =
					    $dir
					  . $cluster_name
					  . "/CLUSTER/"
					  . $first_name . "_to_"
					  . $last_name
					  . "/downstream_border_"
					  . $count_last_motif_score_number{ $most_abundant_last_numbers[$last] }{motif};
					dircopy( $source_dir_meme, "$dest_dir/meme" );
					dircopy( $source_dir_fimo, "$dest_dir/fimo" );
					move( "$dest_dir/meme/binding_sites.fasta",      $dest_dir );
					move( "$dest_dir/meme/promoters.fasta",          $dest_dir );
					move( "$dest_dir/fimo/bs_per_promoter_fimo.csv", $dest_dir );
				}

			}
		}

		# check for "strange" predictions --> "[Warning] …"
		$next_less_abundant = 0;
		for (@most_abundant_first_numbers)
		{
			my $first = $_ - 1;    # revert to zero-based index

			for (@most_abundant_last_numbers)
			{
				my $last = $_ - 1;    # revert to zero-based index

				# cluster is spanning multiple contigs!
				if ( $promoters[$first]{contig} ne $promoters[$last]{contig} )
				{
					push( @warnings,
						    "\n[WARNING] Cluster \""
						  . ( $promoters[$first]{ID} =~ m/^([^\+]+)\+?/ )[0] . " - "
						  . ( $promoters[$last]{ID}  =~ m/\+?([^\+]+)$/ )[0]
						  . "\" is spanning multiple contigs:\n" );
					$next_less_abundant++;

					# search and print position where contig changes
					for ( $first .. $last - 1 )
					{
						if ( $promoters[$_]{contig} ne $promoters[ $_ + 1 ]{contig} )
						{
							push( @warnings,
								    "   "
								  . ( $promoters[$first]{ID} =~ m/^([^\+]+)\+?/ )[0] . " - "
								  . ( $promoters[$_]{ID} =~ m/\+?([^\+]+)$/ )[0]
								  . " --> contig $promoters[$_]{contig}\n" );
							push( @warnings,
								    "   "
								  . ( $promoters[ $_ + 1 ]{ID} =~ m/^([^\+]+)\+?/ )[0] . " - "
								  . ( $promoters[$last]{ID} =~ m/\+?([^\+]+)$/ )[0]
								  . " --> contig $promoters[$_+1]{contig}\n" );
							last;
						}
					}
				}

				# not on multiple contigs, but exactly at the contig border
				elsif ($first - 1 < 0
					or $first + 1 > $#promoters
					or $promoters[ $first - 1 ]{contig} ne $promoters[$first]{contig}
					or $promoters[ $last + 1 ]{contig}  ne $promoters[$last]{contig} )
				{
					push( @warnings,
						    "\n[WARNING] Cluster \""
						  . ( $promoters[$first]{ID} =~ m/^([^\+]+)\+?/ )[0] . " - "
						  . ( $promoters[$last]{ID}  =~ m/\+?([^\+]+)$/ )[0]
						  . "\" is located right next to a contig border:\n" );
					push( @warnings, "   Predictions next to the contig border may lead to partial clusters (cut by the contig border).\n" );
				}

				# not on multiple contigs, but at least near the contig border
				elsif ($first - 10 < 0
					or $first + 10 > $#promoters
					or $promoters[ $first - 10 ]{contig} ne $promoters[$first]{contig}
					or $promoters[ $last + 10 ]{contig}  ne $promoters[$last]{contig} )
				{
					push( @warnings,
						    "\n[WARNING] Cluster \""
						  . ( $promoters[$first]{ID} =~ m/^([^\+]+)\+?/ )[0] . " - "
						  . ( $promoters[$last]{ID}  =~ m/\+?([^\+]+)$/ )[0]
						  . "\" is located near a contig border:\n" );
					push( @warnings, "   Predictions near the contig border may lead to partial clusters (cut by the contig border).\n" );
				}

				# prediction (too) short
				if ( $last - $first + 1 <= 2 )
				{
					push( @warnings,
						    "\n[WARNING] Cluster \""
						  . ( $promoters[$first]{ID} =~ m/^([^\+]+)\+?/ )[0] . " - "
						  . ( $promoters[$last]{ID}  =~ m/\+?([^\+]+)$/ )[0]
						  . "\" has a length of only "
						  . ( $last - $first + 1 )
						  . " promoter(s).\n" );
					$next_less_abundant++;
				}
			}
		}

		# if all predictions of current "abundance level" aren't perfect (too short, spanning multiple contigs)
		# --> go for next less abundant contig borders, if there are any left
		if ( $next_less_abundant == @most_abundant_first_numbers * @most_abundant_last_numbers )
		{
			$i--;
			$j--;

			# do not necessarily decrease abundance of both borders
			# --> check which border has a higher abundance, and revert the decrease in abundance for the other one
			if ( $i >= 0 and $j >= 0 )
			{
				if ( $abundances_first[$i] > $abundances_last[$j] )
				{
					$j++;
				}
				elsif ( $abundances_first[$i] < $abundances_last[$j] )
				{
					$i++;
				}
			}

			# these are already the least abundant contig borders --> stop while loop
			elsif ( $i < 0 and $j < 0 )
			{
				$next_less_abundant = 0;    # last;
			}

			# least abundant upstream border, but still a less abundant downstream border available
			elsif ( $i < 0 )
			{
				$i = 0;                     # keep i, decrease j
			}

			# vice versa
			elsif ( $j < 0 )
			{
				$j = 0;                     # keep j, decrease i
			}
		}
		else
		{
			$next_less_abundant = 0;        # last;
		}
	}

	# print most abundant prediction(s) --> table head
	say "gene names   promoter numbers   length   motif scores*   abundances**" . "\n"
	  . "---------------------------------------------------------------------";

	# print actual prediction(s) --> table body
	say "$_->[0] - $_->[1]   $_->[2] - $_->[3]   $_->[4]   $_->[5]/$_->[6]   $_->[7]/$_->[8] of $_->[9]" for @cluster_predictions_most_abundant;

	# print small explanation
	say "\n" . "*  the lower, the better"                               if $fimo;
	say "\n" . "*  not applicable, binding sites given by SiTaR input " if $sitar;
	say "** the higher, the better";

	# print warnings, if present
	print "@warnings" if @warnings;
}

####################
# Faulty promoters #
####################

sub faulty_promoters { }

# check if there are faulty promoters within cluster predictions
# IF there is at least one faulty promoter
if (@faulty_promoters_short)
{
	say "\n(9) Checking if FAULTY PROMOTERS interfere with cluster prediction(s) …";

	my @unique_first = uniq map { $_->[0] } @cluster_predictions_most_abundant;
	my @unique_last  = uniq map { $_->[1] } @cluster_predictions_most_abundant;

	my @faulty;        # position of ALL faulty features/promoters
	my $very_first;    # position of left border of FIRST prediction in the genome
	my $very_last;     # position of right border of LAST prediction in the genome

	for my $feat ( 0 .. $#features )
	{
		if ( any { /(^\Q$features[$feat]{ID}\E$)|(^\Q$features[$feat]{ID}+\E)|(\Q+$features[$feat]{ID}\E$)/ } @faulty_promoters_short )
		{
			push( @faulty, $feat );
		}
	}

	for my $feat ( 0 .. $#features )
	{
		if ( any { /(^\Q$features[$feat]{ID}\E$)|(^\Q$features[$feat]{ID}+\E)|(\Q+$features[$feat]{ID}\E$)/ } @unique_first )
		{
			$very_first = $feat;
			last;
		}
	}

	for my $feat ( reverse 0 .. $#features )
	{
		if ( any { /(^\Q$features[$feat]{ID}\E$)|(^\Q$features[$feat]{ID}+\E)|(\Q+$features[$feat]{ID}\E$)/ } @unique_last )
		{
			$very_last = $feat;
			last;
		}
	}

	if ( any { $very_first <= $_ and $_ <= $very_last } @faulty )
	{
		say "[WARNING] At least one faulty promoter is located within a cluster prediction!";
		say "";
		print autoformat "\"Faulty\" promoters are promoter sequences which have been excluded from the analysis (ignored). "
		  . "Ignored promoters inside a predicted cluster could have influenced the outcome of the prediction. "
		  . "There are various reasons for such an exclusion, mostly overlapping gene annotations. "
		  . "Check the annotation of the corresponding genes, especially start/stop positions, to avoid faulty promoters.";
	}
	else
	{
		say "no :)";
	}
}

########
# Done #
########

if ($verbose)
{
	my @runtime = localtime( time - $^T );    # @runtime = (seconds, minutes);
	say "\n" . "runtime: $runtime[1] minute(s) and $runtime[0] second(s)";
}
