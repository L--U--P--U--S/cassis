#!/usr/bin/perl

#####################################################
# Run MEME iteratively                              #
# on different promoter sets around the anchor gene #
# to obtain an optimal motif                        #
#####################################################

# Copyright (C) 2015 Leibniz Institute for Natural Product Research and
# Infection Biology -- Hans-Knoell-Institute (HKI)

use 5.14.2;
use strict;
use warnings;

use XML::LibXML;

##############
# Parameters #
##############

my @plus_minus;
for ( my $plus = 0 ; $plus <= 15 ; $plus += 1 )
{
	for ( my $minus = 0 ; $minus <= 15 ; $minus += 1 )
	{
		push( @plus_minus, { plus => $plus, minus => $minus } ) if ( $plus + $minus + 1 >= 4 );
	}
}

my @meme_promoters;      # { ID, sequence }
my @predicted_motifs;    # { motif, score }

###############
# Main method #
###############

sub meme
{
	my ( $dir, $promoter_sequences_file, $cluster_name, $backbone_id, $backbone_contig, $verbose, $meme_parameters_ref, $num_cpus ) = @_;

	# STEP 1: read all promoter sequences from file
	sub read_promoter_sequences { }

	my $sequence_io = Bio::SeqIO->new( -file => $promoter_sequences_file );
	while ( my $sequence = $sequence_io->next_seq() )
	{
		$sequence->display_id() =~ m/__contig_(.+)$/;
		push( @meme_promoters, { ID => $sequence->display_id(), seq => $sequence->seq() } ) if ( $1 eq $backbone_contig );
	}
	die "Anchor's locus does not appear in all_promoter_sequences file.\nStopped"
	  if ( none { /$backbone_id/ } map { $_->{ID} } @meme_promoters );

	# STEP 2: copy promoter sets to seperate files, to be used as MEME input
	# no matter if motif has been predicted previously or not (it's fast, don't worry)
	sub copy_promoter_sets { }

	my @read_from_previous;
	my @set_large_enough;
	my %promoter_file_sizes;
	for ( 0 .. $#plus_minus )
	{
		my $plus     = $plus_minus[$_]{plus};
		my $minus    = $plus_minus[$_]{minus};
		my $meme_dir = "$dir$cluster_name/meme/" . "+" . $plus . "_-" . $minus . "/";
		make_path($meme_dir);    # create directory if it does not exist

		# motif predicted previously?
		if ( -s $meme_dir . "promoters.fasta" and -s $meme_dir . "meme.xml" )
		{
			push( @read_from_previous, 1 );
		}
		else
		{
			push( @read_from_previous, 0 );
		}

		push( @set_large_enough, 1 );    # first, assume we will have at least 4 promoters in the set

		my $meme_promoters_file = $meme_dir . "promoters.fasta";
		open( my $meme_promoters_io, ">", $meme_promoters_file ) or die "Cannot write to promoters file \"$meme_promoters_file\".\n", $!;

		for ( 0 .. $#meme_promoters )
		{
			# look for promoter of anchor gene --> center of "plus_minus"
			# it IS important to check for "$backbone_id+" and "$backbone_id__"
			if ( index( $meme_promoters[$_]{ID}, $backbone_id . "+" ) != -1 or index( $meme_promoters[$_]{ID}, $backbone_id . "__" ) != -1 )
			{
				for ( my $i = $_ - $minus ; $i <= $_ + $plus ; $i++ )
				{
					if ( $i < 0 )    # locus near beginning of contig
					{
						printf( "%s   %7s   %s\n", "promoter set", "+" . $plus . "_-" . $minus, "exceeds upstream contig border" ) if ($verbose);
						$i = 0;      # set next locus to first one of the contig

						# after modifying $i, check again if we have at least 4 promoter sequences in the set
						if ( $plus + $i + 1 < 4 )    # compare to "$plus + $minus + 1 >= 4" from @plus_minus init
						{
							$set_large_enough[-1] = 0;
						}
					}
					elsif ( $i > $#meme_promoters )    # locus near end of contig
					{
						printf( "%s   %7s   %s\n", "promoter set", "+" . $plus . "_-" . $minus, "exceeds downstream contig border" ) if ($verbose);
						last;    # we are at the end of the contig --> nothing more to do here --> jump to last iteration of the for-loop
					}

					( my $seq_fasta = $meme_promoters[$i]{seq} . "\n\n" ) =~ s/(.{60})/$1\n/g;

					if ( $i == $_ )
					{
						print $meme_promoters_io ">" . $meme_promoters[$i]{ID} . "__ANCHOR\n" . $seq_fasta;
					}
					else
					{
						print $meme_promoters_io ">" . $meme_promoters[$i]{ID} . "\n" . $seq_fasta;
					}
				}
			}
		}
		close $meme_promoters_io;

		# promoters.fasta file with same size already present (in different subdir)?
		my $file_size = -s $meme_promoters_file;
		if ( $promoter_file_sizes{$file_size} )
		{
			# promoters.fasta file with same size AND same content already present!
			# remark: size checked before content, because it's much less computational expensive
			if ( compare( $meme_promoters_file, $promoter_file_sizes{$file_size}{file} ) == 0 )
			{
				# mark promoter sets, which are duplicates (same set of sequences, due to proximity to contig border) of previous sets
				# --> already in use for MEME input
				printf( "%s   %7s   %s\n", "promoter set", "+" . $plus . "_-" . $minus, "duplicate of " . $promoter_file_sizes{$file_size}{motif} )
				  if ($verbose);
				$plus_minus[$_] = "";
			}
			else
			{
				# same size (by chance) but different content
				$promoter_file_sizes{$file_size} = { motif => "+" . $plus . "_-" . $minus, file => $meme_promoters_file };
			}
		}
		else
		{
			# promoters.fasta with (up to now) unique file size
			$promoter_file_sizes{$file_size} = { motif => "+" . $plus . "_-" . $minus, file => $meme_promoters_file };
		}
	}

	# STEP 3: run MEME on each promoter set
	sub run_meme { }

	my $pm;    # max processes (== number of CPUs) for parallel MEME instances, if at least 2 CPUs, no forks if less CPUs
	( $num_cpus >= 2 ) ? ( $pm = new Parallel::ForkManager($num_cpus) ) : ( $pm = new Parallel::ForkManager(0) );

	# progress bar (sort of), but updates -- thanks to \r -- the current line, instead of printing to the next one
	my $progress = 1;
	my $total_plus_minus = grep { $_ ne "" } @plus_minus;
	$pm->run_on_finish( sub { print "\rProgress: " . $progress++ . "/$total_plus_minus" if ($verbose); } );

	print "Progress: 0/$total_plus_minus" if ($verbose);
	for ( 0 .. $#plus_minus )
	{
		next if ( !$plus_minus[$_] );          # skip promoter sets which have been marked as duplicates before
		next if ( !$set_large_enough[$_] );    # skip promoter sets which do not have at least 4 promoters

		my $plus     = $plus_minus[$_]{plus};
		my $minus    = $plus_minus[$_]{minus};
		my $meme_dir = "$dir$cluster_name/meme/" . "+" . $plus . "_-" . $minus . "/";

		$pm->start and next;                                                            # fork
		if ( !$read_from_previous[$_] )                                                 # motif not predicted previously? --> run MEME
		{
			my @args = ( $meme_dir . "promoters.fasta", "-oc", $meme_dir, @$meme_parameters_ref );
			system( "meme", @args );
		}
		$pm->finish;                                                                    # exit child process
	}
	$pm->wait_all_children;                                                             # wait for all child processes to finish
	say "\n" if ($verbose);

	# STEP 4: analyse MEME ouput
	# no matter if motif has been predicted previously or not (it's fast, don't worry) - we need the information anyway
	sub analyse_meme_output { }
	printf( "%7s   %s\n", "motif", "e-value (motif score)" ) if ($verbose);

	for ( 0 .. $#plus_minus )
	{
		next if ( !$plus_minus[$_] );
		next if ( !$set_large_enough[$_] );

		my $promoter_set = $_;
		my $plus         = $plus_minus[$_]{plus};
		my $minus        = $plus_minus[$_]{minus};

		my $meme_xml_file = "$dir$cluster_name/meme/" . "+" . $plus . "_-" . $minus . "/" . "meme.xml";
		my $meme_bs_file  = "$dir$cluster_name/meme/" . "+" . $plus . "_-" . $minus . "/" . "binding_sites.fasta";

		my $anchor_sequence_id;
		my $motif_score;
		my @sequences;
		my $check_result;

		# load MEME xml output file
		my $xml = XML::LibXML->new()->load_xml( location => $meme_xml_file );
		my $stopping_reason = $xml->findnodes('//model/reason_for_stopping')->to_literal();

		# NO motif found for given e-value cutoff :(
		if ( index( $stopping_reason, "Stopped because motif E-value > " ) != -1 )
		{
			$check_result = sprintf( "   %-24s", "exceeds e-value cut-off" );
		}

		# motif(s) found :)
		elsif ( $stopping_reason =~ m/Stopped because(.+)motifs(.+)[found|reached]/ )    # allow some flexibility for different MEME versions
		{
			# save id of input sequence of anchor gene (promoter)
			for ( $xml->findnodes('//sequence') )
			{
				$anchor_sequence_id = $_->findvalue('./@id') if ( index( $_->findvalue('./@name'), "ANCHOR" ) != -1 );
			}

			# save motif score (e-value)
			$motif_score =
			  $xml->findnodes('//motifs/motif/@e_value')->to_literal()->value()
			  ;    # need to use value(), otherwise sorting the e-values via "<=>" will complain

			# save actual motif sequences
			for my $site ( $xml->findnodes('//motifs/motif/contributing_sites/contributing_site/site') )
			{
				my $sequence = $site->findnodes('./letter_ref/@letter_id')->to_literal();
				$sequence =~ s/letter_//g if ( index( $sequence, "letter_" ) == 0 );
				push( @sequences, $sequence );
			}

			# only accept motifs which occur (among others) in the anchor genes promoter
			if ($anchor_sequence_id)
			{
				$check_result = sprintf( "   %-24s", $motif_score );
				push( @predicted_motifs, { motif => "+$plus" . "_" . "-$minus", score => $motif_score } );

				open( my $bs_io, ">", $meme_bs_file ) or die "Cannot write to binding sites file \"$meme_bs_file\".\n", $!;
				say $bs_io ">" . "+$plus" . "_" . "-$minus";
				print $bs_io join( "\n", @sequences );
				close $bs_io;
			}
			else
			{
				$check_result = sprintf( "   %-24s", "not in anchor's promoter" );
			}
		}

		# unexpected reason for stopping :$
		else
		{
			die "[ERROR] Reason for stopping (MEME): \"$stopping_reason\".\nStopped";
		}

		die "[ERROR] Couldn't parse the MEME xml output file \"$meme_xml_file\" properly.\nStopped"
		  if !$check_result;

		# print result of the analysis of the MEME output file
		if ( $read_from_previous[$promoter_set] )
		{
			printf( "%7s%s   *\n", "+" . $plus . "_-" . $minus, $check_result ) if ($verbose);
		}
		else
		{
			printf( "%7s%s\n", "+" . $plus . "_-" . $minus, $check_result ) if ($verbose);
		}
	}

	say "(* No new prediction, using existing MEME files.)" if ( sum(@read_from_previous) and $verbose );

	return @predicted_motifs;
}

return 1;
