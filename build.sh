#!/bin/bash

# Perl modules
# sudo apt-get install --no-upgrade libperl-dev cpanminus
# sudo cpanm --skip-satisfied PAR::Packer LWP::UserAgent Bio::Seq Data::Dumper File::Basename File::Copy::Recursive File::Which Getopt::Long IPC::System::Simple List::AllUtils Parallel::ForkManager Roman Scalar::Util Sort::Naturally Sys::CpuAffinity Text::CSV Text::CSV_XS Text::CSV_PP

# build CASSIS binary
sed -i 's/usage: cassis\.pl/usage: cassis/' cassis.pl
sed -i 's/cassis\.pl --dir/cassis --dir/' cassis.pl
perl cassis.pl > /dev/null || exit
pp -M Text::CSV_XS -M Text::CSV_PP -o cassis cassis.pl cassis_meme.pl || exit

# pack CASSIS binary
rm "CASSIS Linux 64bit.7z"
7z a "CASSIS Linux 64bit.7z" cassis README COPYING

# clean up
rm cassis
sed -i 's/usage: cassis/usage: cassis\.pl/' cassis.pl
sed -i 's/cassis --dir/cassis\.pl --dir/' cassis.pl
