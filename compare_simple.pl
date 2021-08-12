#!/usr/bin/perl
# All-vs-all comparison between "simple" files (chr, pos, genotype)
# outputs a matrix of fraction discordant bases between two samples.

#Copyright 2021 Jamie K. Teer
#
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


use warnings;
use strict;
use Getopt::Long;
use File::Basename;

my %Opt;
process_commandline();

my %CHROM = ( "X" => 23, "Y" => 24, "MT" => 25, "M" => 25 );
my %regions;    # {chr}->[start, end], stored as 1-based
my @sample_names;  # names of each sample
my @sample_gens;   # [sample]->[genotypes (or NA) at each position in bed]
my @output;        # [sample1]->[sample2]  2d matrix of comparison results (fraction discordant)
my $raw;           # Filehandle for raw output

# Use NULL for missing, as it has value 0 for later bitwise operations
# Use capital letters for genotypes. They all have bit 7 set, and so an AND will always return a nonzero value
my %geno_map = (
            $Opt{missing} => "\0",
            "AA" => "A",
            "AC" => "M",
            "CA" => "M",
            "AG" => "R",
            "GA" => "R",
            "AT" => "W",
            "TA" => "W",
            "CC" => "C",
            "CG" => "S",
            "GC" => "S",
            "CT" => "Y",
            "TC" => "Y",
            "GG" => "G",
            "GT" => "K",
            "TG" => "K",
            "TT" => "T",
          );


if ($Opt{raw_output}) {
    open $raw, ">$Opt{raw_output}" or die $!;
}

#################
## Load bed regions to memory
#################
my ($last_chr, $last_pos);
open BED, $Opt{bed} or die $!;
while (<BED>) {
    chomp;
    my @l = split /\s+/, $_;
    $l[0] =~ s/^chr//;
    $l[0] = ($CHROM{$l[0]}) ? $CHROM{$l[0]} : $l[0];
    if (!$last_chr || $l[0] != $last_chr) {
        $last_pos = 0;
        $last_chr = $l[0];
    }

    if ($l[1] < $last_pos ) {
        die "Bedfile out of order or not merged: $l[0] $l[1] $l[2] Exiting.\n";
    }
    else {
        push @{$regions{$l[0]}}, [$l[1]+1, $l[2]]; # storing 1-based coords
        $last_pos = $l[2];
    }
}
close BED or die $!;

#################
## Load sample data
#################

my $start_t_sample = time();

SAMPLE: for my $s_file (@ARGV) {
    if (! -s $s_file) {
        warn "File $s_file doesn't exist; skipping.\n";
        next SAMPLE;
    }
    
    my $name = fileparse($s_file, ".smp.snv");
    push @sample_names, $name;
    open my $SF, $s_file or die $!;
    my @this_sample_genotypes;
    my $line_ar = get_line($SF);
    if (! defined $line_ar) {
        die "No lines for sample $s_file?\n";
    }

    CHR: for my $c (sort {$a <=> $b} keys %regions) {
        REG: for my $r ( @{$regions{$c}} ) {
            POS: for my $p ($r->[0]..$r->[1]) {

                # Sample file has ended - fill in missing value.
                if (! defined $line_ar) {
                    push @this_sample_genotypes, $geno_map{$Opt{missing}};
                    next POS;
                }
                
                # Must test for chr, position at each position, so empty value can be added
                if ($line_ar->[0] < $c || ($line_ar->[0] == $c && $line_ar->[1] < $p) ) {
                    # Sample file is behind region. Advance sample file. This shouldn't happen.
                    warn "Sample file has bases not in region? $line_ar->[0]:$line_ar->[1]  $c:$p\n";
                    $line_ar = get_line($SF);
                    redo POS;
                }
                elsif ($line_ar->[0] > $c || ($line_ar->[0] == $c && $line_ar->[1] > $p) ) {
                    # Sample file ahead of region. Add empty value and advance region.
                    push @this_sample_genotypes, $geno_map{$Opt{missing}};
                    next POS;
                }
                elsif ($line_ar->[0] == $c && $line_ar->[1] == $p) {
                    # Sample file aligned with region.  Add value from file.
                    push @this_sample_genotypes, $geno_map{$line_ar->[2]};
                    $line_ar = get_line($SF);
                }
                else {
                    warn "Shouldn't get here: chrom: $c _ $line_ar->[0]  pos: $p _ $line_ar->[1]\n";
                }
            }
        }
    }
    push @sample_gens, join("", @this_sample_genotypes);
    if ($Opt{raw_output}) {
        print $raw (join "\t", $name, @this_sample_genotypes);
        print $raw "\n";
    }
    close $SF or die $!;
}
if ($Opt{raw_output}) {
    close $raw or die $!;
}

my (%s_counts, $pos_count);
for (my $i=0; $i<=$#sample_names; $i++) {
    #print "$sample_names[$i]\t", scalar @{$sample_gens[$i]}, "\n";
    $s_counts{length($sample_gens[$i])}++;
}
if (length keys %s_counts != 1) {
    die "Sample loading error - number of positions different.\n";
}
else {
    $pos_count = (keys %s_counts)[0];
}

my $end_t_sample = time();
printf STDERR ("%s\t%d\n", "Load time:", ($end_t_sample - $start_t_sample));


#################
# All vs all comparison
#################
my $start_t = time();
for (my $i=0; $i <= $#sample_names; $i++) {

    my $a = $sample_gens[$i];
    #my $start_t_row = time();

    for (my $j=$i; $j <= $#sample_names; $j++) {
        
        my $b = $sample_gens[$j];

        #########
        # Use bitwise operations for FAST string comparisons

        # XOR between strings results in 0 at matching chars. We then count NULL chars (value of 0).
        my $match_count = ($a ^ $b) =~ tr/\0//;

        # Missing values are encoded as NULL. Two missing values will match, but we don't want to count them.
        #  So, count number of 0 values using OR.
        my $match_missing = ($a | $b) =~ tr/\0//;

        # We need the total number of bases minus the number of bases that were missing in one or both samples.
        #   To determine the number of missing chars (encoded as NULL), use AND. The characters encoding 
        #   genotypes all have bit 6 set, and so will all AND to a non-zero value.
        my $tot_zero = ($a & $b) =~ tr/\0//;

        # Usable bases is equal to the initial total minus the number of missing bases.
        my $tot = length($a) - $tot_zero;

        # The number of discrepant bases is equal to the total usable bases minus the true matching.
        #   True matching are total matches minus missing matches.
        my $dis = $tot - ($match_count - $match_missing);

        # finally, calculate a rate
        $output[$i][$j] = sprintf("%0.3f", ($dis / $tot));
    }
    #my $end_t_row = time();
    #printf STDERR ("%s\t%d\n", "Row time:", ($end_t_row - $start_t_row));

}
my $end_t = time();
printf STDERR ("%s\t%d\n", "Comparison time:", ($end_t - $start_t));

#################
# Print output
#################
print  join( "\t", "", @sample_names), "\n";
foreach my $l (@output) {
    {
        no warnings; # we know there are missing values as we only calculated one-way comparison.
        print join( "\t", shift @sample_names, @$l), "\n";
    }
}

print "Samples analyzed: ", scalar @sample_gens,"\n";


sub get_line {
    my $fh = shift;
    if (defined (my $l = <$fh>)) {
        chomp $l;
        my @l = split /\t/, $l;
        $l[0] =~ s/^chr//;
        $l[0] = ($CHROM{$l[0]}) ? $CHROM{$l[0]} : $l[0];
        return \@l;
    }
    else {
        return;
    }
}


sub process_commandline {
    #Defaults
    %Opt = ( "missing" => "NA",
           );

    GetOptions(\%Opt, qw(
                          bed=s
                          missing=s
                          raw_output=s
              ));

    usage() if (! $Opt{bed});

    if (@ARGV < 2) {
        warn "Can only compare two or more samples!\n";
        usage();
    }

}

sub usage {
    print STDERR << "EOF";
usage: $0 --bed <bedfile> <file_1> <file_2> ... <file_n>

--bed:  A standard bedfile with the regions to be compared. Every one of the positions in this
        file will be examined.
        Although chromosomes need not be sorted, the regions within a chromosome MUST BE SORTED!!!
        File must be merged !!!
        !!! Note that this only works with hg19 style chroms (1-22, X, Y, MT) !!!

<file_n1>: A file in "simple" format (tab-delim): chr pos genotype

Optional:

--missing:  String to use for missing values. Default = "NA";

--raw_output: File in which to output the genotype for each sample (row) at each position (col)
EOF
exit;
}
