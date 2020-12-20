#Removing duplicates from paired *fastq.gz files

#!/usr/bin/env perl

use strict;
use warnings;

my $in = $ARGV[0];
my $out = $ARGV[1];
my %ids;
my %seqs;
my $NR=0;
open IN, "$in";
open OUT, ">$out";
my $run = 0;
while (my $line = <IN>) {
	$NR++;
	chomp $line;
	if ($NR%4==1) {
		my $header = $line;
		chomp $header;
		my $seq = <IN>;
		chomp $seq;
		my $rev = reverse $seq;
		$rev =~ tr/ACGTacgt/TGCAtgca/;
		$seq = (sort($rev, $seq))[0];
		$NR++;
		unless (exists ($seqs{$seq})) {
			$seqs{$seq} = 0;
			$run++;
			print OUT $header."\n".$seq."\n";
			my $i = <IN>;
			chomp $i;
			$NR++;
			print OUT $i."\n";
                        $i = <IN>;
                        chomp $i;
                        print OUT $i."\n";
			$NR++;
			if ($NR%400000==0) {
				my $reads = $NR/4;
				print "$reads reads analyzed, $run kept\n";
			}
		}
	}
}
	
