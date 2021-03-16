#Script for removal of duplicates in fastq files after lastZ alignment
#Created by GH Perry lab members

#!/usr/bin/env perl

use strict;
use warnings;

my $arg = join " ", @ARGV;
chomp $arg;

unless ($arg) {
	die "\nUsage: fastq_filter.pl [options] infile.fq \nFor more on usage, try fastq_filter.pl --help\n\n";
}


my $in = (split("\ ", $arg))[-1];
$arg = substr $arg, 0, ((length $arg)-(length $in));
unless (-f $in) {
        die "\nCouldn't find $in, check filename\n\n";
}

my $cg = 1;
my $minlen = 0;
my $maxlen;
my $rmdup = 0;
my $lowcomp = 0;
my $o;
my $p5 = 0;
my $p7 = 0;
my $cgfilt = "N/A";
my $dupfilt = "N/A";
my $lenfilt = "N/A";
my $maxfilt = "N/A";

my @args = split "-", $arg;
shift @args;
for my $i (@args) {
        my $flag = (split(/\s/, $i))[0];
        my $value = (split(/\s/, $i))[1];
        unless ($flag eq "cg" or $flag eq "minlen" or $flag eq "rmdup" or $flag eq "lowcomp" or $flag eq "o" or $flag eq "maxlen" or $flag eq "p5" or $flag eq "p7") {
		die "\nNO!\n\n";
	}
	if ($flag eq "cg") {
		$cg = $value;
		$cgfilt = 0;
	}
	if ($flag eq "minlen") {
		$minlen = $value;
		$lenfilt = 0;
	}
	if ($flag eq "maxlen") {
		$maxlen = $value;
		$maxfilt = 0;
	}
	if ($flag eq "rmdup") {
		$rmdup = 1;
		$dupfilt = 0;
	}
	if ($flag eq "lowcomp") {
		$lowcomp=1;
	}
	if ($flag eq "p5") {
		$p5 = 1;
	}
	if ($flag eq "p7") {
		$p7 = 1;
	}
	if ($flag eq "o") {
		$o = $value;
	}
}
 
unless ($o) {
	die "\nPlease specify an output prefix\n\n";
}

if ($in=~/\.gz$/) {
        open(IN, sprintf("zcat %s |", $in));
}
else {
        open IN, $in;
}

open LOG, ">$o.filter.log";
open OUT, ">$o.fastq";

my $NR=0;
my $keep=0;
my $check=0;
my %dups;
while (my $line = <IN>) {
	chomp $line;
	$NR++;
	if ($NR%4==1) {
		$check=0;
		unless ($line =~ /^@/) {
			print LOG "File error at line $NR of $in, aborting.\n";
			die "\nOH NO! $in appears to be truncated or squelched at line $NR\n\n";
		}
		my $header = $line;
		my $seq = <IN>;
		chomp $seq;
		$NR++;
=pod
		if ($p5==1) {
		}
		if ($p7==1) {
		}
=cut
		if ($minlen > 0) {
			if ((length $seq) < $minlen) {
				$check=1;
				$lenfilt++;
			}
		}
                if ($maxlen) {
			if ((length $seq) > $maxlen) {
                                $check=1;
                                $maxfilt++;
                        }
                }
		if ($rmdup == 1) {
#			my $rev = reverse $seq;
#			$rev =~ tr/ACGTacgt/TGCAtgca/;
#			my $query = (sort($rev, $seq))[0];
			my $query = $seq;
			if (exists ($dups{$query})) {
				$check=1;
				$dupfilt++;
			}
			else {
				$dups{$query} = 0;
			}
		}
		if ($cg < 1) {
			my @GC = $seq =~ /[GgCc]/g;
			my $thresh = @GC/(length $seq);
			if ($thresh > $cg) {
				$check=1;
				$cgfilt++;
			}
		}
		if ($check==0) {
			$keep++;
			print OUT "$header\n$seq\n";
			my $i = <IN>;
			chomp $i;
			$NR++;
			print OUT $i."\n";
			$i = <IN>;
			chomp $i;
			$NR++;
			print OUT $i."\n";
		}
	}
	if ($NR%400000==0) {
		my $reads=$NR/4;
		print "$reads reads scanned, $keep kept.\n";
		print LOG "$reads reads scanned, $keep kept.\n";
	}
}

my $reads = $NR/4;
my $readperc = $keep/$reads;

print "\nFinished fastq filtering. Scanned $reads total reads and kept $keep ($readperc\%)\n";
print LOG "\nFinished fastq filtering. Scanned $reads total reads and kept $keep ($readperc\%)\n"; 
if ($minlen>0) {
	print "Minimum length filtered: ".(substr 100*$lenfilt/$reads, 0, 5)."\%\n";
	print LOG "Minimum length filtered: ".(substr 100*$lenfilt/$reads, 0, 5)."\%\n";
}
if ($maxlen) {
	print "Maximum length filtered: ".(substr 100*$maxfilt/$reads, 0, 5)."\%\n";
	print LOG "Maximum length filtered: ".(substr 100*$maxfilt/$reads, 0, 5)."\%\n";
}
if ($cg<1) {
	print "CG filtered: ".(substr 100*$cgfilt/$reads, 0, 5)."\%\n";
	print LOG "CG filtered: ".(substr 100*$cgfilt/$reads, 0, 5)."\%\n";
}
if ($rmdup==1) {
	print "Duplicate filtered: ".(substr 100*$dupfilt/$reads, 0, 5)."\%\n";
	print LOG "Duplicate filtered: ".(substr 100*$dupfilt/$reads, 0, 5)."\%\n";
}
