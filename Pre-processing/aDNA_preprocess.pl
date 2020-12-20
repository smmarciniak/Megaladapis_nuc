#Takes two *fastq.gz files, trims and merges them using MergeReadsFastQ_cc and SplitMerged2Bwa both created by Martin Kircher (2009)
#perl script created by Perry lab members (PSU)
#Used to generate prelim seq stats

#!/usr/bin/perl

use strict;
use warnings;

#Takes 2 short read fastq files, gzipped or not, gets them ready for merging, and writes a shell script for qsub. Assumes demultiplexing is done.

my $r1 = $ARGV[0];	#r1.fastq(.gz)
my $r2 = $ARGV[1];	#r2.fastq(.gz)
my $reads = $ARGV[2];	#reads to prep for each queued job
my $out = $ARGV[3];	#output prefix

unless (@ARGV==4) {
	die "\nRequired: Read1.fastq(.gz), Read2.fastq(.gz), read n per job, output prefix\n\n";
}

unless (-f $r1) {
	die "\nCouldn't find $r1\n\n";
}

unless (-f $r2) {
        die "\nCouldn't find $r2\n\n";
}

my $gz = 0;
if ($r1 =~ /.gz$/ and $r2 =~ /.gz$/) {
	$gz = 1;
}
elsif ($r1 =~ /.gz$/ or $r2 =~ /.gz$/) {
	die "\nIt appears that one read file is compressed and the other is not. Please standardize\n\n";
}

if ($gz==1) {
	open(R1, sprintf("zcat %s |", $r1));
	open(R2, sprintf("zcat %s |", $r2));
}

else {
	open R1, "$r1";
	open R2, "$r2";
}

my $NR=0;
my @lens1;
my @lens2;
my @lenq1;
my @lenq2;
my $splitn = 0;
my $stem = ("0" x (4-(length $splitn))).$splitn;
open OUT, ">$out.split.$stem.fastq";
while (!eof(R1) and !eof(R2)) {
	$NR++;
	if ($NR>1 and ($NR%(4*$reads))==1) {
		close OUT;
		$splitn++;
		$stem = ("0" x (4-(length $splitn))).$splitn;
		open OUT, ">$out.split.$stem.fastq";
	}
	my $w1 = <R1>;
	my $w2 = <R2>;
        chomp $w1;
        chomp $w2;
	if ($NR%4==1) {
		my $header = (split(/\s/, $w1))[0];
		print OUT $header."\n";
	}
	elsif ($NR%4==2) {
		push (@lens1, length $w1);
		if (@lens1 > 2) {
			shift @lens1;
			unless ($lens1[0]==$lens1[1]) {
				die "\nUneven sequence lines in $r1!\n\n";
			}
		}
		push (@lens2, length $w2);
                if (@lens2 > 2) {
                        shift @lens2;
                        unless ($lens1[0]==$lens2[1]) {
                                die "\nUneven sequence lines in $r2!\n\n";
                        }
                }
		print OUT $w1.$w2."\n";
	}
	elsif ($NR%4==3) {
		my $qualhead = (split(/\s/, $w1))[0];
		print OUT "$qualhead\n";
	}
	else {
                push (@lenq1, length $w1);
                if (@lenq1 > 2) {
                        shift @lenq1;
                        unless ($lenq1[0]==$lenq1[1]) {
                                die "\nUneven quality lines in $r1!\n\n";
                        }
                }
                push (@lenq2, length $w2);
                if (@lenq2 > 2) {
                        shift @lenq2;
                        unless ($lenq1[0]==$lenq2[1]) {
                                die "\nUneven quality lines in $r2!\n\n";
                        }
                }
                print OUT $w1.$w2."\n";
	}
}
close OUT;
open FINAL, ">$out.clean.pl";
print FINAL "#!/usr/bin/perl\n

`cat $out.split.????.r1.fastq > $out.r1.fastq`;
`cat $out.split.????.r2.fastq > $out.r2.fastq`;
`cat $out.split.????.SR.fastq > $out.SR.fastq`;
`rm $out.split.????.r1.fastq`;
`rm $out.split.????.r2.fastq`;
`rm $out.split.????.SR.fastq`;
`rm $out.split.????.fastq`;
`rm $out.split.????.merge.qsub`;
`rm $out.split.????.merge.qsub.e*`;
`rm $out.split.????.merge.qsub.o*`;
`rm $out.clean.pl`;
`rm $out.masterQ.sh`;
";

open MASTERQ, ">$out.masterQ.sh";
my $revstart = $lens1[0]+1;
for my $iter (0 .. $splitn) {
	$stem = ("0" x (4-(length $iter))).$iter;
	my $hours = int(4*($reads/1000000));
	open OUT, ">$out.split.$stem.merge.qsub";
	print OUT "#PBS -l walltime=$hours:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=1gb
cd \$PBS_O_WORKDIR
LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/gpfs/group/ghp3/default/libraries_modules/merge/lib/
export LD_LIBRARY_PATH
module load gcc
/gpfs/group/ghp3/default/scripts/MergeReadsFastQ_cc -p -r $revstart -f 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -s 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA' < $out.split.$stem.fastq \\
| /gpfs/group/ghp3/default/scripts/SplitMerged2Bwa.py -p -f $out.split.$stem.r1.fastq -r $out.split.$stem.r2.fastq -m $out.split.$stem.SR.fastq \n";
	close OUT;
	print MASTERQ "qsub $out.split.$stem.merge.qsub\n";
}
`sh $out.masterQ.sh`;
#print "\nFinished. Use \"sh $out.masterQ.sh\" to send your merging jobs to the cluster.\n";
#print "When all merging is finished, use \"perl $out.clean.pl\" to clean up intermediates and combine merged files\n\n";
