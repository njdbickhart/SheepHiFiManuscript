#!/usr/bin/perl
# this script is designed to downsample several Illumina fastq files using a progressive downsampling algorithm
# The goal is to prepare several subset fastq files for MASH profile estimation.
# Downsample thresholds will be in 5% intervals

use strict;
use FileHandle;

my $usage = "perl $0 <fasta file 1> ...\n";

chomp(@ARGV);
unless(scalar(@ARGV) >= 1){
	print $usage;
	exit;
}

open(my $LOG, "> downsampleFasta.log");

# Calculate read numbers
my $totalReads = 0;
foreach my $file (@ARGV){
	print {$LOG} "Working on file: $file\n";
	print "Working on file: $file\n";
	$totalReads += GetFastaReadCount($file);
}

print {$LOG} "Total reads:\t$totalReads\n";
print "Total reads:\t$totalReads\n";

# Create output streams
my @ratios = ("0.90", "0.80", "0.70", "0.60", "0.50", 
	"0.40", "0.30", "0.20", "0.10");
my @fhs; my %counts;
foreach my $r (@ratios){
	my $t = $r;
	$t =~ s/^0\.//; 
	push(@fhs, FileHandle->new("> downsample.$t.fasta"));
}

# Start processing the files
foreach my $file (@ARGV){
	print {$LOG} "Sketching from file: $file\n";
	print "Sketching from file: $file\n";
	open(my $IN, "< $file");
	my $head = '';
	my $seq = '';
	while(my $read = <$IN>){
		if($read =~ /^>/){
			chomp($read);
			if($seq eq ''){
				# First entry!
				$head = $read;
				next;
			}
			chomp($seq);
			
			# Random counter out here to reduce attrition from sequential tests
			my $randm = rand();
			recursiveSelection($randm, -1, \@ratios, \@fhs, $head, $seq);
		
			for(my $x = 0; $x < scalar(@ratios); $x++){
				if($randm < $ratios[$x]){
					# Read counter update
					$counts{$ratios[$x]} += 1;
				}
			}
			$head = $read;
			$seq = ""; 
		}else{
			$seq .= $read;
		}
	}
	close($IN);

	# Last read of the file 
	if(length($seq) > 1){
		my $randm = rand();
                recursiveSelection($randm, -1, \@ratios, \@fhs, $head, $seq);
		
		for(my $x = 0; $x < scalar(@ratios); $x++){
                	if($randm < $ratios[$x]){	
				$counts{$ratios[$x]} += 1;
			}
		}
	}

}

# Log the read counters
foreach my $keys (sort {$b <=> $a} keys(%counts)){
	my $trueRatio = $counts{$keys} / $totalReads;
	print {$LOG} "$keys\t$counts{$keys}\t$trueRatio\n";
	print "$keys\t$counts{$keys}\t$trueRatio\n";
}
close $LOG;

exit;

# Returns the smallest ratio that this read passes
sub recursiveSelection{
	my ($randm, $startIdx, $ratios, $fhs, $head, $seq) = @_;
	if($startIdx + 1 >= scalar(@{$ratios})){
		return $randm; # We reached the end of the ratio list!
	}
	my $r = $ratios->[$startIdx + 1];
	
	if($randm < $r){
		# Passed the test
		print {$fhs->[$startIdx + 1]} "$head\n$seq\n";
		return recursiveSelection($randm, $startIdx + 1, $ratios, $fhs, $head, $seq);
	}
	return $randm;
}

sub GetFastaReadCount{
	my ($file) = @_;
	open(my $IN, "< $file");
	my $counter = 0;
	while(my $line = <$IN>){
		if($line =~ /^>/){
			$counter++;
		}
	}
	seek($IN, 0, 0);
	close $IN;
	return $counter;
}
