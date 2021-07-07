#!/usr/bin/perl
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p priority
#SBATCH -q msn
#SBATCH --mem=30000
use strict;
use File::Basename;

chomp(@ARGV);

if(scalar(@ARGV) < 2){
        print "usage: $0 <input fna folder> <temp prefix> <output file and path>\n";
        exit()
}

my @files = `ls $ARGV[0]/*.fa`;
my @temps;
foreach my $f (@files){
        chomp($f);
        my $t = basename($f);
        $t =~ s/\.fna//;
        chomp($t);
        print "~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -s 100000 -o mash/$ARGV[1]\_$t\_k21 $f\n";
        system("~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -s 100000 -o mash/$ARGV[1]\_$t\_k21 $f");
        push(@temps, "mash/$ARGV[1]\_$t\_k21.msh");
}

print("~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash paste $ARGV[2] mash/$ARGV[1]\_*.msh\n");
system("~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash paste $ARGV[2] mash/$ARGV[1]\_*.msh");

foreach my $t (@temps){
        system("rm $t");
}
