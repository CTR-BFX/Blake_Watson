#!/usr/bin/perl

use strict;

my $Interval  = $ARGV[0]; # e.g. 5000 for 5Kb
my $GenomeIDX = $ARGV[1]; # e.g. Mus_musculus.GRCm38.dna.chromosome.all.fa.fai";



#print "Interval Size = ", $ARGV[0], "\n";

my %Chromosomes;
open(IDX, $GenomeIDX) || die "Can't open $GenomeIDX for reading: $!\n";
while(<IDX>)
  {
    my ($chr,$len) = split(/\s+/, $_);

    $Chromosomes{$chr} = $len; 
  }
close IDX;


my $i;
foreach my $chr (sort keys %Chromosomes)
  {
     for($i=0; $i<=$Chromosomes{$chr}; $i+=$Interval)
        {
          my $end = $i+$Interval;
          if($i+$Interval > $Chromosomes{$chr}){ $end = $Chromosomes{$chr}; }
          my $start = $i+1;
          if($i == 0){ $start = $i; }

          if($chr == "MT"){ $chr = "M"; }       
  
          print "chr$chr\t$start\t$end\n";
        }
  }
