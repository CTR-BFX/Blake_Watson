#!/usr/bin/perl

use strict;


my @TotalCoverage;

my $WindowSize = $ARGV[0];

my @files = <*chr.srtd.${WindowSize}.bed>;

print "echo \"Region";
foreach my $file (@files) 
  {

     my $totcov = `cat $file | awk '{ sum+=\$4} END {print sum}' | tr -d '\n'`;

     push(@TotalCoverage, $totcov);

     my $file_name = $file;
     $file_name    =~ s/.chr.srtd.*.bed//;
     print "\t" . $file_name;
  }
print "\" > CTR_edw23_0003_MeDIPS_Region_Coverage_Combo_${WindowSize}.table\n";



my $cnt = 0;
my $command;
foreach my $file (@files) 
  {
    $command .= "paste <(awk '{print \$1\"_\"\$2\"_\"\$3\"\t\"1000000*(\$4/$TotalCoverage[$cnt]) }' $file) ", if($cnt == 0);

    $command .= "<(awk '{print 1000000*(\$4/$TotalCoverage[$cnt])}' $file) ", if($cnt > 0);

    $cnt++;
  }


$command .= "| grep -v '^chrM' | column >> CTR_edw23_0003_MeDIPS_Region_Coverage_Combo_${WindowSize}.table";

print "\n\n";
print $command;
print "\n";

