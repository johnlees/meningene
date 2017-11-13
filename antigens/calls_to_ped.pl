#!/usr/bin/perl -w

use strict;
use warnings;

# usage
# ./calls_to_ped.pl pspC_predictions.txt ../mapping/snps.fam pspC 3000000

my $call_file = $ARGV[0];
my $fam_file = $ARGV[1];
my $output_prefix = $ARGV[2];
my $offset = $ARGV[3]; # To make sure positions do not overlap

# Read fam into hash by sample id
open(FAM, $fam_file) || die("Could not open fam $fam_file\n");
my %fam;
while (my $line_in = <FAM>)
{
   chomp $line_in;
   my @fam_fields = split(/\s+/, $line_in);

   @{$fam{$fam_fields[0]}} = @fam_fields;
}

# Read calls into hash by sample id
open(CALLS, $call_file) || die("Could not open calls $call_file\n");
my %poss_alleles;
my %calls;
while (my $line_in = <CALLS>)
{
   chomp $line_in;
   my @call_fields = split(/\s+/, $line_in);
   @{$calls{$call_fields[0]}} = @call_fields;

   # Work out number of alleles
   for (my $i = 1; $i < scalar(@call_fields); $i++)
   {
      if (!defined($poss_alleles{$i-1}{$call_fields[$i]}))
      {
         $poss_alleles{$i-1}{$call_fields[$i]} = 1;
      }
   }
}

open(PED, ">$output_prefix.ped") || die("Could not write to $output_prefix.ped\n");
foreach my $sample (sort keys %calls)
{
   my $ped_line = join(" ", @{$fam{$sample}});
   foreach my $type (sort keys %poss_alleles)
   {
      foreach my $allele (sort keys %{$poss_alleles{$type}})
      {
         if ($calls{$sample}[$type+1] eq $allele)
         {
            $ped_line .= " A A";
         }
         else
         {
            $ped_line .= " C C";
         }
      }
   }
   print PED $ped_line . "\n";
}

open(MAP, ">$output_prefix.map") || die("Could not write to $output_prefix.map\n");
my $i = $offset;
foreach my $type (sort keys %poss_alleles)
{
   foreach my $allele (sort keys %{$poss_alleles{$type}})
   {
      print MAP join(" ", "26", "$output_prefix.$type.$allele", "0", "$i") . "\n";
      $i++;
   }
}

exit(0);

