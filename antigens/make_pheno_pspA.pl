#!/usr/bin/perl -w

use strict;
use warnings;

my $pspC_file = $ARGV[0];
my $link = $ARGV[1];
my $bolt = $ARGV[2];

# Usage
# ./make_pheno.pl X_transposed.txt pneumo_link.txt phenotypes.txt

my @pspA_vals = (1, 2, 3, 4);

open(LINK, $link) || die("Could not open link file $link: $!\n");

my %links;
while (my $line_in = <LINK>)
{
   chomp $line_in;

   my ($hpgen, $lane) = split("\t", $line_in);
   $lane =~ m/^(\d+_\d+)#(\d+)$/;
   $lane = "$1_$2";

   $links{$lane} = $hpgen;
}

close LINK;

open(BOLT, $bolt) || die("Could not open bolt file $bolt: $!\n");
my $header = <BOLT>;

my %hpgen_names;
while (my $line_in = <BOLT>)
{
   chomp $line_in;

   my ($full_hpgen, @gubbins) = split(/\s+/, $line_in);
   if ($full_hpgen =~ m/^\d+_.\d+_(hpgen\d+)$/)
   {
      $hpgen_names{$1} = $full_hpgen;
   }
   else
   {
      last;
   }
}

close BOLT;

my @colnames;
foreach my $val (@pspA_vals)
{
   push(@colnames, "PSPA_$val");
}
print join(" ", "FID", "IID", @colnames) . "\n";

open(PSPC, $pspC_file) || die("Could not open pspC file $pspC_file: $!\n");
while (my $line_in = <PSPC>)
{
   chomp $line_in;

   my ($lane, $pspA_allele) = split(/\s+/, $line_in);

   my $sample_id = $hpgen_names{$links{$lane}};

   if (defined($sample_id))
   {
      my @factors;

      foreach my $val (@pspA_vals)
      {
         if ($pspA_allele == $val)
         {
            push(@factors, 1);
         }
         else
         {
            push(@factors, 0);
         }
      }

      print join(" ", $sample_id, $sample_id, @factors) . "\n";
   }
}

exit(0);

