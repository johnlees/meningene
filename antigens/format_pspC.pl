#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;

# Usage ./format_pspC.pl gene_order.txt pspC_genes.fa output_prefix
# Outputs a database for use with srst2, and translated genes for use with
# blastp

my $gene_order_file = $ARGV[0];
my $gene_file = $ARGV[1];
my $output_prefix = $ARGV[2];

open(ORDER, $gene_order_file) || die("Could not open $gene_order_file: $!\n");

my @gene_names;
while (my $line_in = <ORDER>)
{
   chomp $line_in;
   push(@gene_names, $line_in);
}

close ORDER;

my $seqin = Bio::SeqIO->new( -format => 'fasta' , -file => $gene_file);
my $srst2_out = Bio::SeqIO->new( -format => 'fasta', -file => ">$output_prefix.dna.fa");
my $blastp_out = Bio::SeqIO->new( -format => 'fasta', -file => ">$output_prefix.aa.fa");

my $i = 1;
while (my $seqobj = $seqin->next_seq())
{
   my $srst2_obj = Bio::Seq->new( -display_id => "0__pspC__$gene_names[$i-1]__$i",
                                  -seq => $seqobj->seq());
   $srst2_out->write_seq($srst2_obj);

   my $blastp_obj = $seqobj->translate(-complete=>1);
   $blastp_obj->display_id($gene_names[$i-1]);
   $blastp_out->write_seq($blastp_obj);

   $i++;
}


exit(0);

