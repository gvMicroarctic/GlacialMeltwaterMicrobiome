#!/usr/bin/perl

use strict;

#curate viral contigs

#checkv output: viral_gene and host_gene
open(IN1, "./viral_contig/checkv/contamination.tsv") or die $!;

my %index;
my $line = 0;
while(defined(my $input = <IN1>)) {
	if ($line > 0) {
		my @info = split(/\t/, $input);
		my $contig = $info[0];
		my $length = $info[1];
		my $viral_genes = $info[3];
		my $host_genes = $info[4];
		$index{$contig}{'LENGTH'} = $length;
		$index{$contig}{'VIRAL_GENES'} = $viral_genes;
		$index{$contig}{'HOST_GENES'} = $host_genes;
	}
	$line++;
}
close(IN1);

#virsorter2 output: score and hallmark
open(IN2, "./viral_contig/virsorter2_2/final-viral-score.tsv") or die $!;

$line = 0;
while(defined(my $input = <IN2>)) {
	if ($line > 0) {
		my @info = split(/\t/, $input);
		my $contig = $info[0];
		my $score = $info[3];
		my $hallmark = $info[6];
		$index{$contig}{'SCORE'} = $score;
		$index{$contig}{'HALLMARK'} = $hallmark;
	}
	$line++;
}
close(IN2);

#screen contigs
my $keep = 0;
my $manual = 0;
my $discard = 0;

open(my $fileK , ">./viral_contig/contig_keep.txt") or die $!;
open(my $fileM , ">./viral_contig/contig_manual.txt") or die $!;
open(my $fileD , ">./viral_contig/contig_discard.txt") or die $!;

my %manual_contig;

foreach my $contig (keys %index) {
	
	if (($index{$contig}{'VIRAL_GENES'} > 0) or (($index{$contig}{'VIRAL_GENES'} == 0) and ($index{$contig}{'HOST_GENES'} == 0)) or (($index{$contig}{'VIRAL_GENES'} == 0) and ($index{$contig}{'SCORE'} >= 0.95)) or (($index{$contig}{'VIRAL_GENES'} == 0) and ($index{$contig}{'HALLMARK'} > 2))) { #keep
		$keep++;	
		print $fileK "$contig\n";
	} elsif (($index{$contig}{'VIRAL_GENES'} == 0) and ($index{$contig}{'HOST_GENES'} == 1) and ($index{$contig}{'LENGTH'} >= 10000)) { #manual check
		$manual++;
		$manual_contig{$contig} = '';
		print $fileM "$contig\n";
	} else { #discard
		$discard++;
		print $fileD "$contig\n";
	}
	
}
close($fileK);
close($fileM);
close($fileD);

print "Keep: $keep\n";
print "Manual curation: $manual\n";
print "Discard: $discard\n";
