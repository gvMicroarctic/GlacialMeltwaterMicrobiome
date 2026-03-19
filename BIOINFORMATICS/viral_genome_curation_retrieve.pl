#!/usr/bin/perl

use strict;

#curated viral contigs - retrieve contigs

open(IN1, "./viral_contig/contig_keep.txt") or die $!;
open(IN2, "./viral_contigs_unique.fasta") or die $!;
open(my $file, ">./viral_contig/viral_contigs_unique_final.fasta") or die $!;

my %index;

#get keep contigs
while(defined(my $input = <IN1>)) {
	chomp($input);
	$index{$input} = "";
}
close(IN1);

#get contig sequences
my $inside = 0;
my $count = 0;
while(defined(my $input = <IN2>)) {
	chomp($input);
	if ($input =~ /^>/) {
		my ($contig) = $input =~ />(.*)/;
		if (defined($index{$contig})) {
			$inside = 1;
			$count++;
			delete($index{$contig});
		} else {
			$inside = 0;
		}
	}
	if ($inside == 1) {
		print $file "$input\n";
	}
}
close(IN2);
close($file);

print "Number of retrieved contigs: $count\n";
