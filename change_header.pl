#!/usr/bin/perl
use strict;
use warnings;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];
my $sid=$ARGV[2];

my %header_rec=();
open(IN,"<$infile") or die "can not open $infile";
open(OUT,">$outfile") or die "can not open $outfile";
while(my $line=<IN>){
	chomp $line;
	if($line=~m/^\#/){
		if($line=~m/^##/){
			my @headers=split(/\,/,$line);
			if(not exists $header_rec{$headers[0]}){
				print OUT "$line\n";
			}
			$header_rec{$headers[0]}=0;
		}elsif($line=~m/^#CHROM/){
			my @recs=split(/\t/,$line);
			my $new_line=join("\t",(@recs[0..8],$sid));
			print OUT "$new_line\n";
		}
	}else{
		print OUT "$line\n";
	}
	
	
}