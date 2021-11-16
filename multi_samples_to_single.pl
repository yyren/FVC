#!/usr/bin/perl

## function: split the VCF file with multiple samples into multiple VCF files with single sample according the sampleID in header line
## out: ${sid}.split.vcf
## usage e.g: 
## perl multi_samples_to_single.pl input.vcf /data/out_vcf


use strict;
use warnings;

my $infile=$ARGV[0];
my $out_dir=$ARGV[1];

$out_dir=~s/\/$//;
#print "$out_dir\n";
my $vcf_header=`awk '/^#CHROM/ {print}' $infile`;

$vcf_header=~s/\r|\n//g;
my @header_recs=split(/\t/, $vcf_header);
my @sampleIds=@header_recs[9..$#header_recs];

## checking and remove existing files ##
for(my $i=0;$i<=$#sampleIds;$i++){
	#print "3: $sampleIds[$i]";
	my $out_file=$out_dir.'/'.$sampleIds[$i].'.split.vcf';
	if(-e $out_file){
		`rm $out_file`;
		
	}
}


my %out;
for(my $i=0;$i<=$#sampleIds;$i++){
	my $out='OUT'.$i;
	my $out_file=$out_dir.'/'.$sampleIds[$i].'.split.vcf';
	open($out{$i}, ">>$out_file") or die "can not open $out";
}

open(IN,"<$infile") or die "can not open $infile";
while(my $line=<IN>){
	chomp $line;
	if($line=~m/^#/){
		#print "$#sampleIds\n";
		for(my $i=0;$i<=$#sampleIds;$i++){
			my $out_file=$out_dir.'/'.$sampleIds[$i].'.split.vcf';
			if($line=~m/^##/){
				if(($line!~m/^##INFO/) and ($line!~m/^##FORMAT/)and ($line!~m/^##FILTER/)){
					$out{$i}->print("$line\n");
				}
			}else{
				my @recs=split(/\t/,$line);
				my $idx=$i+9;
				my $header=join("\t",@recs[0..8,$idx]);
				$out{$i}->print("$header\n");
			}
		}
	}else{
		my @recs=split(/\t/,$line);
		for(my $i=0;$i<=$#sampleIds;$i++){
			my $out_file=$out_dir.'/'.$sampleIds[$i].'.split.vcf';
			my $idx=$i+9;
			my $new_line=join("\t",(@recs[0..5],'.','.','.',$recs[$idx]));
			$out{$i}->print("$new_line\n");
		}
	}
}
close IN;

for(my $i=0;$i<=$#sampleIds;$i++){
	close $out{$i};
}










