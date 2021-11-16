#!/usr/bin/perl
# Description: separate the file according to the variant type
# Author: Yongyong Ren

my $infile=$ARGV[0];
my $out_snp=$ARGV[1];
my $out_indel=$ARGV[2];
my $type=$ARGV[3];

open(IN,"<$infile") or die "can not open $infile:$!";
open(OUT_SNP,">$out_snp") or die "can not open $out_snp:$!";
open(OUT_INDEL,">$out_indel") or die "can not open $out_indel:$!";
while(my $line=<IN>){
	chomp $line;
	if($line=~m/^#/){
		print OUT_SNP "$line\n";
		print OUT_INDEL "$line\n";
	}else{
		my @recs=split(' ',$line);
		$ref=$recs[3];
		$alt=$recs[4];
		if($type eq 'feature'){
			$ref=$recs[4];
			$alt=$recs[5];
		}
		if(($alt=~m/\*/) or ($alt=~m/\,/)){
			my @alts=split(',',$alt);
			
			if(((length($ref)!=length($alts[1])) and ($alts[1] ne '*')) or ((length($ref)!=length($alts[0])) and ($alts[0] ne '*'))){
				print OUT_INDEL "$line\n";
			}else{
				print OUT_SNP "$line\n";
			}
		}else{
			if(length($ref)==length($alt)){
				print OUT_SNP "$line\n";
			}else{
				print OUT_INDEL "$line\n";
			}
		}
	}

}
close IN;
close OUT_SNP;
close OUT_INDEL;