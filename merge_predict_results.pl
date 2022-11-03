#!/usr/bin/perl
use strict;
use warnings;

my $snp_predict_results=$ARGV[0];
my $indel_predict_results=$ARGV[1];
my $in_vcf_file=$ARGV[2];
my $out_file=$ARGV[3];

my $numb=0;
my (%fnvc_predict_snp,%fnvc_predict_indel)=();

open(IN,"<$snp_predict_results") or die "can not open $snp_predict_results\n";
while(my $line=<IN>){
	chomp $line;
	my @recs=split(/\t/,$line);
	if($numb!=0){
		my $id=join("\_",@recs[0,1,3,4]);
		$fnvc_predict_snp{$id}=$recs[6];
	}
	$numb++;
}
close IN;

$numb=0;
open(IN,"<$indel_predict_results") or die "can not open $indel_predict_results\n";
while(my $line=<IN>){
	chomp $line;
	my @recs=split(/\t/,$line);
	if($numb!=0){
		my $id=join("\_",@recs[0,1,3,4]);
		$fnvc_predict_indel{$id}=$recs[6];
	}
	$numb++;
}
close IN;

my $infor_endline='F';
my $fnvc_filter='##FILTER=<ID=FVC_Filtered,Description="The probability of the variant to be true is less than 0.5">';
my $fnvc_infor='##INFO=<ID=FVC,Number=1,Type=Float,Description="FNVC prediction results"';
open(IN, "<$in_vcf_file") or die "can not open $in_vcf_file";
open(OUT,">$out_file") or die "can not open $out_file";
while(my $line=<IN>){
	chomp $line;
	if($line=~m/^#/){
		if($line=~m/^##contig/ and $infor_endline eq 'F'){
			$infor_endline='T';
			print OUT "$fnvc_filter\n$fnvc_infor\n$line\n";
		}else{
			print OUT "$line\n";
		}
	}else{
		my @recs=split(/\t/,$line);
		my $id=join("\_",@recs[0..1,3..4]);
		my $predict_value=1.0;
		if(exists $fnvc_predict_indel{$id}){
			$predict_value=$fnvc_predict_indel{$id};
		}
		if(exists $fnvc_predict_snp{$id}){
			$predict_value=$fnvc_predict_snp{$id};
		}
		$recs[7]=$recs[7].';FVC='.$predict_value;
		if($predict_value>=0.5){
			$recs[6]='PASS'
		}else{
			$recs[6]='FVC_Filtered'
		}
		my $new_line=join("\t",@recs);
		print OUT "$new_line\n";
	}
}
close IN;
close OUT;
