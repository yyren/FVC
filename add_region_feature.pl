#!/usr/bin/perl
use strict;
use warnings;

my $snp = $ARGV[0];
my $indel = $ARGV[1];
my $out_snp=$ARGV[2];
my $out_indel=$ARGV[3];


my $seq_len = 10;#1000000000
my $hash_seed = 7;

my (%all_rec,%snp_rec,%indel_rec) = ();
print"Start store snp infor...\n";
open(SNP,"<$snp") or die "can not open $snp";
while(my $line = <SNP>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	$snp_rec{$recs[1]}{$recs[2]}='snp';
	$all_rec{$recs[1]}{$recs[2]}=0;
}
close SNP;


print"Start store indel infor...\n";
my %group_rec=();
open(INDEL,"<$indel") or die "can not open $indel";
while(my $line = <INDEL>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	$indel_rec{$recs[1]}{$recs[2]}='indel';
	$all_rec{$recs[1]}{$recs[2]}=0;
}
close INDEL;

print"Start merge snp and indel infor...\n";
foreach my $chr(keys %all_rec){
	foreach my $pos(keys %{$all_rec{$chr}}){
		my $group = &getIndex($seq_len,$hash_seed,$pos);
		my @pos_arr=();
		$pos_arr[0]=$pos;
		if(exists $group_rec{$chr}{$group}){
			push @{$group_rec{$chr}{$group}},@pos_arr;
		}else{
			@{$group_rec{$chr}{$group}}=@pos_arr;
		}
	}
}

print"Start out snp infor...\n";
open(SNP,"<$snp") or die "can not open $snp";
open(OSNP,">$out_snp") or die "can not open $out_snp";
while(my $line = <SNP>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	my $group=&getIndex($seq_len,$hash_seed,$recs[2]);
	
	my @positions=@{$group_rec{$recs[1]}{$group}};
	my $tag_label=0;
	if($#positions>=1){
		my $type_temp='snp';
		for(my $i=0;$i<=$#positions;$i++){
			if(($positions[$i]>=$recs[2]-10) and ($positions[$i]<=$recs[2]+10)){
				if(exists $indel_rec{$recs[1]}{$positions[$i]}){
					$tag_label=3;#snp+indel
					last;
				}elsif($positions[$i] !=$recs[2]){
					$tag_label=1;#snp+snp
				}
			}
		}
	}
	
	print OSNP "$line $tag_label\n";
}
close SNP;
close OSNP;


print"Start out indel infor...\n";
open(INDEL,"<$indel") or die "can not open $indel";
open(OINDEL,">$out_indel") or die "can not open $out_indel";
while(my $line = <INDEL>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	my $group=&getIndex($seq_len,$hash_seed,$recs[2]);
	
	my @positions=@{$group_rec{$recs[1]}{$group}};
	my $tag_label=0;
	if($#positions>=1){
		my $type_temp='indel';
		for(my $i=0;$i<=$#positions;$i++){
			if(($positions[$i]>=$recs[2]-10) and ($positions[$i]<=$recs[2]+10)){
				if($positions[$i] !=$recs[2]){
					if(exists $indel_rec{$recs[1]}{$positions[$i]}){
						$tag_label=2;#indel+indel
						last;
					}else{
						$tag_label=3;#indel+snp
					}
				}elsif(exists $snp_rec{$recs[1]}{$positions[$i]}){
					$tag_label=3;#indel+snp
				}
			}
		}
	}
	
	print OINDEL "$line $tag_label\n";
}
close INDEL;
close OINDEL;



sub getIndex()
{
	my $len = shift;
	my $seed = shift;
	my $num = shift;
	my $complete = substr('0'x($len - length($num)).$num,0,$seed);
	return int($complete);
}