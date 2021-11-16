#!/usr/bin/perl
use strict;
use warnings;

my $vcf = $ARGV[0];
my $bed = $ARGV[1];
my $inbed = $ARGV[2];
my $out_of_bed = $ARGV[3];

my $seq_len = 10;#1000000000
my $hash_seed = 7;

my %chr = ();
my %header_rec=();

open BED,"<$bed";
while(my $line = <BED>)
{
	chomp $line;
	my @arr = split/\t/,$line;
	if($arr[0]=~/^chr/){
		$arr[0]=~s/chr//;
	}
	if($arr[0]=~/^MT/){
		$arr[0]=~s/MT/M/;
	}
	my $low = &getIndex($seq_len,$hash_seed,$arr[1]);
	my $high = &getIndex($seq_len,$hash_seed,$arr[2]);
	if($low == $high)
	{
		if(exists $chr{$arr[0]}{$low})
		{
			push @{$chr{$arr[0]}{$low}},@arr[1..2];
		}
		else
		{
			@{$chr{$arr[0]}{$low}} = @arr[1..2];
		}
	}
	else
	{
		for(my $i=$low;$i<=$high;$i++)
		{
			if(exists $chr{$arr[0]}{$i})
			{
				push @{$chr{$arr[0]}{$i}},@arr[1..2];
			}
			else
			{
				@{$chr{$arr[0]}{$i}} = @arr[1..2];
			}
		}
	}
}
close BED;

open(VCF,"<$vcf") or die "can not open $vcf:$!";
open(OUT,">$inbed") or die "can not open $inbed:$!";
open(OUT2,">$out_of_bed") or die "can not open $out_of_bed:$!";
my $i = 0;
my $mid = 0;
my $curr_chr = 1;
my $curr_block = -1;
while(my $line = <VCF>)
{
	chomp $line;
	if($line =~m/^#/){
		my @headers=split(/\,/,$line);
		if(not exists $header_rec{$headers[0]}){
			print OUT "$line\n";
			print OUT2 "$line\n";
		}
		$header_rec{$headers[0]}=0;
	}else{
		$i++;
		if($i % 100000 == 0)
		{
			print"$i\n";
		}
		
		my @arr = split/\t/,$line;
		if($arr[0] =~m/^chr/)
		{
            #print"$arr[0]\n";
			$arr[0] =~s/^chr//;
			#$line=~s/^chr//;
		}
		if($arr[0]=~/^MT/){
			$arr[0]=~s/MT/M/;
			#$line=~s/^MT/M/;
		}
		my $block = &getIndex($seq_len,$hash_seed,$arr[1]);
		if($arr[0] ne $curr_chr)
		{
			$curr_chr = $arr[0];
			$curr_block = $block;
			$mid = 0;
		}
		else
		{
			if($block != $curr_block)
			{
				$mid = 0;
				$curr_block = $block;
			}
		}
		
		my @position = ();
		if(not exists $chr{$arr[0]}{$block})
		{
			print OUT2 "$line\n";
		}
		else
		{
			@position = @{$chr{$arr[0]}{$block}};
			my $low = $mid;
			my $high = $#position;
			my $left_pos = 0;
			my $in_bed = 0;
			if($arr[1] < $position[$low] || $arr[1] > $position[$high])
			{
				#print"abnormal: $position[$low]	$position[$high]	$arr[1]\n";
				print OUT2 "$line\n";
				next;
			}
			while($low <= $high)
			{
				#print"$low	$high	$arr[1]\n";
				$mid = int(($low + $high)/2);
				if($arr[1] == $position[$mid] || $arr[1] == $position[$mid+1])
				{
					$in_bed = 1;
					$left_pos = $mid;
					last;
				}
				elsif($arr[1] > $position[$mid] && $arr[1] < $position[$mid+1])
				{
					$left_pos = $mid;
					last;
				}
				elsif($arr[1] < $position[$mid])
				{
					$high = $mid - 1;
				}
				elsif($arr[1] > $position[$mid])
				{
					$low = $mid + 1;
				}
				
			}
			if($left_pos % 2 == 0|| $in_bed == 1)
			{
				print OUT "$line\n";
				#print $line,"\n";
				#exit;
			}
			else
			{
				print OUT2 "$line\n";
				
			}
		}

	}
}
close VCF;
close OUT;
close OUT2;

sub getIndex()
{
	my $len = shift;
	my $seed = shift;
	my $num = shift;
	my $complete = substr('0'x($len - length($num)).$num,0,$seed);
	return int($complete);
}