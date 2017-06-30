#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use Array::Utils qw(:all);
use Data::Dumper;
use POSIX;

my @sets=("WTCHG_350045_236_F2_s","WTCHG_350045_240_F2_s");

my $binsize=20;
open(all,"<", "AsiSI_gammaplus94_genic_strand.bed")or die("unable to open all");

my %coords;
my %coords_HR;
my %coords_NHEJ;



while(my $line=<all>)
{
chomp($line);
my @array=split(/\t/, $line);
print $array[5]."\n";
	if($array[5]eq"+")
	{
	push @{$coords{$array[0]}}, [$array[1]-2000, $array[1]+2000, $array[5],$array[0], $array[3]];
		print $array[1]." ".$array[5]." ".$array[0]." ".$array[3]."\n";
	}
	elsif($array[5]eq"-")
	{
	push @{$coords{$array[0]}}, [$array[2]-2000, $array[2]+2000, $array[5],$array[0], $array[3]];
		print $array[2]." ".$array[5]." ".$array[0]." ".$array[3]."\n";

	}
	

}
close(all);

for(my $rr=0; $rr<=$#sets; $rr++)
{
my $dataset=$sets[$rr];
print $sets[$rr]."\n";


my %wig;

my $q=0;
my $NF = qx(wc -l $dataset.bedgraph);
open(F,"<",$dataset.".bedgraph")or die("unable to open F");
<F>;
while(my $line=<F>)
{
$q+=1/$NF;
print $q." F\n";

chomp($line);
my @array=split(/\t/, $line);

for(my $u=$array[1]; $u<$array[2]; $u++)
{
$wig{$array[0]}[$u]=$array[3];
}

}
close(F);


system("mkdir antisense");
open(sense_F, ">", "metagenes_final/sense_F_".$sets[$rr]."_AsiSI_hg38_intersect_94")or die("unable to open sense_F");
open(anti_F, ">", "metagenes_final/anti_F_".$sets[$rr]."_AsiSI_hg38_intersect_94")or die("unable to open anti_F");

my $key;
foreach $key (keys  %coords)
{

{
		if (exists $wig{$key})
		{
		
print $key." ".$#{$coords{$key}}."\n";
	for (my $u=0; $u<=$#{$coords{$key}}; $u++)
	{
		if ($coords{$key}[$u][2]eq"+")
		{

					for(my $i=$coords{$key}[$u][0]; $i<=$coords{$key}[$u][1]; $i++)
					{
					if ($i<=$#{$wig{$key}} && $i>=0)
					{
					my $a=$wig{$key}[$i]+0;
               		print sense_F $a."\t";
					}
					else
					{
					print sense_F "0\t";
					}
					}
					print sense_F "\n";
					
					
					for(my $i=$coords{$key}[$u][0]; $i<=$coords{$key}[$u][1]; $i+=$binsize)
					{
						my $sum_bin=0;
						for(my $j=0; $j<$binsize; $j++)
						{
							if ($i+$j<=$#{$wig{$key}} && $i+$j>=0)
							{
							$sum_bin+=$wig{$key}[$i+$j]/$binsize;	
							}	
						}
						for(my $j=0; $j<$binsize; $j++)
						{
							if ($i+$j<=$#{$wig{$key}} && $i+$j>=0)
							{
							print sense_F_bin $sum_bin."\t";
							}
							else
							{
							print sense_F_bin "0\t";
							}
						}
					}
					print sense_F_bin "\n";
					
		}
				
		elsif ($coords{$key}[$u][2]eq"-")
		{

					for(my $i=$coords{$key}[$u][1]; $i>=$coords{$key}[$u][0]; $i--)
					{
					if ($i<=$#{$wig{$key}} && $i>=0)
					{
					my $a=-$wig{$key}[$i]+0;
               		print anti_F $a."\t";
					}
					else
					{
					print anti_F "0\t";
					}
					}
					print anti_F "\n";
					
					for(my $i=$coords{$key}[$u][1]; $i>=$coords{$key}[$u][0]; $i-=$binsize)
					{
						my $sum_bin=0;
						for(my $j=0; $j<$binsize; $j++)
						{
							if ($i-$j<=$#{$wig{$key}} && $i-$j>=0)
							{
							$sum_bin+=$wig{$key}[$i-$j]/$binsize;	
							}	
						}
						for(my $j=0; $j<$binsize; $j++)
						{
							if ($i-$j<=$#{$wig{$key}} && $i-$j>=0)
							{
							print anti_F_bin $sum_bin."\t";
							}
							else
							{
							print anti_F_bin "0\t";
							}
						}
					}
					print anti_F_bin "\n";
					
		}
				
		print $u."\n";
	}
		}
}
}
close(sense_F);
close(anti_F);
close(sense_F_bin);
close(anti_F_bin);

undef(%wig)
}