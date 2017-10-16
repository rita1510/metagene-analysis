#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use Array::Utils qw(:all);
use Data::Dumper;
use POSIX;
my @sets=("WTCHG_417720_289","WTCHG_417720_290","WTCHG_417720_293","WTCHG_417720_294","WTCHG_350045_236","WTCHG_350045_240");
my $binsize=20;




open(all,"<", "AsiSI_gammaminus94_genic_strand.bed")or die("unable to open all");

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
my $NR = qx(wc -l $dataset/R2_s.bedgraph);
open(R,"<",$dataset."/R2_s.bedgraph")or die("unable to open R");
<R>;
while(my $line=<R>)
{
$q+=1/$NR;
print $q." R\n";

chomp($line);
my @array=split(/\t/, $line);

for(my $u=$array[1]; $u<$array[2]; $u++)
{
$wig{$array[0]}[$u]=-$array[3];
}

}
close(R);


open(sense_R, ">", "metagenes_ELIFE/sense_R_".$sets[$rr]."_AsiSI_hg38_gammaminus94")or die("unable to open sense_R");
open(anti_R, ">", "metagenes_ELIFE/anti_R_".$sets[$rr]."_AsiSI_hg38_gammaminus94")or die("unable to open anti_R");

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
					my $a=-$wig{$key}[$i]+0;
               		print anti_R $a."\t";
					}
					else
					{
					print anti_R "0\t";
					}
					}
					print anti_R "\n";
					
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
							print anti_R_bin $sum_bin."\t";
							}
							else
							{
							print anti_R_bin "0\t";
							}
						}
					}
					print anti_R_bin "\n";
		}
			
		elsif ($coords{$key}[$u][2]eq"-")
		{

					for(my $i=$coords{$key}[$u][1]; $i>=$coords{$key}[$u][0]; $i--)
					{
					if ($i<=$#{$wig{$key}} && $i>=0)
					{
					my $a=$wig{$key}[$i]+0;
					print sense_R $a."\t";
					}
					else
					{
					print sense_R "0\t";
					}
					}
					print sense_R "\n";
					
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
							print sense_R_bin $sum_bin."\t";
							}
							else
							{
							print sense_R_bin "0\t";
							}
						}
					}
					print sense_R_bin "\n";
				
				
				
		}
		print $u."\n";
	}
		}
}
}
close(anti_R);
close(sense_R);
close(anti_R_bin);
close(sense_R_bin);



undef(%wig)
}