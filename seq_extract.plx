#this is for extracting the sequence around the target site. there will be error if the chrmosome idenitfier from the genome file and target file are different.
use strict;
use warnings;

print "enter the file name contain chromosome identifier and site information:\t";
my $infile=<>;

print "enter genome file(.fna):\t";
my $file=<>;

print "half region length:\t";
my $len=<>;

#read in target file
my %info;
open IN, "$infile" or die $!;
while(<IN>)
  {
  	chomp;
  	if($_=~/([^\t]+)\t([^\t]+)/)
  	  {
  	  	 $info{$1}->{$2}=1;
  	  }
  	else
  	  {
  	    print "error $_\n";
  	  }
  }
close IN;

#read in genom file
my $chr;
my %gme;
open (my $in, "$file") or die $!;
while(<$in>)
  {
    chomp;
    if($_=~/>([^\s]+)\s/)
      {
        $gme{$1}->{seq}=0;
        $chr=$1;
      }  
    elsif($_=~/>([^\t]+)\t/)
      {
        $gme{$1}->{seq}=0;
        $chr=$1;
      }    
    if($_!~/>/)
      {
        $gme{$chr}->{seq}.=$_;
      }
  }
close $in;

#transform to uppercase
foreach $chr(keys %gme)
  { 
    $gme{$chr}->{seq}=~tr/[a-z]/[A-Z]/;
  }

#get the sequence and write
open OUT, ">target_seq.fna" or die $!;
foreach $chr(sort keys %info)
  {
    foreach my $site (sort keys %{$info{$chr}})
      {
      	print OUT ">$chr\t$site\n";
      	my $tmp=length $gme{$chr}->{seq};
      	if($site-$len > 0 and $site+$len+1 <= $tmp)#normal
      	  {
    	      $info{$chr}->{$site}=substr $gme{$chr}->{seq},$site-$len,2*$len+1;
    	    }
    	  elsif($site-$len <= 0)#start
    	    {
    	      $info{$chr}->{$site}=substr $gme{$chr}->{seq},1,$site+$len;
    	      $info{$chr}->{$site}=~tr/0//;
    	    }
    	  elsif($site-$len > 0 and $site+$len+1 > $tmp)#end
    	    {
    	      $info{$chr}->{$site}=substr $gme{$chr}->{seq},$site-$len,$tmp-$site+$len;
    	    }
    	  else
    	    {
    	      print "error\t$chr\t$site\n";	
    	    }
    	  print OUT "$info{$chr}->{$site}\n";
      }	
  }
close OUT;

