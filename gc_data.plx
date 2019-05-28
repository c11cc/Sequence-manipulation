#this is for calculating the GC content of target length sequence around the target genome site. for GC box that has an odd value, it calculates the exact region. for GC box that has an even value, it calculate the region of length +1(box length/2 after and before the site and the site) 
use strict;
use warnings;

print "enter the name of the file containg genome and site information:\t";
my $infile=<>;

print "enter genome file:\t";
my $file=<>;

print "half region length:\t";
my $len=<>;

print "GC box length:\t";
my $box=<>;

my %info;
open IN, "$infile" or die $!;
while(<IN>)
  {
  	chomp;
  	if($_=~/([^\t]+)\t([^\t]+)/)
  	  {
  	  	 $info{$1}->{$2}->{0}=1;
  	  }
  	else
  	  {
  	    print "error $_\n";
  	  }
  }
close IN;

my $chr;
my %gme;
open (my $in, "$file") or die $!;
while(<$in>)
  {
    chomp;
    if($_=~/>([^\s]+)\s/)
      {
        $gme{$1}=0;
        $chr=$1;
      }  
    elsif($_=~/>([^\t]+)\t/)
      {
        $gme{$1}=0;
        $chr=$1;
      }    
    if($_!~/>/)
      {
        $gme{$chr}.=$_;
      }
  }
close $in;

#transform to uppercase
foreach $chr(keys %gme)
  { 
    $gme{$chr}=~tr/[a-z]/[A-Z]/;
  }

foreach $chr(sort keys %info)
  {
    foreach my $site (sort keys %{$info{$chr}})
      {
      	my $tmp= length $gme{$chr},;
      	if($site-$len > 0 and $site+$len+1 <= $tmp)#sites from the first complete box to the last
      	  {
    	      $info{$chr}->{$site}->{0}=substr $gme{$chr},$site-$len,2*$len+1;#get sequence
    	    }
    	  elsif($site-$len <= 0)#start region of the genome
    	    {
    	      $info{$chr}->{$site}->{0}=substr $gme{$chr},1,$site+$len;    	      
    	    }
    	  elsif($site-$len > 0 and $site+$len+1 > $tmp)#end region of the genome
    	    {
    	      $info{$chr}->{$site}->{0}=substr $gme{$chr},$site-$len,$tmp-$site+$len;
    	    }
    	  else
    	    {
    	      print "error\t$chr\t$site\n";	
    	    }
    	  get_gc(length $info{$chr}->{$site}->{0},$chr,$site);  
      }	
  }

sub get_gc
  {
  	my $slen=shift;#length
  	my $schr=shift;#genome
  	my $ssite=shift;#site
  	my $gc=0;
  	my $sj=0;
  	while($sj<=int($box/2))#GC content of next box length of each base,start region
      {
        $gc=0;
      	foreach my $k(1..$sj+int($box/2))#initialize,start
      	  {
      	    $gc++ if (substr $info{$schr}->{$ssite}->{0},$k,1)=~/[GC]/; 	
      	  }	      	  
      	$info{$schr}->{$ssite}->{1}.=($gc/($sj+int($box/2)))."\t" if($gc);
      	$info{$schr}->{$ssite}->{1}.="0\t" if(!$gc);  
      	$sj++;
      }
    while($sj<=$slen-1-int($box/2))#GC content of next box length of each base,k to k+1
      {
      	$gc=0;
      	foreach my $k($sj-int($box/2)..$sj+int($box/2))#
      	  {
      	    $gc++ if (substr $info{$schr}->{$ssite}->{0},$k,1)=~/[GC]/; 	
      	  }
      	$info{$schr}->{$ssite}->{1}.=($gc/$box)."\t" if($gc);
      	$info{$schr}->{$ssite}->{1}.="0\t" if(!$gc);
      	$sj++;
      }    
    foreach $sj($slen-int($box/2)..$slen-1)
      {
      	$gc=0;
        foreach my $k($sj..$slen-1)
      	  {
      	    $gc++ if (substr $info{$schr}->{$ssite}->{0},$k,1)=~/[GC]/; 	
      	  }	
      	$info{$schr}->{$ssite}->{1}.=($gc/($slen-1-$sj+int($box/2)))."\t" if($gc);	
      	$info{$schr}->{$ssite}->{1}.="0\t" if(!$gc);
      }
  }


#write out result
open OUT, ">gc_$infile" or die $!;
print OUT "chr\tsite\t";
foreach my $i (0..2*$len)
  {
  	print OUT "loc",$i-$len,"\t";
  }
print OUT "\n";
foreach $chr(sort keys %info)
  {
    foreach my $site (sort keys %{$info{$chr}})
      {
      	print OUT "$chr\t$site\t$info{$chr}->{$site}->{1}\n";
      }    
  }
close OUT;