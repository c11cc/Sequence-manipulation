#this is for calculating the GC content of target length sequence around the target genome site. 
use strict;
use warnings;

print "enter the name of the file containing genome and site information:\t";
my $infile=<>;

print "enter genome file:\t";
my $file=<>;

print "half region length:\t";
my $len=<>;

print "GC box length:\t";
my $box=<>;
chomp $box;

my $remain;
if($box%2)
  {
    $box=int($box/2);	
    $remain=1;
  }
else
  {
    $box=$box/2;
  }

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
      	my $na;
      	my $tmp= length $gme{$chr},;
      	if($site-$len > 0 and $site+$len+1 <= $tmp)#sites from the first complete box to the last
      	  {
    	      $info{$chr}->{$site}->{0}=substr $gme{$chr},$site-$len,2*$len+1;#get sequence
    	      get_gc(length $info{$chr}->{$site}->{0},$chr,$site);
    	    }
    	  elsif($site-$len <= 0)#start region of the genome
    	    {
    	      $info{$chr}->{$site}->{0}=substr $gme{$chr},1,$site+$len;    	      
    	      get_gc(length $info{$chr}->{$site}->{0},$chr,$site);
    	      foreach my $i(0..$len-$site)
    	        {
    	          $na.="na\t";
    	        }
    	      $info{$chr}->{$site}->{1}=$na.$info{$chr}->{$site}->{1};
    	    }
    	  elsif($site-$len > 0 and $site+$len+1 > $tmp)#end region of the genome
    	    {
    	      $info{$chr}->{$site}->{0}=substr $gme{$chr},$site-$len,$tmp-$site+$len;
    	      get_gc(length $info{$chr}->{$site}->{0},$chr,$site);
    	      foreach my $i(0..$len+$site-$tmp)
    	        {
    	          $na.="na\t";
    	        }
    	      $info{$chr}->{$site}->{1}.=$na;    	          	 
    	    }
    	  else
    	    {
    	      print "error\t$chr\t$site\n";	
    	    }  
      }	
  }

#calculate gc content
sub get_gc 
  {
  	my $slen=shift;#length
  	my $schr=shift;#genome
  	my $ssite=shift;#site
  	my $gc=0;
  	my $sj=-1;#relative location
  	my $n=0;
    my $tmpv;
  	while(++$sj<=$box)#GC content of next box length of each base,start region
      {
        $gc=0;
        $n=0;
        $tmpv=$sj+$box if($remain);
        $tmpv=$sj+$box-1 if(!$remain);
      	foreach my $k(0..$tmpv) #initialize,start
      	  {
      	    $gc++ if (substr $info{$schr}->{$ssite}->{0},$k,1)=~/[GC]/; #the base is g/c
      	    $n++ if(substr $info{$schr}->{$ssite}->{0},$k,1)=~/N/;	#the base is n
      	  }	    
      	$info{$schr}->{$ssite}->{1}.=($gc/($tmpv+1-$n))."\t" and next if($gc and $tmpv+1-$n);
      	$info{$schr}->{$ssite}->{1}.="0\t" and next if(!$gc and $tmpv+1-$n);  
      	$info{$schr}->{$ssite}->{1}.="na\t" and next if(!($tmpv+1-$n));
      }
    $sj--;
    while(++$sj<=$slen-1-$box)#GC content of next box length of each base,k to k+1
      {
      	$gc=0;
      	$n=0;
        $tmpv=$sj+$box if($remain);
        $tmpv=$sj+$box-1 if(!$remain);      	
      	foreach my $k($sj-$box..$tmpv)#
      	  {
      	    $gc++ if (substr $info{$schr}->{$ssite}->{0},$k,1)=~/[GC]/;
      	    $n++ if(substr $info{$schr}->{$ssite}->{0},$k,1)=~/N/; 	
      	  }
      	$info{$schr}->{$ssite}->{1}.=($gc/(2*$box-$n))."\t" and next if($gc and 2*$box != $n and !$remain);
      	$info{$schr}->{$ssite}->{1}.=($gc/(2*$box+1-$n))."\t" and next if($gc and 2*$box+1 != $n and $remain);
      	$info{$schr}->{$ssite}->{1}.="0\t" and next if(!$gc and 2*$box !=$n and !$remain);
      	$info{$schr}->{$ssite}->{1}.="0\t" and next  if(!$gc and 2*$box+1!=$n and $remain);
      	$info{$schr}->{$ssite}->{1}.="na\t" and next if( 2*$box ==$n and !$remain);
      	$info{$schr}->{$ssite}->{1}.="na\t" and next if( 2*$box+1==$n and $remain);
      } 
    $sj--;   
    foreach $sj($slen-$box..$slen-1)
      {
      	$gc=0;
      	$n=0;
        $tmpv=$slen if($remain);
        $tmpv=$slen-1 if(!$remain);      	
        foreach my $k($sj-$box..$tmpv)
      	  {
      	    $gc++ if(substr $info{$schr}->{$ssite}->{0},$k,1)=~/[GC]/; 
      	    $n++ if(substr $info{$schr}->{$ssite}->{0},$k,1)=~/N/;	
      	  }	
      	$info{$schr}->{$ssite}->{1}.=($gc/($slen-$sj+$box-$n))."\t" and next if($gc and $slen-$sj+$box-$n);	
      	$info{$schr}->{$ssite}->{1}.="0\t" and next if(!$gc and $slen-$sj+$box-$n);
      	$info{$schr}->{$ssite}->{1}.="na\t" and next if(!($slen-$sj+$box-$n));
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