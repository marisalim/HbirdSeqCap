use warnings;
use strict;
use File::Basename;


die (qq/

Usage: average_thetas.pl <folder with all per site thetas>
 
\n/) if !@ARGV;


my $orig_folder;
  
  if ($ARGV[0] =~ m/\/$/ ){
    $orig_folder = $ARGV[0]; 
  }
  else {
    $orig_folder = $ARGV[0] . "/";
  }


my @file = <$orig_folder*>;
foreach (@file) {
my $file = $_;
my $lib = $1 if basename($file) =~ m/(\S+)\./;
my $prefix = basename($file);
open (IN, "<", $file);
my $pi;
my $S;
my $count;
while (<IN>) {
  chomp (my @line = split /\s+/, $_);
 
    $pi += exp ($line[3]);
    $S += exp ($line[2]);
    
  $count++;
}
close IN;

print "$prefix:average pi ", $pi/$count, "\n";
print "$prefix:average thetaW " , $S/$count, "\n";

}
