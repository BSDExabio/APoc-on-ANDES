#!/usr/bin/perl -w

use strict;
use File::Basename;

my $input_file = $ARGV[0];

my $proj_dir = $ENV{"PROJWORK"};


my $dataroot = "/gpfs/alpine/world-shared/bif135/species/R_rubrum/apoc";

my $subdir = "apoc_pdb70_210721";

open SET, "<$input_file";


while(<SET>) {
  next if /^#/;
  chomp;
  my @fields = split(' ', $_);
  my $name = $fields[0];

  my $datapath = "$dataroot/$name/$subdir";
  my $datafile = "$datapath/$name\_aln.dat.gz";

  unless ( -s $datafile ) {
    print "$name Error: cannot find the results $datafile\n";
    next;
  }
  else{
    #print "$name done\n";
  }

  my $num = `less $datafile |grep "### Alignment" |wc -l`;
  chomp $num;
  if( $num != 100 ) {
     print "$name Warning: found only $num lines\n";
  }
}

close SET;

