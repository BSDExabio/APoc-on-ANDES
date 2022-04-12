#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $help;
my $sort;
my $column = 2;

##################################################################
GetOptions(
             's|sort'       => \$sort,
             'c|column=i'   => \$column,
             'h|help'       => \$help
          );

##################################################################
printHelp() if( $help );

if ( $#ARGV != 1 ) {
  printHelp();
}

my $setfile = $ARGV[0];
my $numJob  = $ARGV[1];

printHelp() unless ( -s $setfile );
my @set = read_inp_file( $setfile );

my $numentry = scalar @set;

if ( $numJob >= 1000 ) { die "Insane job numbers! Please check it...\n"; }
if ( $numJob > $numentry ) { die "Job number is larger than the number of entries\n"; }

if( $sort ) {
    my @sorted_set = sort { $b->{'len'}<=>$a->{'len'} } @set;
    @set = @sorted_set;
}

my @jobfiles = ();
for(my $i=1; $i<=$numJob; $i++) {
  my $filename = "$setfile.$i";
  open my $file, '>', "$filename" or die "Error: could not open $filename\n";
  push (@jobfiles, $file);
}


my $num = scalar @set;
print "Found $num entries in the $setfile, spliting them into $numJob subset\n";
my $round = 0;
for(my $counter=0; $counter < $numentry; $counter++) {
  my $jobIndex = $counter % $numJob;
  $round++ if( $jobIndex == 0 );
  if( $round % 2  == 0 ) {
    $jobIndex = $numJob - $jobIndex - 1;
  }
  my $file = $jobfiles[$jobIndex];
  #print "$counter $set[$counter]->{ln} $jobIndex\n";
  print $file "$set[$counter]->{'ln'}\n";
}

for (@jobfiles) { close $_; }

##################################################################

#####################################################
sub read_inp_file {
  my $file = shift @_;

  my @set = ();
  open my $fh, "<$file" or die "Error: could not open $file\n";
  while(<$fh>){
    chomp;
    next if /name/;
    next if /^#/;
    next unless /\S/;

    my @fields = split(' ',$_);
    my $name   = $fields[0];
    my $len    = $fields[$column-1] if( $sort );

    if( $sort ) {
      push( @set, {'len'=>$len, 'ln'=>$_} );
    }
    else {
      push( @set, {'ln'=>$_} );
    }

  }
  close $fh;

  return @set;
}


sub printHelp {
  print "\nSplit an input set into subsets\n\n";
  print "$0 <options> <file> <number of subsets>\n\n";
        print "Options: -s Sort set, default using numbers supplied in column 2\n";
        print "Options: -c <column>, change column for sorting\n";
        print "         -h Print this help message\n";
  exit;
}
