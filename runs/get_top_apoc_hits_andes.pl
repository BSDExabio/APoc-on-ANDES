#! /usr/bin/env perl
## Retrieve the top ranked hits from an APoc scoring output file

use warnings;
use strict;
use File::Path;
use File::Basename;
use Getopt::Long;

my $proj_dir = $ENV{"PROJWORK"};
#my $proj_dir = $ENV{"HOME"};


my $apoc_aln_root = "/gpfs/alpine/world-shared/bif135/species/R_rubrum/apoc";

my $apoc_subdir   = "apoc_pdb70_210721";
my $templ_lst_file  = "$proj_dir/bip198/pdb70_hhm/pdb70_210721.lst";
my $lib_des_file    = "$proj_dir/bip198/pdb70_hhm/pdb70_210721_name";


my $top = 1;
my $index;
my $help;

my $metric = 'score';    #### sort the hits using alignment TM-score

GetOptions(
            'a|aln_root=s' => \$apoc_aln_root,
            't|top=i'      => \$top,
            'h|help'       => \$help
          );

printHelp() if( $help );

my $lst_file = $ARGV[0];

unless( defined $lst_file ) {
  print "Error: a list file must be supplied and exists!";
  printHelp();
  exit 1;
}

my %lst = read_set_file( $lst_file );

my ($base, $ext) = fileparse( $lst_file, qr/\.\w*/ );
my $out_file = "$base\_${apoc_subdir}_top$top.dat";


#my %lib = read_set_file( $templ_lst_file );
#my %lib_des = read_lib_des_file( $lib_des_file );

my %top_aln = (); my $header;
my $total = scalar keys %lst;
my $counter = 0;
for my $query (sort keys %lst) {
  showProgress( ++$counter, $total, $query, 10, 25, '=' );

  my $aln_path = "$apoc_aln_root/$query/$apoc_subdir";
  my $aln_file = "$aln_path/$query\_aln.dat.gz";
  my $sco_file = "$aln_path/$query\_sco.dat.gz";

  unless( -s $sco_file ) {
    print "Error: $sco_file not found\n";
    next;
  }

  $header = `zless $sco_file |grep -P "^tname"`;


  my $hits = read_dt_file( $sco_file );
  my $num_hits = scalar @$hits;

  my @sorted_hits = sort{$b->{$metric}<=>$a->{$metric}} @$hits;
  my $qlen = $lst{$query};
  my $rank = 1;
  for (@$hits) {
    my $ln = $_->{'ln'};
    $rank++;

    #my $des = $pdb_des{$templ};
    #print $out_fh "$query\t$qlen\t$ln\n";
    $top_aln{$query} = {'rank'=>$rank, 'ln'=>$_->{'ln'}, 'metric'=>$_->{$metric}};
    last if $rank > $top;
  }

}

print "Writting the results...";
open my $out_fh, ">$out_file";
$counter = 1;
for my $query (sort {$top_aln{$b}->{'metric'}<=>$top_aln{$a}->{'metric'}} keys %top_aln) {
  my $entry = $top_aln{$query};
  my $qlen  = $lst{$query};

  if( $counter == 1 ) {
    print $out_fh "Query       \tQlen\t$header";
  }
  print $out_fh "$query\t$qlen\t$entry->{'ln'}\n";
  $counter++;
}
close $out_fh;
print "\n";


#####################################################################################
############################         Subroutine        ##############################
#####################################################################################
sub read_set_file {
  my $file = shift @_;

  my %set = ();
  open my $fh, "<$file" or die "Error: could not open $file\n";
  while(<$fh>){
    chomp;
    next if /name/;
    next if /^#/;

    my @fields = split(' ',$_);
    my $name = $fields[0];
    my $len  = $fields[1];
    $set{$name} = $len;
  }
  close $fh;

  return %set;
}

##### read fasta sequence ######
sub read_seq_file {
  my $file = shift @_;

  my %seq = ();
  open my $inp_fh, "<$file" or die "Error: could not open seq file $file\n";
  local $/ = "\>";
  while(<$inp_fh>) {
    chomp;
    next unless /\S+/;

    my @lines = split(/\n/, $_);
    my @cols  = split(' ', shift @lines );
    my $name   = $cols[0];
    #my $length = $cols[1];

    my $seq1 = join( ''  , @lines );

    ######### Check presence of non-standard amino acids in sequence
    if ($seq1 =~ /[BJOUXZ]/) {
      print "$name Warning: non-standard amino acids detected\n";
      next;
    }

    my $len = length( $seq1 );
    $seq{$name} = $seq1;
  }

  return %seq;
}


sub read_dt_file {
  my $file = shift @_;

  my @set = ();
  my $fh;
  if ( $file =~ /\.gz$/ ) {
    open $fh, "gunzip -c $file |" or die "Error: can’t open pipe to $file\n";
  }
  else {
    open $fh, "<$file" or die "Error: could not open seq file $file\n";
  }

  while(my $ln = <$fh>){
    chomp $ln;
    next if $ln =~ /tname/;
    next if $ln =~ /^#/;

    my @fields  = split(' ', $ln);
    my $templ   = $fields[0];
    my $alnlen  = $fields[1];
    my $rmsd    = $fields[2];
    my $score   = $fields[3];
    my $seqid   = $fields[4];

    $seqid =~ s/\%//;
    #next unless( exists $lib{$templ} );  ### only consider selected templates

    push( @set, { 'score'=>$score, 'alnlen'=>$alnlen, 'rmsd'=>$rmsd,
                  'templ'=>$templ, 'seqid'=>$seqid, 'ln'=>$ln } );
  }
  close $fh;

  return \@set;
}

#####################################################
sub read_aln_file {
  my $file = shift @_;

  my @set = ();
  my $fh;
  if ( $file =~ /\.gz$/ ) {
    open $fh, "gunzip -c $file |" or die "Error: can’t open pipe to $file\n";
  }
  else {
    open $fh, "<$file" or die "Error: could not open seq file $file\n";
  }

  local $/ = "### Alignment";
  while(<$fh>){
    chomp;
    next if /^#/;
    next unless /\S+/;

    my @lines = split(/\n/, $_);
    my $header = shift @lines;
    my $aln = join("\n", @lines);
    if( $header =~ /(\d+)\s+to:\s+(\S+.*)$/ ) {
      my @fields = split(' ', $2);
      my $rank   = $1;
      my $templ  = $fields[0];
      push( @set, { 'rank'=>$rank, 'templ'=>$templ, 'aln'=>$aln } );
    }
    else {
      print "Warning: template name not found, ignored\n";
      next;
    }
  }
  close $fh;

  return @set;
}

#####################################################
sub read_lst_file {
  my $file = shift @_;

  my @set = ();
  open my $fh, "<$file" or die "Error: could not open $file\n";
  while(<$fh>){
    chomp;
    next if /name/;
    next if /^#/;

    my @fields = split(' ',$_);
    my $name = $fields[0];
    push( @set, $name );
  }
  close $fh;

  return @set;
}

#####################################################
sub read_lib_des_file {
  my $file = shift @_;

  my %set = ();
  open my $fh, "<$file" or die "Error: could not open $file\n";
  while(<$fh>){
    chomp;
    next if /^#/;

    my @fields = split(' ',$_);
    my $name = $fields[1];
    my $num = scalar @fields;
    next unless( $fields[0] eq 'NAME');
    my $des = '';
    for(my $i=2; $i<$num; $i++) {
      $des .= " $fields[$i]";
    }

    $set{$name} = $des;
  }
  close $fh;

  return %set;
}




sub printHelp {

  print "Combining apoc alignments from split runs, usage:\n\n";
  print "\t$0 <options> <list file>\n\n";
  print "Options: -a <apoc alignment root directory> \n";
  print "         -n <top N alignments to make models from>\n";
  print "         -h print this help message\n";
  exit;
}


######################################################
# progress bar
######################################################

sub showProgress {
    my ( $got, $total, $name, $name_len, $width, $char ) = @_;
    $width ||= 25;
    $char  ||= '=';
    $name_len ||= '4';
    my $num_width = length $total;
    local $| = 1;


    printf "|%-${width}s| %-${name_len}s %${num_width}s of %s (%.2f%%)\r",
        $char x (($width-1)*$got/$total). '>', $name, $got, $total, 100*$got/$total;

    if( $got == $total ) { print "\n"; }
}
