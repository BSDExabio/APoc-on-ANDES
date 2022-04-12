#! /usr/bin/env perl
## Input: alignment from apoc runs on splits of the template set
## Output: combine them into a single output file

use warnings;
use strict;
use File::Path;
use File::Basename;
#use Compress::Zlib;
use Getopt::Long;
#use Statistics::Descriptive;


my $home_dir = $ENV{"HOME"};

my $apoc_aln_root = "/gpfs/alpine/world-shared/bif135/species/desulfovibrio/apoc/apoc_pdb70_210721";
if( defined $ARGV[1] ) {
  $apoc_aln_root = $ARGV[1];
}


my $templ_lst_file = "/gpfs/alpine/proj-shared/bip198/pdb70_hhm/pdb70_210721.lst";
my $lib_des_file   = "/gpfs/alpine/proj-shared/bip198/pdb70_hhm/pdb70_210721_name";

my $max_num_hits = 1000;

my $top = 100;
my $index;
my $help;

my $metric = 'score';   #### sort the hits using predicted TM-score
#my $metric = 'score';    #### sort the hits using alignment score

GetOptions(
            'a|aln_root=s' => \$apoc_aln_root,
            'n|num=i'      => \$top,   ### number of top entries to print alignment out
            'h|help'       => \$help
          );

printHelp() if( $help );

my $lst_file = $ARGV[0];

unless( defined $lst_file ) {
  print "Error: a list file must be supplied and exists!";
  printHelp();
  exit 1;
}

my @lst = ();
if( -s $lst_file ) {
  @lst = read_lst_file( $lst_file );
}
else {
  push( @lst, $ARGV[0] );
}

my %lib = read_set_file( $templ_lst_file );
my %lib_des = read_lib_des_file( $lib_des_file );

my $total = scalar @lst;
my $counter = 0;
for my $entry (@lst) {
  my $query = $entry;
  showProgress( ++$counter, $total, "$query", 10, 25, '=' );

  #my $aln_path = "$apoc_aln_root/$query";
  my $aln_path = "$apoc_aln_root";
  my @sa_sco_lst = <$aln_path/$query\_apoc_[0-9]*.gz>;

  unless( scalar @sa_sco_lst > 0 ) {
    print "Info: no score file found, skipped\n";
    next;
  };

  my @score = (); my @sa_data = (); my @aln_data = ();
  my $header = '';
  for my $sco_file (@sa_sco_lst) {
    #print "Info: processing $sco_file\n";
    $header = `zless $sco_file |grep -P "^(tname)"`;
    my $hits = read_dt_file( $sco_file );
    for (@$hits) {
      push( @score, $_->{'score'} );
      push( @sa_data, $_ );
    }

    my $index = '';
    $index = $1 if ( $sco_file =~ /$query\_apoc_([0-9]*).gz/ );
    my $aln_file = "$aln_path/$query\_aln_$index.gz";

    next unless( -s $aln_file );
    my @aln_lst = read_aln_file( $aln_file );
    for (@aln_lst) {
      push( @aln_data, $_ );
    }
  }

  #my $stat = Statistics::Descriptive::Full->new();
  #$stat->add_data( \@score );
  #my $mean = $stat->mean();
  #my $sd = $stat->standard_deviation();

  my $out_dir = "$apoc_aln_root";
  mkpath $out_dir unless( -d $out_dir );

  my $out_score_file = "$out_dir/$query\_sco.dat.gz";
  open my $out_score_fh, "| /bin/gzip -c >$out_score_file" or die "Error: could not write $out_score_file\n";
  chomp $header;
  print $out_score_fh "$header\tDescription\n";
  my %top_aln = ();
  my $rank = 1;
  for my $entry (sort {$b->{$metric}<=>$a->{$metric}} @sa_data) {
    #my $zscore = ($entry->{'score'} - $mean) / $sd;
    my $templ  = $entry->{'templ'};
    my $des    = $lib_des{$templ};
    printf $out_score_fh "%-10s %-4d %s %s %5.1f%%\t|%s\n", $templ,
      $entry->{'alnlen'}, $entry->{'rmsd'}, $entry->{'score'}, $entry->{'seqid'}*100, $des;
    if( $rank <= $top ) {
      $top_aln{$templ} = {};
      for my $attr (keys %$entry) {
         $top_aln{$templ}->{$attr} = $entry->{$attr};
      }
      $top_aln{$templ}->{'rank'}   = $rank;
      #$top_aln{$templ}->{'zscore'} = $zscore;
    }
    $rank++;
    last if( $rank == $max_num_hits );
  }
  close $out_score_fh;

  #### find top ranked alignments to print out
  next if( $top == 0 );   ### no alignment requested
  next if( scalar @aln_data == 0 ); ### no alignment
  for my $entry (@aln_data) {
    my $templ = $entry->{'templ'};
    if( exists $top_aln{$templ} ) {
      #$top_aln{$templ}->{'predtms'} = $entry->{'predtms'};
      $top_aln{$templ}->{'aln'} = $entry->{'aln'};
      $top_aln{$templ}->{'trans'} = $entry->{'trans'};
      $top_aln{$templ}->{'rot'} = $entry->{'rot'};
    }
  }
  my $out_aln_file = "$out_dir/$query\_aln.dat.gz";
  open my $out_aln_fh, "| /bin/gzip -c >$out_aln_file" or die "Error: could not write $out_aln_file\n";
  for my $templ (sort {$top_aln{$a}->{'rank'}<=>$top_aln{$b}->{'rank'}} keys %top_aln) {
    my $ent = $top_aln{$templ};
    printf $out_aln_fh "### Alignment %d to: %s naln=%d score=%.4f sid=%.1f%%\n",
      $ent->{'rank'}, $templ, $ent->{'alnlen'}, $ent->{'score'}, $ent->{'seqid'}*100;
    print $out_aln_fh $ent->{'trans'}, "\n";
    print $out_aln_fh $ent->{'rot'}, "\n";
    print $out_aln_fh $ent->{'aln'}, "\n";
  }
  close $out_aln_fh;
}



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
    next unless( exists $lib{$templ} );  ### only consider selected templates

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

  local $/ = "### ";
  while(<$fh>){
    chomp;
    next if /^#/;
    next unless /\S+/;

    my @lines = split(/\n/, $_);
    my $header    = shift @lines;
    my $trans_vec = shift @lines;
    my $rot_mat   = shift @lines;

    my $aln = join("\n", @lines);
    if( $header =~ /(\S+)\s+(\S+.*)$/ ) {
      my @fields = split(' ', $2);
      my $templ  = $1;
      push( @set, { 'templ'=>$templ, 'aln'=>$aln, 'trans'=>$trans_vec, 'rot'=>$rot_mat } );
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
