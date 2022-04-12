#! /usr/bin/env perl
## use APoc to perform global structure alignment against a pdb library, i.e., templates.
## output scores and alignment
## use mutliprocessing to take advantage of multiple cores

use warnings;
use strict;
use File::Path;
use File::Basename;
use Getopt::Long;
use Benchmark;
use POSIX ":sys_wait_h";

### defining variable namespace to be overwritten by command line arguments if provided
my $prog;
my $seq_flag;
my $mea_flag=1;

my $prealn=1;
my $help;
my $num_cpu = 1;

my $home_dir = $ENV{'HOME'};
my $hostname = $ENV{'HOSTNAME'};
my $que_dat_root = "/gpfs/alpine/world-shared/bif135/species/desulfovibrio/af2c_mono";
my $out_root     = "/gpfs/alpine/world-shared/bif135/species/desulfovibrio/apoc";

################################################################################################
my $tem_dat_root = "/gpfs/alpine/proj-shared/bip198/pdb70_hhm";
my $tem_lst_file = "pdb70_210721_block.lst";
#my $tem_lst_file = "block_test.lst";
my $tem_blk_file = "pdb70_210721_block.pdb";

my $bin_path   = "/gpfs/alpine/bip198/proj-shared/rbd_work/APoc/apoc/bin";
#my $combine_apoc_output = "/ccs/home/mugao/projects/apoc/scripts/combine_apoc_results_andes.pl";
#my $splitset = "splitset.pl";
my $combine_apoc_output = "/gpfs/alpine/bip198/proj-shared/rbd_work/APoc/scripts/combine_apoc_results_andes.pl";
my $splitset = "/gpfs/alpine/bip198/proj-shared/rbd_work/APoc/scripts/splitset.pl";

#################################################################################################
GetOptions('p|program=s' => \$prog,
           's|seq'       => \$seq_flag,
           'm|measure=i' => \$mea_flag,
           'ps|prealn=i' => \$prealn,
           'c|num_cpu=i' => \$num_cpu,
           'out_root=s'  => \$out_root,
           'dat_root=s'  => \$que_dat_root,
           'h|help'      => \$help);

my $SWITCH = "-v 2";
my $measure;
my $mode;

if( defined $seq_flag ) {
  $SWITCH .= ' -sod';
  $mode = 'seq';
  $measure = 'tms';
}
else {
  $measure = 'pss';
  $mode = 'nseq';
}

if( $prealn==0 ) {
  $SWITCH .= " -fa 0";   #### no global alignment
}

if( not defined $prog ) {
  $prog = 'apoc';
}


my @query_lst = ();
if( -s $ARGV[0] ) {
  @query_lst = read_inp_file( $ARGV[0] );
}
else {
  for (@ARGV) { push( @query_lst, $_ ); }
}



###### randomly generate a work path ######
my $user = `whoami`;
chomp $user;
my $rannum = int(rand(1000000));
#my $work_path = "/tmp/$user/aln-$rannum";
my $work_path = "/dev/shm/$user/aln-$rannum";

##########################################


mkpath $work_path unless( -d $work_path );
chdir $work_path or die "Error: could not enter $work_path";

`cp $tem_dat_root/$tem_blk_file .`;
`cp $tem_dat_root/$tem_lst_file . && $splitset $tem_lst_file $num_cpu`;



############################################
for my $query (@query_lst) {
  print "Info: working on $query at $hostname\n";
  #my $query_pdb_file = "$que_dat_root/$query.pdb";
  my $query_pdb_file = "$que_dat_root/$query/ranked_0.pdb";   ### top 1 AlphaFold2 model name

  ############## putput path  ########
  my $out_path = "$out_root/$query/apoc_pdb70_210721";
  mkpath $out_path unless ( -d $out_path );

  my $forks = 0;
  for(my $i=1; $i<=$num_cpu; $i++) {
    my $pid = fork;
    if( not defined $pid ) {
      die "fork failed: $!\n";
    }
    if ($pid) {
      $forks++;
    }
    else {
      my $lst_file = "$tem_lst_file.$i";
      die unless ( -s $lst_file );

      my $start_time = Benchmark->new;
      my $index = sprintf "%03d", $i;
      my $score_file = "$query\_apoc\_$index";   ### record APoc scores
      #my $trans_file = "$query\_trans\_$index";  ### record the best transformation matrix
      my $align_file = "$query\_aln\_$index";    ### record alignment

      ###--------------------->   start APoc scanning template library   <-----------------------###
      my $aln_output = `$bin_path/$prog -block $tem_blk_file -lt $lst_file $query_pdb_file $SWITCH`;
      my ( $glb_aln_out, $pkt_aln_out ) = readPKAlnOut( $aln_output );

      open  OUT, ">$score_file" or die "Error: could not write $score_file\n";
      print OUT "tname  malnlen  mrmsd   mscore  mseqid\n";

      #open  TRA, ">$trans_file" or die "Error: could not write $trans_file\n";
      #print TRA "tname\tmscore\ttranslation_vec\trotation_mat\n";

      open  ALN, ">$align_file" or die "Error: could not write $align_file\n";

      for my $mout (@$glb_aln_out) {
        my $mstruct1 = $mout->{struct1};
        my $mstruct2 = $mout->{struct2};

        my ($mbase1, $mpath1, $mext1) = fileparse( $mstruct1, qr/\.\w*$/ );
        printf OUT "%-8s %5d %7.3f %8.5f %6.3f\n", $mbase1, $mout->{alnlen}, $mout->{rmsd}, $mout->{score}, $mout->{sid};

        my $trans_vec = $mout->{'trans_vec'};
        my $rot_mat   = $mout->{'rot_mat'};
        #printf TRA "%s\t%7.3f\t%s\t%s\n", $mbase1, $mout->{score}, join( ',', @$trans_vec ),
        #  join( ',', @{$$rot_mat[0]}, @{$$rot_mat[1]}, @{$$rot_mat[2]} );

        my $aln = $mout->{'aln'};
        printf ALN "### %-8s %5d %7.3f %8.5f %6.3f\n", $mbase1, $mout->{alnlen}, $mout->{rmsd}, $mout->{score}, $mout->{sid};
        printf ALN "%s %s\n", 'Translation:', join( ',', @$trans_vec );
        printf ALN "%s %s\n", 'Rotation:   ', join( ',', @{$$rot_mat[0]}, @{$$rot_mat[1]}, @{$$rot_mat[2]} );
        print ALN "Index Ch1 Resid1  AA1 Ch2 Resid2  AA2 Distance Cos(theta)\n";
        for (@$aln) {
          print ALN $_->{'ln'}, "\n";
        }
      }
      ###--------------------->   end of scanning template library   <---------------------------###

      my $end_time = Benchmark->new;
      my $tot_td   = timediff($end_time, $start_time);
      if( timestr($tot_td) =~ /=\s+(\S+)\s+CPU\)/ ) {
        my $tot_time = $1;
        printf OUT "### Total run time %.3e seconds\n", $tot_time;
      }

      close OUT;
      #close TRA;
      close ALN;

      #`gzip $score_file && mv $score_file.gz $out_path/`;
      #`gzip $align_file && mv $align_file.gz $out_path/`;
      `gzip -f $score_file`;
      `gzip -f $align_file`;

      print "Info: finished aligning $query against $lst_file at $hostname\n";
      exit;
    }
  }

  for (1 .. $forks) {
     my $pid = wait();
  }

  print "Info: combining results of $query ...\n";
  `echo $query > test.lst`;
  my $out = `$combine_apoc_output test.lst .`;
  #print "$out";
  `mv $query\_sco.dat.gz $out_path/`;
  `mv $query\_aln.dat.gz $out_path/`;
  print "Info: combining results of $query done\n";
}

`rm -rf $work_path` if ( -d $work_path );

print "Info: waiting for all child process gracefully exit at $hostname...\n";
1 while waitpid(-1, WNOHANG) > 0;  ### avoid potential segmentation fault caused by child process clean-up not done.

exit;





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

    push( @set, $name )
  }
  close $fh;

  return @set;
}

sub readListFile {
  my $file = shift @_;

  my %set = ();
  open my $fh, "<$file" or die "Error: could not open $file\n";
  while(<$fh>){
    chomp;
    next if /name/;
    next if /^#/;

    my @fields = split(' ',$_);
    my $pdb_file = $fields[0];
    my $pocket   = $fields[1];
    my $blk_file = $fields[2];

    my @cols = split(/\./, $pocket);
    my $mono = $cols[0];

    $set{$pdb_file} = { 'pkt'=>$pocket, 'blk'=>$blk_file };
  }
  close $fh;

  return %set;
}

sub  mk_pkt_file {
  my ( $query_blk, $query_pkt, $query ) = @_;

  open my $inp_fh, "<$query_blk" or die "Error: could not open file $query_blk\n";
  open my $out_fh, ">$query.pdb" or die "Error: could not write $query.pdb\n";
  my $write_flag = 0;
  while(<$inp_fh>) {

    if( /^PDBSTA\s+$query/ ) {
      $write_flag = 1;
    }
    elsif( /^PDBEND/ ) {
      $write_flag = 0;
    }

    print $out_fh $_ if( $write_flag );
  }
  close $out_fh;
  close $inp_fh;
}


sub readPKAlnOut {
  my $aln_out = shift @_;

  my @glb_aln_out = ();
  my @pkt_aln_out = ();

  my $sep = ">>>>>>>>>>>>>>>>>>>>>>>>>";
  my @rec = split(/$sep/, $aln_out);
  for my $record (@rec) {
    if ( $record =~ /Global alignment/ ) {
      my @lines  = split(/\n/, $record );
      my $result = readMyAlnOut( @lines );
      push( @glb_aln_out, $result );
    }
    elsif ( $record =~ /Pocket alignment/ ) {
      my @lines  = split(/\n/, $record );
      my $result = readMyAlnOut( @lines );
      push( @pkt_aln_out, $result );
    }
  }

  return ( \@glb_aln_out, \@pkt_aln_out );
}



#### read myalign output ####
sub readMyAlnOut {
  my @lines = @_;

  my ($struct1, $struct2, $len1, $len2, $pk1, $pk2);
  my ($alnlen, $rmsd, $score, $sid, $best_init, $pvalue, $zscore);
  my ($tot_seqsc, $tot_seqsc_pos, $run_time);
  my $num_sim_aa = 0;
  my @aln = ();
  my @trans_vec = ();
  my @rot_mat = ();

  my $read_flag = 0;

  my $num = scalar @lines;
  for(my $i=0; $i<$num; $i++) {
    my $line = $lines[$i];
    #print "$line\n";
    if( $line =~ /Structure 1:\s*(\S+).*Length =\s*(\d+) AAs, Full/ ) {
      $struct1 = $1;
      $len1 = $2;
    }
    elsif ( $line =~ /Structure 2:\s*(\S+).*Length =\s*(\d+) AAs, Full/ ) {
      $struct2 = $1;
      $len2 = $2;
    }
    elsif ( $line =~ /Structure 1:\s*(\S+).*Length =\s*(\d+) AAs, Pocket:(\S+)/ ) {
      $struct1 = $1;
      $len1 = $2;
      $pk1  = $3;
    }
    elsif ( $line =~ /Structure 2:\s*(\S+).*Length =\s*(\d+) AAs, Pocket:(\S+)/ ) {
      $struct2 = $1;
      $len2 = $2;
      $pk2  = $3;
    }
    elsif ( $line =~ /^Aligned length=\s*(\d+), RMSD=\s*(\d+\.\d+), \w{2}-score=(\d+\.\d+), ID=(\d+\.\d+)/ ) {
      $alnlen = $1;
      $rmsd   = $2;
      $score  = $3;
      $sid    = $4;
    }
    elsif ( $line =~ /^\w{2}-score =\s*(\-*\d+\.\d+), P-value = \s*(\S+), Z-score =\s*(\S+)/ ) {
      $score  = $1;
      $pvalue = $2;
      $zscore = $3;
    }
    elsif ( $line =~ /^\w{2}-score =\s*(\d+\.\d+)/ ) {
      $score    = $1;
    }
    elsif ( $line =~ /^Number of aligned residues  =\s*(\d+)/ ) {
      $alnlen = $1;
    }
    elsif ( $line =~ /^RMSD =\s*(\d+\.\d+), Seq identity  =\s*(\d+\.\d+)/ ) {
      $rmsd = $1;
      $sid  = $2;
    }
    elsif( $line =~ /Best alignment search: initial =\s*(\d+)/ ) { $best_init = $1;}

    elsif( $read_flag and $line =~ /^\s*\d+/ ) {
      my $ch1   = substr($line,9, 1);
      (my $res1 = substr($line,11,4)) =~ s/\s//g;
      (my $ic1  = substr($line,15,1)) =~ s/\s//g;
      my $aa1   = substr($line,20,1);
      my $ch2   = substr($line,25,1);
      (my $res2 = substr($line,27,4)) =~ s/\s//g;
      (my $ic2  = substr($line,31,1)) =~ s/\s//g;
      my $aa2   = substr($line,36,1);
      (my $dist  = substr($line,39,6))=~ s/\s//g;
      (my $note = substr($line,57,4))=~ s/\s//g;

      #if( $dist < 5 ) {
      my $ares1 = "$ch1.$res1$ic1.$aa1";
      my $ares2 = "$ch2.$res2$ic2.$aa2";

      ### excluding redundant residues, deal with an issue in Jeff's pdb structures
      my $last = $aln[-1];
      next if( defined $last and ($ares1 eq $last->{'1'} or $ares2 eq $last->{'2'}) );

      push( @aln, { '2'=>$ares2, '1'=>$ares1, 'ln'=>$line } );
      $num_sim_aa++ if( $note eq '*' or $note eq ':' );
      #print "$ares1 -- $ares2 -- $dist\n";
      #}
    }
    elsif( $line =~ /rotation matrix/ ) {
      $i++;
      for(my $j=0; $j<3; $j++) {
         $i++;
         my @fields = split(' ',$lines[$i]);
         push( @trans_vec, $fields[1] );
         push( @rot_mat, [ ($fields[2], $fields[3], $fields[4]) ] );
      }
    }
    elsif( $line =~ /Sequence similarity.*\: sum =\s*(\-?\d+),\s+sum_pos =\s*(\d+)/ ) {
      $tot_seqsc = $1;
      $tot_seqsc_pos = $2;
    }
    elsif( $line =~ /Running time:\s*(\S+)\s+seconds/ ) {
      $run_time = $1;
    }

    $read_flag =1 if ( $line =~ /^ Index Ch1 Resid1/ );
  }

  #### transformation from struct 2 to struct 1
  my @rev_trans_vec = ();
  my @rev_rot_mat = ();

  if( scalar @rev_trans_vec and scalar @rev_rot_mat ) {
    for(my $i=0; $i<3; $i++) {
      $rev_rot_mat[$i] = [];
      for( my $j=0; $j<3; $j++) {
         $rev_rot_mat[$i][$j] = $rot_mat[$j][$i];    ##### unitary matrix
      }
    }

    my $x = -$trans_vec[0];
    my $y = -$trans_vec[1];
    my $z = -$trans_vec[2];

    $rev_trans_vec[0] = $x*$rev_rot_mat[0][0] + $y*$rev_rot_mat[0][1] + $z*$rev_rot_mat[0][2];
    $rev_trans_vec[1] = $x*$rev_rot_mat[1][0] + $y*$rev_rot_mat[1][1] + $z*$rev_rot_mat[1][2];
    $rev_trans_vec[2] = $x*$rev_rot_mat[2][0] + $y*$rev_rot_mat[2][1] + $z*$rev_rot_mat[2][2];
  }

  my %result = ( 'pv'=>$pvalue, 'zs'=>$zscore, 'alnlen'=>$alnlen, 'rmsd'=>$rmsd, 'score'=>$score,
     'sid'=>$sid, 'best_init'=>$best_init, 'aln'=>\@aln, 'rot_mat'=>\@rot_mat, 'trans_vec'=>\@trans_vec,
     'seqsc'=>$tot_seqsc, 'seqsc_pos'=>$tot_seqsc_pos, 'struct1'=>$struct1, 'struct2'=>$struct2,
     'len1'=>$len1, 'len2'=>$len2, 'pkt1'=>$pk1, 'pkt2'=>$pk2, 'run_time'=>$run_time, 'numsimaa'=>$num_sim_aa,
         );

  return \%result;
}
