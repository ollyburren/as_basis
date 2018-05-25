use strict;
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;


## we can speed things up if we presort the files
#OUTDIR=/home/ob219/scratch/as_basis/gwas_stats/processed_sorted/
#for f in `\ls /home/ob219/scratch/as_basis/gwas_stats/processed_new/*.tab`;do
#  outfile=$OUTDIR$(basename $f)
# echo "(head -n 1 $f && tail -n +2 $f | sort -k1) > $outfile"
#done

my $COMM_BIN = '/usr/bin/comm';

my $IN_DIR = "/home/ob219/scratch/as_basis/gwas_stats/processed_sorted";
my $MANIFEST = "/home/ob219/git/as_basis/manifest/as_manifest_feb_2018.csv";

# read in manifest

open(MAN,$MANIFEST) || die "Cannot open $MANIFEST\n";
my @mhead;
my @man;
while(<MAN>){
  #print "$_\n";
  chomp;
  if(!@mhead){
    @mhead = split(",",$_);
    next;
  }
  my %rhash=();
  my @tmp = split(",",$_);
  for(my $i=0;$i<@tmp;$i++){
    $rhash{$mhead[$i]} = $tmp[$i];
  }
  push @man,\%rhash if $rhash{'include'} eq 'Y';
}
close(MAN);

my @tfiles;
for (my $i=0,;$i<2;$i++){
  my $tmp = File::Temp->new(UNLINK=>0);
  $tfiles[$i] = $tmp->filename();
  $tmp->close();
}

#if not sorted do this but slower
my $cmd_template1 = '%s -12 <(cut -f1 %s | sort) <(cut -f1 %s | sort) > %s';
my $cmd_template2 = '%s -12 <(cut -f1 %s | sort) <(cat %s) > %s';
#my $cmd_template1 = '%s -12 <(tail -n +2 %s | cut -f1) <(tail -n +2 %s | cut -f1) > %s';
#my $cmd_template2 = '%s -12 <(tail -n +2 %s | cut -f1) <(cat %s) > %s';
my @files;
push @files,"$IN_DIR/".$_->{'file'} foreach @man;
my $cmd = sprintf($cmd_template1,($COMM_BIN,$files[0],$files[1],$tfiles[1]));
print "echo \"$cmd\"\n";
#print "$cmd\n";
my $ofile;
for(my $i=2;$i<@files;$i++){
  $cmd = sprintf($cmd_template2,($COMM_BIN,$files[$i],$tfiles[($i % 2)-1],$tfiles[($i % 2)]));
  print "echo \"$cmd\"\n";
  #print "$cmd\n";
  #system($cmd);
  $ofile = $tfiles[($i % 2)];
}
print "cp $ofile ~/tmp/intersect.txt\n";
print "rm $_\n" foreach @tfiles;
