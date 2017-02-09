#!/usr/bin/perl

# CHSL Hackathon
# Oct 24th, 2016

use lib "./";
use strict;
use warnings;
use feature 'state';
use Fcntl;
use Getopt::Std;
use Data::Dumper;
use File::Basename;
use Cwd 'abs_path';
use Math::Random qw(random_poisson random_uniform_integer);

# Check dependencies
my $absPath = dirname(abs_path($0));
die "DWGSIM executable not found\n" if  (!-e "$absPath/dwgsim");
die "SURVIVOR executable not found\n" if  (!-e "$absPath/SURVIVOR");
die "SURVIVOR parameter list not found\n" if  (!-e "$absPath/parameter");
die "bfr executable not found\n" if  (!-e "$absPath/bfr");
die "pigz executable not found\n" if  (!-e "$absPath/pigz");
die "samtools executable not found\n" if  (!-e "$absPath/samtools");
# Check dependencies end

&main;
0;

# atexit
my @fileHandlersAtExit = ();
END {
  foreach(@fileHandlersAtExit)
  {
    close $_;
  }
}
# atexit end

sub main
{
  &usage if (@ARGV < 1);
  my %opts = (h=>undef, d=>2, r=>undef, p=>undef, b=>undef,
              e=>0.01, E=>0.01, i=>350, s=>35, x=>1000, f=>50, t=>1500, m=>10 );
  getopts('hd:r:p:b:e:E:i:s:x:f:t:m:', \%opts);
  &usage if (defined $opts{h});

  #Check options
  die "Number of haplotypes should be between 1 to 3\n" if ($opts{d} < 1 || $opts{d} > 3); #3 is a soft limit
  die "Please provide a reference genome with -r\n" if (not defined $opts{r});
  die "Please provide a output prefix with -p\n" if (not defined $opts{p});
  die "Please provide a barcodes file with -b\n" if (not defined $opts{b});
  die "Reference genome $opts{r} not exist\n" if (!-e "$opts{r}");
  die "Reference genome index $opts{r}.fai not exist\n" if (!-e "$opts{r}.fai");
  die "Barcodes file $opts{b} not exist\n" if (!-e "$opts{b}");
  die "Please provide a output prefix for this read simulation job with -p\n" if (not defined $opts{p});
  die "The value of -e should be set between 0 and 1\n" if ( $opts{e} < 0 || $opts{e} > 1 );
  die "The value of -E should be set between 0 and 1\n" if ( $opts{E} < 0 || $opts{E} > 1 );
  die "The value of -i should be set between 350 and 400\n" if ( $opts{i} < 350 || $opts{i} > 400 );
  die "The value of -s should be set between 35 and 40\n" if ( $opts{s} < 35 || $opts{s} > 40 );
  #die "The value of -x should be set between 800 and 1200\n" if ( $opts{x} < 800 || $opts{x} >1200 );
  #die "The value of -f should be set between 50 and 100\n" if ( $opts{f} < 50 || $opts{f} >100 );
  #die "The value of -t should be set between 1300 and 1700\n" if ( $opts{t} < 1300 || $opts{t} >1700 );
  #die "The value of -m should be set between 5 and 15\n" if ( $opts{m} < 5 || $opts{m} >15 );
  warn "$opts{p}.status exists\n" if (-e "$opts{p}.status");
  #Check options end

  #Global variables
  &Log("$opts{p}.status"); #Initializing Log routine
  my $arrayLimit = 2147483648;
  my $readLenghtWithBarcode = 151;
  my $genomeSize = 0;
  my %faidx = ();
  my @boundary = ();
  my %readPositionsInFile = ();
  my @barcodes = ();
  my $numBarcodes = 0;
  #Global variables end

  #Generate copies of haplotypes
  {
    for(my $i = 0; $i < $opts{d}; ++$i)
    {
      # SURVIVOR command
      # ./SURVIVOR 1 ref.fa parameter 0 ./alt.fa
      if(-e "$opts{p}.survivor.$i.fasta" && -e "$opts{p}.survivor.$i.insertions.fa" && -e "$opts{p}.survivor.$i.bed")
      { &Log("SURVIVOR round $i done already"); next; }
      &Log("SURVIVOR round $i start");
      &Log("Running: $absPath/SURVIVOR 1 $opts{r} parameter 0 $opts{p}.survivor.$i");
      system("$absPath/SURVIVOR 1 $opts{r} parameter 0 $opts{p}.survivor.$i 1>/dev/null");
      if(!-e "$opts{p}.survivor.$i.fasta")
      { &LogAndDie("SURVIVOR round $i error on missing $opts{p}.survivor.$i.fasta"); }
      if(!-e "$opts{p}.survivor.$i.insertions.fa")
      { &LogAndDie("SURVIVOR round $i error on missing $opts{p}.survivor.$i.insertions.fa"); }
      if(!-e "$opts{p}.survivor.$i.bed")
      { &LogAndDie("SURVIVOR round $i error on missing $opts{p}.survivor.$i.bed"); }
      &Log("SURVIVOR round $i end");
    }
  }
  #Generate copies of haplotypes end

  #Generate reads for haplotypes
  {
    my $readsPerHaplotype = int($opts{x} * 1000000 / $opts{d}) * 2;
    for(my $i = 0; $i < $opts{d}; ++$i)
    {
      # dwgsim command
      # ./dwgsim -N 1000 -e 0.02 -E 0.02 -d 350 -s 35 -1 151 -2 151 -S 0 -c 0 ref.fa ./test
      if(-e "$opts{p}.dwgsim.$i.12.fastq")
      { &Log("DWGSIM round $i done already"); next; }
      &Log("DWGSIM round $i start");
      &Log("$absPath/dwgsim -N $readsPerHaplotype -e $opts{e} -E $opts{E} -d $opts{i} -s $opts{s} -1 $readLenghtWithBarcode -2 $readLenghtWithBarcode -S 0 -c 0 -m /dev/null $opts{p}.survivor.$i.fasta $opts{p}.dwgsim.$i");
      system("$absPath/dwgsim -N $readsPerHaplotype -e $opts{e} -E $opts{E} -d $opts{i} -s $opts{s} -1 $readLenghtWithBarcode -2 $readLenghtWithBarcode -S 0 -c 0 -m /dev/null $opts{p}.survivor.$i.fasta $opts{p}.dwgsim.$i");
      if(!-e "$opts{p}.dwgsim.$i.12.fastq")
      { &LogAndDie("DWGSIM round $i error on missing $opts{p}.dwgsim.$i.12.fastq"); }
      &Log("DWGSIM round $i end");
    }
  }
  #Generate reads for haplotypes end

  #Load reference genome index
  {
    &Log("Load faidx start");

    for(my $i = 0; $i < $opts{d}; ++$i)
    {
      &Log("Load faidx haplotype: $i");
      &Log("$absPath/samtools faidx $opts{p}.survivor.$i.fasta");
      system("$absPath/samtools faidx $opts{p}.survivor.$i.fasta");
      $genomeSize = &LoadFaidx(\%{$faidx{$i}}, \@boundary, "$opts{p}.survivor.$i.fasta");
      &LogAndDie("Failed loading genome index $opts{p}.survivor.$i.fasta.fai") if ($genomeSize == 0);
    }
    #print Dumper %faidx;
    &Log("Load faidx end");
  }
  #Load reference genome index end

  #Load read positions into memory
  {
    &Log("Load read positions start");
    for(my $i = 0; $i < $opts{d}; ++$i)
    {
      &Log("Load read positions haplotype $i");
      $readPositionsInFile{$i} = ();
      open my $fh, "$opts{p}.dwgsim.$i.12.fastq" or &LogAndDie("Error opening $opts{p}.dwgsim.$i.12.fastq");
      my $l1; my $l2; my $l3; my $l4; my $l5; my $l6; my $l7; my $l8;
      my $newFpos;
      my $fpos = tell($fh); &LogAndDie("Fail to tell file position") if $fpos == -1;
      while($l1=<$fh>)
      {
        $l2=<$fh>; $l3=<$fh>; $l4=<$fh>; $l5=<$fh>; $l6=<$fh>; $l7=<$fh>; $l8=<$fh>;
        $newFpos = tell($fh);
        #if(!eof($fh))
        #{
        #  my $c = getc($fh); seek($fh, $newFpos, 0);
        #  die "$c\n$l1$l2$l3$l4$l5$l6$l7$l8" if $c ne "@";
        #}
        #@chr1_111758675_111758819_0_1_0_0_2:0:0_3:0:0_0/1
        $l1=~/@(chr\d+)_(\d+?)_/;
        my $gCoord = &GenomeCoord2Idx(\%{$faidx{$i}}, "$1", $2);
        if($gCoord < 0 || $gCoord >= $genomeSize)
        { &LogAndDie("$1 $2 $gCoord $fpos"); }
        push @{$readPositionsInFile{$i}{$gCoord}}, $fpos;
        $fpos = $newFpos;
      }
      close $fh;
    }
    &Log("Load read positions end");
  }
  #Load read positions into memory end

  #Load barcodes
  {
    &Log("Load barcodes start");
    open my $fh, "$opts{b}" or &LogAndDie("Barcodes file $opts{b} not exist");
    @barcodes = <$fh>; chomp(@barcodes);
    $numBarcodes = scalar(@barcodes);
    close $fh;
    &Log("Load barcodes end");
  }
  #Load barcodes end

  #Simulate reads
  {
    &Log("Simulate reads start");
    open my $fq1Outputfh, "| $absPath/bfr -b 256M | $absPath/pigz -c > $opts{p}.10Xsimulate.1.fq.gz" or &LogAndDie("Error opening $opts{p}.10Xsimulate.1.fq.gz");
    open my $fq2Outputfh, "| $absPath/bfr -b 256M | $absPath/pigz -c > $opts{p}.10Xsimulate.2.fq.gz" or &LogAndDie("Error opening $opts{p}.10Xsimulate.2.fq.gz");

    # depthPerMol * molLength * #molPerPartition * Partitions = reads * length
    # ? * 50k * 10 * 1.5M = 1000M * 270
    # ? = 0.36x
    # readsPerParition = depthPerMol * molLength * #molPerPartition / length
    # ? = 0.36x * 50k * 10 / 270
    # ? = 666.6
    #
    my $readsPerMolecule = int(0.499 + ($opts{x} * 1000000) / ($opts{t} * 1000) / $opts{m});
    &Log("readsPerMolecule: $readsPerMolecule");

    # For every Haplotype
    for(my $i = 0; $i < $opts{d}; ++$i)
    {
      &Log("Simulating on haplotype: $i");
      my $readsCountDown = int($opts{x} * 1000000 / $opts{d});
      &Log("readsCountDown: $readsCountDown");
      open my $readsInputfh, "$opts{p}.dwgsim.$i.12.fastq" or &LogAndDie("Fail opening $opts{p}.dwgsim.$i.12.fastq");

      while($readsCountDown > 0)
      {
        #Pick a barcode
        my $selectedBarcode = "";
        while(1)
        {
          my $idx = int(rand($numBarcodes));
          next if $barcodes[$idx] eq "";
          $selectedBarcode = $barcodes[$idx];
          $barcodes[$idx] = "";
          last;
        }

        my $numberOfMolecules = &PoissonMoleculePerPartition($opts{m}-1)+1;
        #&Log("numberOfMolecules: $numberOfMolecules");
        my $readsToExtract = $numberOfMolecules * $readsPerMolecule;
        #&Log("readsToExtract: $readsToExtract");
        for(my $j = 0; $j < $numberOfMolecules; ++$j)
        {
          #Pick a starting position
          my $startingPosition = int(rand($genomeSize));
          #&Log("startingPosition: $startingPosition");
          #Pick a fragment size
          my $moleculeSize  = &PoissonMoleculeSize($opts{f}*1000);

          #Check and align to boundary
          my $lowerBoundary; my $upperBoundary;
          &bSearch($startingPosition, \@boundary, \$lowerBoundary, \$upperBoundary);
          if(($startingPosition + $moleculeSize) > $upperBoundary)
          {
            my $newMoleculeSize = $upperBoundary - $startingPosition;
            if($newMoleculeSize < 1000) #dirty thing
            {
              --$j;
              next;
            }
            $readsToExtract = int($readsToExtract * $newMoleculeSize / $moleculeSize);
            $moleculeSize = $newMoleculeSize;
          }

          #Get a list of read positions
          my @readPosToExtract = random_uniform_integer($readsToExtract, $startingPosition, $startingPosition+$moleculeSize-1);
          loop1: foreach my $gCoord (@readPosToExtract)
          {
            my $filePosToExtract;
            loop2: while($gCoord >= $startingPosition)
            {
                if(not defined $readPositionsInFile{$i}{$gCoord}) { --$gCoord; next loop2; }
                loop3: for(my $k = 0; $k < scalar(@{$readPositionsInFile{$i}{$gCoord}}); ++$k)
                {
                  if(${$readPositionsInFile{$i}{$gCoord}}[$k] != -1)
                  {
                    $filePosToExtract = ${$readPositionsInFile{$i}{$gCoord}}[$k];
                    ${$readPositionsInFile{$i}{$gCoord}}[$k] = -1;
                    last loop2;
                  }
                }
                --$gCoord;
            }
            next if not defined $filePosToExtract;

            #Extract reads and output
            #&Log("filePosToExtract: $filePosToExtract");
            if( seek($readsInputfh, $filePosToExtract, 0) != 1 ) { &LogAndDie("Seek failed on $filePosToExtract");}
            my $l1 = <$readsInputfh>; my $l2 = <$readsInputfh>; my $l3 = <$readsInputfh>; my $l4 = <$readsInputfh>;
            my $l5 = <$readsInputfh>; my $l6 = <$readsInputfh>; my $l7 = <$readsInputfh>; my $l8 = <$readsInputfh>;
            #print STDERR "$l1$l2$l3$l4$l5$l6$l7$l8";

            #Attach barcode
            my $read1OrRead2 = int(rand(2));
            if($read1OrRead2 == 0)
            {
              $l2 = "$selectedBarcode".substr($l2, 16);
              $l4 = "AAAFFFKKKKKKKKKK".substr($l4, 16);
            }
            elsif($read1OrRead2 == 1)
            {
              $l6 = "$selectedBarcode".substr($l6, 16);
              $l8 = "AAAFFFKKKKKKKKKK".substr($l8, 16);
            }
            else { die "[Attach barcode] Should never reach here."; }

            #Output reads
            print $fq1Outputfh "$l1$l2$l3$l4";
            print $fq2Outputfh "$l5$l6$l7$l8";
            --$readsCountDown;
            if($readsCountDown % 10000 == 0)
            { &Log("$readsCountDown reads remaining"); }
          }
        }
      }

      close $readsInputfh;
    }

    close $fq1Outputfh;
    close $fq2Outputfh;
    &Log("Simulate reads end");
  }
  #Simulate reads end

  0;
}

sub usage {
  die(qq/
    Usage:   $0 -r <reference> -p <output prefix> -b <barcodes> [options]

    Other options:
    -d <int>    Haplotypes to simulate [2]
    -e <float>  Per base error rate of the first read [0.010]
    -E <float>  Per base error rate of the second read [0.010]
    -i INT      Outer distance between the two ends for pairs [350]
    -s INT      Standard deviation of the distance for pairs [35]
    -x INT      Number of million reads in total to simulated [1000]
    -f INT      Mean molecule length in kbp [50]
    -t INT      n*1000 partitions to generate [1500]
    -m INT      Average # of molecules [10]

    /);
}

# Log routine
sub Log
{
  state $statusFH;
  if(not defined $statusFH)
  {
    open $statusFH, ">>$_[0]" or die "Error opening $_[0].\n";
    push @fileHandlersAtExit, $statusFH;
  }
  my $time = localtime;
  print $statusFH "$time: $_[0]\n";
  print STDERR "$time: $_[0]\n";
}

sub LogAndDie
{
  &Log(@_);
  die $!;
}
# Log routine end

sub LoadFaidx
{
  my $faidx = shift;
  my $boundary = shift;
  my $fn = shift;
  open my $fh, "$fn.fai" or &LogAndDie("Error opening faidx: $fn.fai");
  my $accumulation = 0;
  while(<$fh>)
  {
    chomp;
    my @a = split;
    $$faidx{acc}{"$a[0]"} = $accumulation;
    $$faidx{size}{"$a[0]"} = $a[1];
    push @$boundary, $accumulation;
    $accumulation += $a[1];
  }
  push @$boundary, $accumulation;
  close $fh;
  return $accumulation;
}

sub getChrSize { return ${$_[0]}{size}{$_[1]}; }

sub getChrStart { return ${$_[0]}{acc}{$_[1]}; }

sub GenomeCoord2Idx { return ${$_[0]}{acc}{$_[1]} + $_[2]; }

sub bSearch {
  my ( $elem, $list, $lowerLimit, $upperLimit ) = @_;
  my $max = $#$list;
  my $min = 0;

  my $index;
  while ( $max >= $min ) {
    $index = int( ( $max + $min ) / 2 );
    if    ( $list->[$index] < $elem ) { $min = $index + 1; }
    elsif ( $list->[$index] > $elem ) { $max = $index - 1; }
    else                              { last; }
  }
  if($elem >= $list->[$index]) { $$lowerLimit = $list->[$index]; $$upperLimit = $list->[$index+1]; }
  elsif($elem < $list->[$index]) { $$lowerLimit = $list->[$index-1]; $$upperLimit = $list->[$index]; }
  else {  die "bSearch: Should never reach here"; }
}

sub PoissonMoleculePerPartition
{
  state $mu = $_[0];
  state $i = 10000;
  state $pool;
  $i = 10000 if($mu != $_[0]);
  if($i == 10000)
  {
    @{$pool} = random_poisson(10000, $_[0]);
    $i = 0;
  }
  return ${$pool}[$i++];
}

sub PoissonMoleculeSize
{
  state $mu = $_[0];
  state $i = 10000;
  state $pool;
  $i = 10000 if($mu != $_[0]);
  if($i == 10000)
  {
    @{$pool} = random_poisson(10000, $_[0]);
    $i = 0;
  }
  return ${$pool}[$i++];
}

0;

