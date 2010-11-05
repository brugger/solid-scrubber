#!/usr/bin/perl 
# 
# Digging into the true quality of solid sequencing
# 
# 
# Kim Brugger (13 Oct 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

use DBI;

my $dbase = 0;
my $dbh = DBI->connect('DBI:mysql:demo_01_s4', 'easih_ro') || die "Could not connect to database: $DBI::errstr";


my $chr_file = "chrX_hg18.fa";

my ($chr, $seq) = readfasta( $chr_file );

my %opts;
getopts('b:s:e:S:', \%opts);
usage() if ( $opts{h});



my $start  = 55050070;
my $end    = 55050080;
$start     = 40391631;
$end       = 40391651;

$start     = 71410411;
$end       = 72410411;
#$end       = 71410431;

$start     = $opts{ s } if ($opts{ s });
$end       = $opts{ e } if ($opts{ e });

my $region = "chrX:$start-$end";

my $bam_file = "../../A03.303.004.000/Octs_s4/MARIS/XLMR.Demo_01_s4_frag.bam";
$bam_file   = "/ifs/data/chris/project/XLMR/Demo01_S4_PE_no_saet.bam";
$bam_file   = "/ifs/data/chris/project/XLMR/Demo01_S4_PE_with_saet.bam";

$bam_file   = "/ifs/data/chris/project/XLMR/Demo01_S4_F3_no_saet.bam";

$bam_file  = $opts{b} if ($opts{b});


my $snps = readin_snps( $opts{S}) if ( $opts{S});

open (my $bam, "samtools view $bam_file $region | ") || die "Could not open 'stream': $!\n";

my $full = 1;
my @res;
my @di_res;


while(<$bam>) {

  chomp;
  my @F = split("\t");
  my ($id, $flags, $chr, $pos, $cigar, $sequence, $csfasta) = ($F[0], $F[1], $F[2], $F[3], $F[5], $F[9], $F[13]); 

  next if ($cigar =~ /[IDNP]/ || $flags & 0x0004);

  if ($flags & 0x0040 ) {
    $id .= "_F3";
  }
  else {
    $id .= "_F5-P2";
  }

  if ( $dbase ) { 
    
    $id =~ s/XLMR.Demo_01_s4:1://;;
    $id =~ s/:/_/g;

    my $sth = $dbh->prepare("SELECT csfasta, qual FROM data WHERE id='$id'");
    $sth->execute();
    ($csfasta, undef)  = $sth->fetchrow_array;
    
#    next if ( !$csfasta);
    die "no csfasta for $id, $flags\n" if ( !$csfasta);

    
    $cigar = "1H$cigar";

#    print "$id $csfasta \n";
#    exit;
  }
  else {
    if ( $csfasta !~ /CS:Z:/) {
      my $i = 14;
      while( $i < @F ) {
	$csfasta = $F[$i++];
	last  if ($csfasta =~ /CS:Z:/);
      }
    }
    next if ( $csfasta !~ /CS:Z:/);    

    $csfasta =~ s/CS:Z://;
    
  }

#   my $fasta = csfasta2fasta($csfasta);
#   if ($flags & 0x0010 ) {
#     $fasta =~ tr/ACGT/TGCA/;
#     $fasta = reverse($fasta);
#   }
#   $fasta    = patch_alignment($fasta, $cigar);

  my $hlen = length($sequence);

  # remove the primer base, and the colour that depends on it.
  $csfasta  = substr($csfasta, 2);
  if ($flags & 0x0010 ) {
    $csfasta = reverse( $csfasta);
  }

  if ( ! $full ) {
    $csfasta  = patch_alignment($csfasta, $cigar);
    $csfasta  = substr($csfasta, 0, -1);
  }
  else {

    if ( $cigar =~ /^(\d+)H/) {
      $pos -= $1;
    }
    $hlen = length( $csfasta ) + 1;
  }

  my $gseq = substr($seq, $pos-1, $hlen);
  my $gcsfasta = substr(fasta2csfasta($gseq), 1);
  # odd hackity-hack thingy as the genome colour string is sometimes one colour to long(!!!????!!!)
  $gcsfasta = substr($gcsfasta, 0, -1) if ( length( $gcsfasta)> length($csfasta));
  my $identical = $csfasta eq $gcsfasta || 0;
  
  if ( 0  &&  !$identical ) {
    print "$id $pos $cigar (seq, genome)\n";
    print align($csfasta, $gcsfasta);
  }

  
#  print "$id, $chr, $pos, $cigar , $sequence, $csfasta, $fasta --> $gseq [$gcsfasta]\n";
#  print "$id, $chr, $pos, $csfasta, $gcsfasta ( $identical ) $hlen $cigar\n";

  my @gcsf = split("", $gcsfasta);
  my @csf = split("", $csfasta);
  for(my $i = 0; $i<@csf;$i++) {
    next if ($pos - $start + $i <0);

    $res[$pos - $start + $i]{ref} = $gcsf[$i];
    $res[$pos - $start + $i]{pos} = $pos  + $i;
    $res[$pos - $start + $i]{$csf[$i]}++;

    if ( $i + 1 < @csf ) {

      $di_res[$pos - $start + $i]{ref} = $gcsf[$i].$gcsf[$i+1];
      $di_res[$pos - $start + $i]{pos} = $pos  + $i;
      $di_res[$pos - $start + $i]{ $csf[$i].$csf[$i+1] }++;
    }

  }
}


foreach my $pos ( @res ) {

  next if ( ! $$pos{pos});
  last if ($$pos{pos} > $end );

#  print "$$pos{pos} || $$pos{pos} - 1 \n";

  next if ( $snps &&  $$snps{ chrX }{$$pos{pos}} || $$snps{ chrX}{$$pos{pos} + 1 });

  my $s = "chrX:$$pos{pos}\t$$pos{ref}\t";

  my ($total, $right) = (0,0);
  foreach my $colour ( 0,1,2,3 ) {
    $$pos{$colour} = "0" if ( !$$pos{$colour});
    $s .= "$$pos{$colour}\t";
    $total += $$pos{$colour};
    $right += $$pos{$colour} if ( $$pos{ref} == $colour);
    
  }
  next if ( ! $total);
  $s .= sprintf("$total\t%.2f", $right*100/$total);
  print "$s\n";
}


#print Dumper($$snps{ chrX });


exit;

print "           \tref\t"; 
foreach my $c1 ( 0,1,2,3 ) {
  foreach my $c2 ( 0,1,2,3 ) {
    print "$c1$c2\t";
  }
}
print "\n";

foreach my $pos ( @di_res ) {

  next if ( ! $$pos{pos});
  last if ($$pos{pos} > $end );
  print "chrX:$$pos{pos}\t$$pos{ref}\t";

  my ($total, $right) = (0,0);
  foreach my $c1 ( 0,1,2,3 ) {
    foreach my $c2 ( 0,1,2,3 ) {
      my $di_col = "$c1$c2";
      $$pos{ $di_col } = "0" if ( !$$pos{ $di_col });
      print "$$pos{ $di_col }\t";
      $total += $$pos{$di_col};
      $right += $$pos{$di_col} if ( $$pos{ref} == $di_col );
    }
  }
  printf("$total\t%.2f", $right*100/$total);
  print "\n";
}





# 
# 
# 
# Kim Brugger (19 Oct 2010)
sub readin_snps {
  my ( $file ) = @_;

  my %res;

  open(my $f, $file) || die "Could not '$file': $!\n";
  while(<$f>) {
    next if (/\#/);
    chomp;
    my ( $chr, $pos ) = (split("\t"))[0, 1];
    $res{ $chr}{$pos }++;
  }
  close( $f );

  return \%res;
}


# 
# 
# 
# Kim Brugger (13 Oct 2010)
sub align {
  my ( $s1, $s2) = @_;
  my @s1 = split("", $s1);
  my @s2 = split("", $s2);
  my @align;
  my $errors = 0;
  for(my $i = 0; $i < @s1; $i++) {
    if ( defined $s1[$i] && defined $s2[$i] && $s1[$i] eq $s2[$i]) {
      $align[ $i] = "|";
    }
    else {
      $align[ $i] = " ";
      $errors++;
    }
  }
  push @align, sprintf("\terrors; %.2f%%", $errors*100/int(@s1));
  return join("", @s1, "\n", @align, "\n", @s2, "\n")
  
}



# 
# 
# 
# Kim Brugger (13 Oct 2010)
sub fasta2csfasta {
  my ($fasta) = @_;


  my %base2colour = (
    'AA' => '0',
    'AC' => '1',
    'AG' => '2',
    'AT' => '3',
    'CA' => '1',
    'CC' => '0',
    'CG' => '3',
    'CT' => '2',
    'GA' => '2',
    'GC' => '3',
    'GG' => '0',
    'GT' => '1',
    'TA' => '3',
    'TC' => '2',
    'TG' => '1',
    'TT' => '0');
  
  my @letters = split( //, $fasta );
  my $first_base = $letters[0];
  for( my $i = 1; $i < @letters ; $i++ ){
    
    my $encoding = $base2colour{ $first_base.$letters[$i] };
    $first_base = $letters[ $i ];
    $letters[ $i  ] = $encoding;
    
  }

  return ( join("",@letters));;
}





# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub patch_alignment {
  my ( $seq, $cigar ) = @_;

#  return ($seq) if ( $cigar !~ /[DIS]/);
  
  my @seq  = split("", $seq );


  my (@cigar) = $cigar =~ /(\d*\w)/g;

  my $offset = 0;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)


  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d*)(\w)/;
    $length ||= 1;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "D") {
      my @dashes = split("", "-"x$length);
      splice(@seq,  $offset, 0, @dashes);
      $offset += $length;
    }
    elsif ( $type eq "I" || $type eq "S" || $type eq "H") {
      splice(@seq,  $offset, $length);
    }    

  }


  return (join("", @seq));
}



# 
# 
# 
# Kim Brugger (13 Oct 2010)
sub csfasta2fasta {
  my ($csfasta) = @_;

  my %colourspace = (
    "A0" => "A",
    "C0" => "C",
    "G0" => "G",
    "T0" => "T",
    "A1" => "C",
    "C1" => "A",
    "G2" => "A",
    "A2" => "G",
    "A3" => "T",
    "T3" => "A",
    "C2" => "T",
    "T2" => "C",
    "C3" => "G",
    "G3" => "C",
    "G1" => "T",
    "T1" => "G" );

  my @letters = split( //, $csfasta );
  my $first_base = $letters[0];
  for( my $i = 1; $i < @letters ; $i++ ){
    
    my $colour = $letters[$i];
    my $encoding = $first_base.$colour;
    $first_base = $colourspace{ $encoding };
    $letters[ $i ] = $first_base;    
  }
  shift @letters;
  return(join("",@letters));
}



#
# Read the fasta files and puts entries into a nice array
#
sub readfasta {
  my ($file) = @_;

  my $sequence;
  my $header;
  open (my $f, $file) || die "Could not open $file:$1\n";
  while (<$f>) {
    chomp;
    if (/^\>/) {
      if ($header) { # we have a name and a seq
	return ($header, $sequence);
      }
      $header = $_;
      $header =~ s/^\>//;
    }
    else {$sequence .= $_;}
  }
  

  return ($header, $sequence);
}
