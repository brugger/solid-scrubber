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

my $dbase = 1;
my $dbh = DBI->connect('DBI:mysql:demo_01_s4', 'easih_ro') || die "Could not connect to database: $DBI::errstr";

my $chr_file = "chrX_hg18.fa";

my ($chr, $seq) = readfasta( $chr_file );

my %opts;
getopts('b:s:e:S:', \%opts);
usage() if ( $opts{h});

# my $start  = 55050070;
# my $end    = 55050080;
# $start     = 40391631;
# $end       = 40391651;

# $start     = 30146789;
# $end       = 30146809;

# $start     = $opts{ s } if ($opts{ s });
# $end       = $opts{ e } if ($opts{ e });


my $good_trans = legal_transitions();


my $bam_file = "../../A03.303.004.000/Octs_s4/MARIS/XLMR.Demo_01_s4_frag.bam";
$bam_file   = "/ifs/data/chris/project/XLMR/Demo01_S4_PE_no_saet.bam";
$bam_file   = "/ifs/data/chris/project/XLMR/Demo01_S4_PE_with_saet.bam";

$bam_file  = $opts{b} if ($opts{b});


my $snps = readin_snps( $opts{S}) if ( $opts{S});

print `samtools view -H $bam_file`;

my @entries;

while (<>) {
  chomp;
  my ($chr, $pos) = split(":", $_);

  analyse( $chr, $pos - 10, $pos + 10);
  @entries = undef;
}

sub analyse {
  my ( $chr, $start, $end) = @_;

  $chr =~ s/chr//;
  $chr = "chr$chr";
  
  my $region = "$chr:$start-$end";

  
  

  open (my $bam, "samtools view $bam_file $region | ") || die "Could not open 'stream': $!\n";

  my @res;

  while(<$bam>) {

    chomp;
    my @F = split("\t");
    my ($id, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, @opts) = @F;
    
    my $entry = { sam => \@F,
		  id       => $id,
		  flags    => $flags, 
		  chr      => $chr, 
		  pos      => $pos, 
		  end      => $pos + length($sequence) - 1};
    
    my $csfasta = $opts[0];
    
    next if ($cigar =~ /[IDNP]/);
    
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
      
      die "no csfasta for $id, $flags\n" if ( !$csfasta);
      $cigar = "1H$cigar";
    }
    else {
      # ensure that we have the csfasta entry, otherwise fetch it.. A bit clunky, but will do for now.
      if ( $csfasta !~ /CS:Z:/) {
	my $i = 1;
	while( $i < @opts ) {
	  $csfasta = $opts[$i++];
	  
	  last  if ($csfasta =~ /CS:Z:/);
	}
      }
      
      $csfasta =~ s/CS:Z://;
    }
    
    my $hlen = length($sequence);
    $csfasta =~ tr/0123/abcd/;
    
    # remove the primer base, and the colour that depends on it.
    $csfasta  =  substr($csfasta, 2);
    if ($flags & 0x0010 ) {
      $csfasta = reverse( $csfasta);
    }


    $csfasta  = patch_alignment($csfasta, $cigar);
    $csfasta  = substr($csfasta, 0, -1);

    my $gseq = substr($seq, $pos-1, $hlen);
    my $gcsfasta = substr(fasta2csfasta($gseq), 1);
    # odd hackity-hack thingy as the genome colour string is sometimes one colour to long(!!!????!!!)
    $gcsfasta = substr($gcsfasta, 0, -1) if ( length( $gcsfasta)> length($csfasta));
    
#    print "$csfasta\n$gcsfasta\n$cigar\n" if ($id eq $name);
#    print Dumper( $entry) if ( $id eq $name);
    

    my ($singles, $doubles) =  align($csfasta, $gcsfasta, $flags & 0x0010);
    

    $$entry{singles}  = $singles;
    $$entry{doubles}  = $doubles;
  

    push @entries, $entry;
    
  
    my @gcsf = split("", $gcsfasta);
    my @csf = split("", $csfasta);
    for(my $i = 0; $i<@csf;$i++) {
      next if ($pos - $start + $i <0);
      
      $res[$pos - $start + $i]{ref} = $gcsf[$i];
      $res[$pos - $start + $i]{pos} = $pos  + $i;
      $res[$pos - $start + $i]{$csf[$i]}++;
      
    }
  }
  
  my @col_bal;
  
#exit;
#die Dumper(\@res);
  
  foreach my $pos ( @res ) {
    
    next if ( ! $$pos{pos});
    last if ($$pos{pos} > $end );

#  print "$$pos{pos} || $$pos{pos} - 1 \n";
    
    next if ( $snps &&  $$snps{ chrX }{$$pos{pos}} || $$snps{ chrX}{$$pos{pos} + 1 });

    my ($total, $right) = (0,0);
    foreach my $colour ( 'a','b','c','d' ) {
      next if (!$$pos{$colour});
      $total += $$pos{$colour};
      $right += $$pos{$colour} if ( $$pos{ref} eq $colour);
    }
    
    next if ( ! $total);
    push @col_bal, [$$pos{ pos}, $right*100/$total];
  }
  
#print Dumper( \@col_bal );
  
  for(my $i=0; $i< @col_bal - 1; $i++ ) {
    
    if ( $col_bal[$i][1] < 60 && $col_bal[$i + 1][1] < 60) {
#    print "$i ($col_bal[$i][0]) --> $i +1 low coverage\n";
      filter($col_bal[$i][0], $col_bal[$i][0]+1);
    }
  }
  
}


# 
# 
# 
# Kim Brugger (20 Oct 2010)
sub filter {
  my ($start, $end) = @_;

  my ($dropped, $kept) = (0,0);

  foreach my $entry (@entries) {


    if (  $$entry{pos} - 1 <= $start && $$entry{end} >= $end && 
	  ( @{$$entry{singles}} || @{$$entry{doubles}})) {

      foreach my $single (@{$$entry{singles}}) {

	if (($single - 1 + $$entry{pos} <= $start &&
	     $single + 1 + $$entry{pos} >= $end)
	                || 
	    ($single - 1 + $$entry{pos} <= $start &&
	     $single + 1 + $$entry{pos} >= $start)
                        ||
	    ($single - 1 + $$entry{pos} <= $end &&
	     $single + 1 + $$entry{pos} >= $end)) {

	  $dropped++;
	  
	  $$entry{flags} -= 4;
	  goto POST_PRINT;
	}
      }  

      foreach my $double (@{$$entry{doubles}}) {

	if (($$double[0] + $$entry{pos} <= $start &&
	     $$double[1] + $$entry{pos} >= $end)
	                   || 
	    ($$double[0] + $$entry{pos} <= $start &&
	     $$double[1] + $$entry{pos} >= $start)
                           ||
	    ($$double[0] + $$entry{pos} <= $end &&
	     $$double[1] + $$entry{pos} >= $end)) {

	  $dropped++;
	  $$entry{flags} -= 4;
	  goto POST_PRINT;
	}
      }  
    }
    
    if ( $$entry{pos} == $end || $$entry{end}  == $end ) {
      

      goto POST_PRINT;
    }
    my $sam    =  $$entry{sam};
    @$sam[ 1 ] = $$entry{flags};
    
    print join("\t", @$sam) . "\n";
  POST_PRINT:
  }
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
  my ( $s1, $s2, $strand) = @_;

  return ([],[]) if ($s1 eq $s2);

  my @s1 = split("", $s1);
  my @s2 = split("", $s2);


  my (@singles, @align, @doubles);
  my $snp = 0;
  
  for(my $i = 0; $i < @s1; $i++) {
    if ( $s1[$i] eq $s2[$i] ) {
      $align[ $i ] = 1;
      $snp = 0;
    }
    else {
      $align[ $i] = 0;
      # Check and see if this is a double "error" IE a SNP
      if ( $i > 0 && ! $align[ $i - 1 ]) {

	if ( !$snp ) {
	  push @doubles, [$singles[-1], $i];
	  pop @singles;

	}
	$snp = 1;
	$doubles[-1][1] = $i;
      }
      else {
	push @singles, $i;
	$snp = 0;
      }
    }
  }

  my @long;
#  print Dumper( @doubles );

  foreach my $d ( @doubles ) {
    
    my $s1_d = join("",@s1[$$d[0]..$$d[1]]);
    my $s2_d = join("",@s2[$$d[0]..$$d[1]]);

#    die "$s1_d -- $s2_d ==>>> $$good_trans{$s1_d}{$s2_d} $$d[0]..$$d[1] $s1\n";

    if ( ! $$good_trans{$s1_d}{$s2_d}) {
      push @long, $d;      
    }
    
  }
  return( \@singles, \@long);
}




# 
# 
# 
# Kim Brugger (21 Oct 2010)
sub legal_transitions {
  my ( $length ) = @_;

  $length ||= 3;

  my %trans;

  foreach my $length (2, 3) {
    my @keys = makekeys($length);
    
    foreach my $right ( @keys ) {
      my $e = cs2fasta("C$right");
      
      foreach my $key ( @keys ) {
	my $last = cs2fasta("C$key");
	if ( $last eq $e) {
	  $trans{ $right }{ $key }++;
	  #      print "$input  C$key ==> $last eq $e\n";
	}
      }
    }
  }
  
  return \%trans;
}



# 
# 
# 
# Kim Brugger (20 Oct 2010)
sub cs2fasta {
  my ( $cs ) = @_;
  
  my %colourspace = (
    "Aa" => "A",
    "Ca" => "C",
    "Ga" => "G",
    "Ta" => "T",
    "Ab" => "C",
    "Cb" => "A",
    "Gc" => "A",
    "Ac" => "G",
    "Ad" => "T",
    "Td" => "A",
    "Cc" => "T",
    "Tc" => "C",
    "Cd" => "G",
    "Gd" => "C",
    "Gb" => "T",
    "Tb" => "G" );

  my @letters = split( //, $cs );
  my $first_base = $letters[0];
  for( my $i = 1; $i < @letters ; $i++ ){
    
    my $colour = $letters[$i];
    my $encoding = $first_base.$colour;
    $first_base = $colourspace{ $encoding };
    $letters[ $i ] = $first_base;    
  }
  shift @letters;
  return $letters[-1];
  return join("",@letters);

}


sub makekeys {
  my ($length) = @_;

  my @ALFA = split('', qq/abcd/);
  my $al = @ALFA;
  
  $length = shift || 2;
  my @s;
  my $nk = $al**$length;
  
  if ($length > 1) { 
    my @A = makekeys($length-1);
    for (my $i = 0, my $c = 0; $i<$nk; $i+=$al**($length-1), $c++) {
      for (my $j = 0; $j<$al**($length-1); $j++) {
        $s[$i+$j] = $ALFA[$c];  
        $s[$i+$j] .= $A[$j];
      }
    }
  }
  else {
    for (my $i = 0;$i< $al; $i++) {
      $s[$i] = "$ALFA[$i]";
    }
  }
  return @s;
}



# 
# 
# 
# Kim Brugger (13 Oct 2010)
sub fasta2csfasta {
  my ($fasta) = @_;


#   my %base2colour = (
#     'AA' => '0',
#     'AC' => '1',
#     'AG' => '2',
#     'AT' => '3',
#     'CA' => '1',
#     'CC' => '0',
#     'CG' => '3',
#     'CT' => '2',
#     'GA' => '2',
#     'GC' => '3',
#     'GG' => '0',
#     'GT' => '1',
#     'TA' => '3',
#     'TC' => '2',
#     'TG' => '1',
#     'TT' => '0');


  my %base2colour = (
    'AA' => 'a',
    'AC' => 'b',
    'AG' => 'c',
    'AT' => 'd',
    'CA' => 'b',
    'CC' => 'a',
    'CG' => 'd',
    'CT' => 'c',
    'GA' => 'c',
    'GC' => 'd',
    'GG' => 'a',
    'GT' => 'b',
    'TA' => 'd',
    'TC' => 'c',
    'TG' => 'b',
    'TT' => 'a');
  
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
      print STDERR "$cigar\n";
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
