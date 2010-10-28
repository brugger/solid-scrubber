#!/usr/bin/perl -w
# 
# 
# 
# 
# Kim Brugger (Oct 2010), contact: brugger@brugger.dk

use strict;
use Data::Dumper;
use warnings;


my $bam_file = shift;

my $chr_file = "chrX_hg18.fa";
$chr_file    = 'chrX.fa';
my ($chr, $seq) = readfasta( $chr_file );


my (@cs_splits, @ref_ids, @reads, @SNPs);

my $CS_PRE_POS     = 1;
my $CS_PRE_PRE_POS = 0;
my $MIN_SPLIT      = 20; # should be something like 10-15...

#this should be > readlength + the maximum number of expected inserts. 2xreadlength should do the trick
my $FILTER_BUFFER  = 100;

my $current_pos = undef;

my $good_trans = legal_transitions();

print `./samtools view -H $bam_file`;

open (my $bam, "./samtools view $bam_file | ") || die "Could not open 'stream': $!\n";
while(<$bam>) {

  chomp;
  my @F = split("\t");
  my ($id, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, @opts) = @F;
  
  my $entry = { sam    => \@F,
		id     => $id,
		flags  => $flags, 
		chr    => $chr, 
		pos    => $pos, 
		mapq   => $mapq,
		end    => $pos + length($sequence) - 1};
  
  my $csfasta = $opts[0];

  if ( $csfasta !~ /CS:Z:/) {
    my $i = 1;
    while( $i < @opts ) {
      $csfasta = $opts[$i++];
      
      last  if ($csfasta =~ /CS:Z:/);
    }
  }
      
  $csfasta =~ s/CS:Z://;
  $csfasta =~ tr/0123/abcd/;
  $csfasta  = substr($csfasta, 2);
  if ($flags & 0x0010 ) {
    $csfasta = reverse( $csfasta);
  }

  
  next if ($cigar =~ /[IDNP]/);


  if ( $current_pos && $current_pos != $pos )  {
    # the cs_splits array needs to be synced with this new pos. Bring forth the 
    # array as many places. If there is a gap, traverse the whole thing and reset
    # thewhole thing, so we start from fresh.
    
#    print STDERR "stepping from $current_pos =>>  $pos\n";
    for(my $i = 0; $i < $pos - $current_pos; $i++ ) {

      # This is a sliding array, keeping track of all the colour balances, and update the SNP array...
      shift @ref_ids if ( @ref_ids >= 2 );
      
      my ($total, $right) = (0,0);
      foreach my $colour ( 'a','b','c','d' ) {
	next if (!$cs_splits[$i]{$colour});
	$total += $cs_splits[$i]{$colour};
	$right += $cs_splits[$i]{$colour} if ( $cs_splits[$i]{ref} eq $colour);
      }

      if ( $total == 0 ) {
	@ref_ids = ();
	next;
      }

      push @ref_ids, [$cs_splits[$i]{pos}, $right*100/$total];

      if ( $ref_ids[ $CS_PRE_PRE_POS ] && $ref_ids[ $CS_PRE_POS ]  &&
	   $ref_ids[ $CS_PRE_PRE_POS ][ 1 ] < $MIN_SPLIT && 
	   $ref_ids[ $CS_PRE_POS ][ 1 ] < $MIN_SPLIT ) {
	  
#	  push @SNPs, [$ref_ids[ $CS_PRE_PRE_POS ][0], $ref_ids[ $CS_PRE_POS ][0]];
	  push @SNPs, [$ref_ids[ $CS_PRE_PRE_POS ], $ref_ids[ $CS_PRE_POS ]];
	  
      }

      if ( @reads &&  $reads[0]{ pos } + $FILTER_BUFFER <= $pos + $i ) {
	
	
	while ( 1 ) {	   
	  last 	if ( !@SNPs  || $reads[0]{ pos } + $FILTER_BUFFER > $SNPs[0][1][0]  );
	  my @snp = shift @SNPs;
#	  print "Removing SNP :: $snp[0][0][0] - $snp[0][1][0] \n";
	}

#	print "Removing reads mapped before: ". ($pos + $i - $FILTER_BUFFER ). " checking them against: ". @SNPs ." SNPs (read[0]:: $reads[0]{ pos } --> $reads[0]{ end })\n";
#	print "Pre-removal nr of reads: " . @reads ."\n";
	
	my ( $s, $d, $e) = (0,0,0);

	while ( @reads > 0 && $reads[0]{ pos } + $FILTER_BUFFER < $pos + $i ) {
	  my $read = shift @reads;
	  
#	  print "$$read{pos} ==> $$read{ end } \n";

	  goto PRINT if ( ! $$read{singles} && ! $$read{doubles} );
	  
#	  print "Checking againt snps\n";
	  foreach my $snp ( @SNPs ) {
	    my ($snp_start, $snp_end) = ($$snp[0][0], $$snp[1][0]);
	    
	    foreach my $single (@{$$read{singles}}) {

	      if (($single - 1 + $$read{pos} <= $snp_start &&
		   $single + 1 + $$read{pos} >= $snp_end)
		  || 
		  ($single - 1 + $$read{pos} <= $snp_start &&
		   $single + 1 + $$read{pos} >= $snp_start)
		  ||
		  ($single - 1 + $$read{pos} <= $snp_end &&
		   $single + 1 + $$read{pos} >= $snp_end)) {
		
		$$read{flags} -= 4;
		$$read{mapq}   = 2;	
		$s++;
#		print "One single\n";
		
		goto PRINT;
	      }
	    }  

	    foreach my $double (@{$$entry{doubles}}) {
	      
	      if (($$double[0] + $$read{pos} <= $snp_start &&
		   $$double[1] + $$read{pos} >= $snp_end)
		  || 
		  ($$double[0] + $$read{pos} <= $snp_start &&
		   $$double[1] + $$read{pos} >= $snp_start)
		  ||
		  ($$double[0] + $$read{pos} <= $snp_end &&
		   $$double[1] + $$read{pos} >= $snp_end)) {
		
		$$read{flags} -= 4;
		$$read{mapq}   = 2;		
		$d++;
#		print "One double\n";
		
		goto PRINT;
	      }
	    }  
	  
	  
	    if ( $$read{pos} == $snp_end || $$read{end}  == $snp_end ) {
	    
	      $$read{flags} -= 4;
	      $$read{mapq}   = 2;		
	      $e++;
#	      print "One ends\n";
	      goto PRINT;
	    }
	  }
	PRINT:
	  my $sam    = $$read{ sam   };
	  @$sam[ 1 ] = $$read{ flags };
	  @$sam[ 4 ] = $$read{ mapq  };
	  
	  print join("\t", @$sam) . "\n";

	  
	}
	  
	print STDERR "$s, $d, $e\n";
#	print "Post-removal nr of reads: " . @reads ."\n";

	

      }

    }
  }
  $current_pos = $pos;
  
#  print "pos: $current_pos, buffered reads: " . @reads . "\n";

#  print Dumper( \@SNPs) if ( @SNPs );

  my $hlen = length($sequence);
  my $gseq = substr($seq, $pos - 1, $hlen + 1);
  my $gcsfasta = substr(fasta2csfasta($gseq), 1);
  # for some odd reason, I cannot be arsed to figure out right now...
  $gcsfasta = substr($gcsfasta, 0, length($csfasta)) if ( length( $gcsfasta) > length($csfasta));

  $csfasta  = patch_alignment($csfasta, $cigar);

  my ($singles, $doubles) =  align($csfasta, $gcsfasta, $flags & 0x0010);
   
  $$entry{singles}  = $singles;
  $$entry{doubles}  = $doubles;

  push @reads, $entry;
#  print "$csfasta\n";

  my @gcsf = split("", $gcsfasta);
  my @csf = split("", $csfasta);
  for(my $i = 0; $i<@csf;$i++) {
    $cs_splits[$i]{ref} = $gcsf[$i];
    $cs_splits[$i]{pos} = $pos  + $i;
    $cs_splits[$i]{$csf[$i]}++;
  }
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
      $align[ $i ] = 0;
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

#  print "$s1\n".join("", @align)."\n$s2\n\n";

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

#  return $seq;

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
