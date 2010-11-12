#!/usr/bin/perl -w
# 
# SOLiD-scrubber: Removes ambiguous mappings around SNPs for cleaner SNP calling.
# 
# 
# Kim Brugger (Oct 2010), contact: brugger@brugger.dk

use strict;
use Data::Dumper;
use warnings;
use Getopt::Std;

my %opts;
getopts('b:d:R:B:hs:r:UM:m:', \%opts);
my $bam_file = $opts{'b'} || Usage();
my $chr_file = $opts{'R'} || Usage();

die "'$chr_file' does not exist\n" 
    if ( ! -e $chr_file );
die "index does not exist for '$chr_file', please create one with samtools faidx\n"
    if (! -e "$chr_file.fai");

my $samtools    = find_program('samtools');

my $gcsfasta;
my $CS_PRE_POS     = 1;
my $CS_PRE_PRE_POS = 0;
my $MIN_SPLIT      = $opts{'s'} || 60; # should probably be something like 10-15...
my $MIN_DEPTH      = $opts{'d'} || 15; # I guess something between 7 and 20 should be goo
#this should be > readlength + the maximum number of expected inserts. 2xreadlength should do the trick
my $FILTER_BUFFER  = $opts{'B'} || 100;
my $set_unmapped   = $opts{'U'};
my $set_mapq_score = $opts{'M'};
my $min_mapq_score = $opts{'m'} || 0;

if ( ! $set_unmapped && ! defined $set_mapq_score) {
  print STDERR  "Please use either the -M or the -U flag otherwise the reads will not be changed!!!\n";
  Usage();
}


my $current_pos = undef;
my $good_trans = legal_transitions();

print `$samtools view -H $bam_file`;
my ( $s, $d, $e, $b_pre, $b_post) = (0,0,0, 0, 0);

my $region = $opts{'r'};
#my $region = "chrX:2,789,549-2,789,590";
#$region = "chrX:2,789,818-2,789,859";
#$region = "chrX:2,789,818-3,789,818";
#$region = "chrX";

if ( ! $region ) {
  # loop through the regions one by one, but only keep the most current chr in memory.
  
  open(my $spipe, "$samtools view -H $chr_file | ") || die "Could not open '$chr_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      if ( $field =~ /SN:(.*)/) {
	$gcsfasta = readfasta( $chr_file, $1 );
	analyse($1)
      }
    }
  }
}
else {
  $region =~ s/,//g;
  if ($region =~ /^(\w+):\d+-\d+/ || $region =~ /^(\w+):\d+\z/ ) {
    $gcsfasta = readfasta( $chr_file, $1 );
    analyse($region)
  }
  else {
    $gcsfasta = readfasta( $chr_file, $region );
    analyse($region)
  }
}


# 
# 
# 
# Kim Brugger (05 Nov 2010)
sub analyse {
  my ($region) = @_;

  my (@cs_splits, @ref_ids, @reads, @SNPs);

  open (my $bam, "$samtools view $bam_file $region | ") || die "Could not open 'stream': $!\n";
  while(<$bam>) {

    chomp;
    my @F = split("\t");
    my ($id, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, @opts) = @F;
    
    my $entry = { sam    => \@F,
		  id     => $id,
		  flags  => $flags, 
		  chr    => $chr, 
		  cigar  => $cigar,
		  pos    => $pos, 
		  mapq   => $mapq,
		  end    => $pos + length($sequence) - 1};
    

    if ($flags & 0x0004 ) {
      push @reads, $entry;
      next;
    }

    if ( $mapq <= $min_mapq_score ) {
      push @reads, $entry;
      next;
    }


    my $csfasta = $opts[0];
    
    if ( $csfasta !~ /CS:Z:/) {
      my $i = 1;
      while( $i < @opts ) {
	$csfasta = $opts[$i++];
	
	last  if ($csfasta =~ /CS:Z:/);
      }
    }
    
    $csfasta =~ s/CS:Z://;

    my $BWA = 1;
    
    if ( $BWA ) { 
      $csfasta  = substr($csfasta, 3);
      $csfasta =~ tr/0/O/;
    }
    else {
      $csfasta  = substr($csfasta, 2);
      $csfasta =~ tr/0/O/;
      if ($flags & 0x0010 ) {
	$csfasta = reverse( $csfasta);
      }
    }
      
    if ( $current_pos && $current_pos != $pos )  {
      
      
      
      # the cs_splits array needs to be synced with this new pos. Bring forth the 
      # array as many places. If there is a gap, traverse the whole thing and reset
      # thewhole thing, so we start from fresh.
      
#    print STDERR "stepping from $current_pos =>>  $pos\n";
      for(my $i = 0; $i < $pos - $current_pos; $i++ ) {
	
	# This is a sliding array, keeping track of all the colour balances, and update the SNP array...
	shift @ref_ids if ( @ref_ids >= 2 );
	
	next if (! $cs_splits[0]{total} || $cs_splits[0]{total} == 0);
	
	
	my ($right) = (0);
	map { $right += $cs_splits[0]{ $_ } if ( $cs_splits[0]{ $_ } && $cs_splits[0]{ref} eq $_)} ( 'O','1','2','3');
	
	push @ref_ids, [$cs_splits[0]{pos}, $right*100/$cs_splits[0]{total}, $cs_splits[0]{total}];
	shift @cs_splits;
	
	if ( @ref_ids == 2 && 
	     $ref_ids[ $CS_PRE_PRE_POS ] && $ref_ids[ $CS_PRE_POS ]  &&
	     $ref_ids[ $CS_PRE_PRE_POS ][ 1 ] < $MIN_SPLIT && 
	     $ref_ids[ $CS_PRE_POS ][ 1 ]     < $MIN_SPLIT && 
	     $ref_ids[ $CS_PRE_POS ][ 2 ]     > $MIN_DEPTH     &&
	     $ref_ids[ $CS_PRE_PRE_POS ][ 2 ] > $MIN_DEPTH ) {


	  print STDERR "SNP at pos $region:".($current_pos  + $i)." (depth: $ref_ids[ $CS_PRE_POS ][ 2 ]) splits: ($ref_ids[ $CS_PRE_PRE_POS ][1], $ref_ids[ $CS_PRE_POS ][1])\n";
	  push @SNPs, [$ref_ids[ $CS_PRE_PRE_POS ], $ref_ids[ $CS_PRE_POS ]];
	}
	
	if ( @reads &&  $reads[0]{ pos } + $FILTER_BUFFER <= $pos + $i ) {
	  #remove SNPs further behind that we will ever be seeing again...
#	print STDERR "Buffered SNPs: " . @SNPs . "\n". " ($reads[0]{ pos } + $FILTER_BUFFER > $SNPs[0][1][0]); \n";
	  while ( @SNPs ) { 
#	  print STDERR Dumper( $SNPs[0]);
	    
	    last 	if ($reads[0]{ pos } - $FILTER_BUFFER < $SNPs[0][1][0]);
	    my @snp = shift @SNPs;
	  }
#	print STDERR "Buffered SNPs: " . @SNPs . "\n". " ($reads[0]{ pos } + $FILTER_BUFFER > $SNPs[0][1][0]); \n";
#	print STDERR "---------------------------------\n";
#	  exit;
	  
	  #scrub all the reads that are far enough away..
	  while ( @reads > 0 && $reads[0]{ pos } + $FILTER_BUFFER < $pos + $i ) {
	    my $read = shift @reads;
	    
	    scrub( $read, \@SNPs);
	    print_sam( $read );
	  }
	  
#	print "Post-removal nr of reads: " . @reads ."\n";
	  
	}
	
      }
    }
    $current_pos = $pos;
    
    
#  print "pos: $current_pos, buffered reads: " . @reads . "\n";
    
    
    my $hlen = length($csfasta);
    if ( 1 &&  $cigar =~ /[IDNP]/) {
      my $t_cigar = $cigar;
      my $padding = 2;
      $t_cigar =~ s/(\d+)[ID]/ {$padding += $1}/ge;
      $hlen += $padding;
    }
    
    my $gseq = substr($gcsfasta, $pos - 1, $hlen);
    # for some odd reason, I cannot be arsed to figure out right now...
#  $gseq = substr($gseq, 0, length($csfasta)) if ( length( $gseq) > length($csfasta));
    
    ($csfasta, $gseq)  = patch_alignment($csfasta, $gseq, $cigar);
    
    
    my ($singles, $doubles, $a) =  align($csfasta, $gseq, $flags & 0x0010);
    
    $$entry{singles}  = $singles;
    $$entry{doubles}  = $doubles;
    $$entry{a}        = $a;
    
    push @reads, $entry;
#  print "$csfasta\n";
    
    my @gcsf = split("", $gseq);
    my @csf = split("", $csfasta);
    my $gaps = 0;
    for(my $i = 0; $i<@csf;$i++) {
      
      # there is an insert in the reference, so this number needs to 
      # be subtracted to get the real genome position.
      if ($gcsf[ $i ] eq "-") {
	$gaps++;
	next;
      }

      $cs_splits[$i - $gaps]{ ref      } = $gcsf[ $i ];
      $cs_splits[$i - $gaps]{ pos      } = $pos  + $i;
      $cs_splits[$i - $gaps]{ $csf[$i] }++;
      next if ($csf[ $i ] eq "-");
      $cs_splits[$i - $gaps]{ total    }++;

      
    }
  }


# empty the reads buffer..
  while ( @reads ) {
    my $read = shift @reads;
    
    scrub( $read, \@SNPs);
    print_sam( $read );
    
  }
}




# 
# 
# 
# Kim Brugger (28 Oct 2010)
sub print_sam {
  my ( $read ) = @_;

#  return;
  my $sam    = $$read{ sam   };
  @$sam[ 1 ] = $$read{ flags };
  @$sam[ 4 ] = $$read{ mapq  };
  
  print join("\t", @$sam) . "\n";
}


# 
# 
# 
# Kim Brugger (28 Oct 2010)
sub scrub {
  my ( $read, $SNPs ) = @_;

  return if ( $$read{indel} );

  return if ( ! $$read{singles} && ! $$read{doubles} );
	  
  foreach my $snp ( @$SNPs ) {
    my ($snp_start, $snp_end) = ($$snp[0][0], $$snp[1][0]);
#    next if ( ! $snp_start || ! $snp_end);

    if ( $$read{ pos } - 5 > $snp_end) {
#      print STDERR "Pre-Bailing... $$read{pos} -> $$read{end} vs $snp_start => $snp_end\n";
      $b_post++;
      next;
    }


    if ( $$read{ end } + 5 < $snp_start) {
#      print STDERR "Bailing...\n";
      $b_pre++;
      return;
    }

#    next;
      
    foreach my $single (@{$$read{singles}}) {

      if (($single - 1 + $$read{pos} <= $snp_start &&
	   $single + 1 + $$read{pos} >= $snp_end)
	  || 
	  ($single - 1 + $$read{pos} <= $snp_start &&
	   $single + 1 + $$read{pos} >= $snp_start)
	  ||
	  ($single - 1 + $$read{pos} <= $snp_end &&
	   $single + 1 + $$read{pos} >= $snp_end)) {

#	print STDERR "$$read{id} -- $$read{cigar}\n$$read{a}";
	
	
	$$read{flags} -= 4 if ( $set_unmapped );
	$$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);
	$s++;
	return;
      }
    }  

#    next;

    foreach my $double (@{$$read{doubles}}) {
      
      if (($$double[0] + $$read{pos} <= $snp_start &&
	   $$double[1] + $$read{pos} >= $snp_end)
	  || 
	  ($$double[0] + $$read{pos} <= $snp_start &&
	   $$double[1] + $$read{pos} >= $snp_start)
	  ||
	  ($$double[0] + $$read{pos} <= $snp_end &&
	   $$double[1] + $$read{pos} >= $snp_end)) {
	
#	print STDERR "$$read{id} -- $$read{cigar}\n$$read{a}";
#	print $$read{a};

	$$read{flags} -= 4 if ( $set_unmapped );
	$$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);
#	$$read{flags} -= 4;
#	$$read{mapq}   = 2;		
	$d++;
	return;
      }
    }  	  	  
	  
    if ( $$read{pos} == $snp_end || $$read{end}  == $snp_end ) {
	    
#      print STDERR "$$read{id} -- $$read{cigar}\n$$read{a}";
#      print $$read{a};
      $$read{flags} -= 4 if ( $set_unmapped );
      $$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);

#      $$read{flags} -= 4;
#      $$read{mapq}   = 2;		
      $e++;
      return;
    }
  }

  return;

}



# 
# probe vs genome
# 
# Kim Brugger (13 Oct 2010)
sub align {
  my ( $s1, $s2, $strand) = @_;

  return ([],[], "NA\n") if ($s1 eq $s2);

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

  my @long;

  foreach my $d ( @doubles ) {
    
    my $s1_d = join("",@s1[$$d[0]..$$d[1]]);
    my $s2_d = join("",@s2[$$d[0]..$$d[1]]);

    if ( ! $$good_trans{$s1_d}{$s2_d}) {
      push @long, $d;      
    }
    
  }

  if ( 1 ) {
    my $align = join("", @align);
    $align =~ tr/01/ \|/;
#    print STDERR "$s1\n$align\n$s2\n\n" if ( @singles || @long);
    return( \@singles, \@long, "$s1\n$align\n$s2\n\n");
  }

  return( \@singles, \@long, );
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
    "AO" => "A",
    "CO" => "C",
    "GO" => "G",
    "TO" => "T",
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

  my @ALFA = split('', qq/O123/);
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

   my %base2colour = (
     'AA' => 'O',
     'AC' => '1',
     'AG' => '2',
     'AT' => '3',
     'CA' => '1',
     'CC' => 'O',
     'CG' => '3',
     'CT' => '2',
     'GA' => '2',
     'GC' => '3',
     'GG' => 'O',
     'GT' => '1',
     'TA' => '3',
     'TC' => '2',
     'TG' => '1',
     'TT' => 'O');

  my $res = "";
#  for ( my $i = 0;$i< 3000000 ; $i++ ) {
  for ( my $i = 0;$i< length($fasta) - 1 ; $i++ ) {
    $res .= $base2colour{ substr( $fasta, $i, 2)} || ".";
  }

  return $res;
}

sub patch_alignment {
  my ( $read, $ref, $cigar ) = @_;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)

  if ( $cigar !~ /[HDIS]/) {
    $ref = substr($ref, 0, length($read));
    return ($read, $ref);
  }
  
  my @read  = split("", $read );
  my @ref   = split("", $ref );

  my $ref_cigar = $cigar;
  $ref_cigar =~ s/^\d+[HDS]//;

  my (@cigar) = $ref_cigar =~ /(\d+\w)/g;

  my $offset = 0;
  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d+)(\w)/;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "I") {
      my @dashes = split("", "-"x$length);
      splice(@ref, $offset, 0, @dashes);
    }
    elsif ( $type eq "S" || $type eq "H") {
      splice(@ref,  $offset, $length);
    }    
  }

  if (1){
    $offset = 0;
    my (@cigar) = $cigar =~ /(\d+\w)/g;
    foreach my $patch ( @cigar ) {
      my ($length, $type) =  $patch =~ /^(\d+)(\w)/;

      if ( $type eq 'M') {
	$offset += $length;
	next;
      }
      elsif ( $type eq "D") {
	my @dashes = split("", "-"x$length);
	splice(@read,  $offset, 0, @dashes);
	$offset += $length;
      }
      elsif ( $type eq "H" || $type eq "S") {
      splice(@read,  $offset, $length);
      }    
    
    }
  }
  $read = join("", @read);
  $ref  = join("", @ref );

  $ref = substr($ref, 0, length($read));
  
  return ($read, $ref);
}



# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub patch_alignment_old {
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
    "AO" => "A",
    "CO" => "C",
    "GO" => "G",
    "TO" => "T",
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
  my ($file, $region) = @_;  

  my $sequence;
  my $header;
  
  open (my $f, "$samtools faidx $file $region |" ) || die "Could not open $file:$1\n";
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
  

  return (fasta2csfasta($sequence));
}




# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program) = @_;
  
  my $username = scalar getpwuid $<;
  
  my @paths = ("/home/$username/bin/",
	       "./",
	       "/usr/local/bin");
  
  foreach my $path ( @paths ) {
    return "$path/$program" if ( -e "$path/$program" );
  }

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );
  
  die "Could not find '$program'\n";
}


# 
# 
# 
# Kim Brugger (05 Nov 2010)
sub Usage {
  $0 =~ s/.*\///;
  die "USAGE: $0 -b<am file> -R<eference genome (fasta)> -d[ min depth, default=15] -s[ min Split, default=60] -B[uffer, default=100] -M[ set mapq score for offending reads] -U[n set mapped flag for offending reads]\n";

  # tests :::: odd SNP reporting:  10:74879852-74879852
  # large indel: 10:111800742
}
