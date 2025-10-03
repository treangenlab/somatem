#!/usr/bin/env perl

# MANUAL FOR validate_nt_fasta.pl
=pod

=head1 NAME

validate_nt_fasta.pl -- makes sure this is a NT FASTA file

=head1 SYNOPSIS

 validate_nt_fasta.pl --fasta=/Path/to/infile.fasta [--cleaned_output=/Path/to/clean.fasta] [--max_seq_size=1000] [--max_N=10] [--peptide]
                     [--help] [--manual]

=head1 DESCRIPTION

 This script will check that this is a valid FASTA with no repeating sequence IDs
 and that all sequence characters are valid (i.e. A,T,G,C,W,S,M,K,R,Y,B,D,H,V,N)
 
=head1 OPTIONS

=over 3

=item B<-f, --fasta>=FILENAME

Input file in FASTA format. (Required) 

=item B<-co, --cleaned_output>=FILENAME

File to write a clean, valid, FASTA output to. (Default, dont write a clean output, just die if the input is bad)

=item B<-mx, --max_seq_size>=INT

Sequences greater than the value will be invalid (Default, no limit to sequence size)

=item B<-mn, --max_N>=INT

Seuqences with this many N's or more will fail (Default, no limit to the number of Ns)

=item <-p, --peptide>

This is a peptide FASTA file (Default, NT FASTA file)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.



=head1 AUTHOR

Written by Daniel Nasko, 
Center for Bioinformatics and Computational Biology, University of Maryland.

=head1 REPORTING BUGS

Report bugs to dnasko@umiacs.umd.edu

=head1 COPYRIGHT

Copyright 2017 The S2FAST Team.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS WITH NO DEFAULT
my($fasta,$cleaned_output,$peptide,$help,$manual);

## Arguments with defaults
my ($max_seq_size,$max_N) = (1000000000,1000000000); ## Set very high by default

GetOptions (	
				"f|fasta=s"	      => \$fasta,
                                "co|cleaned_output=s" => \$cleaned_output,
                                "mx|max_seq_size=s"   => \$max_seq_size,
                                "mn|max_N=s"            => \$max_N,
				"p|peptide"     =>      \$peptide,
                                "h|help"	=>	\$help,
				"m|manual"	=>	\$manual) || pod2usage({-exitval => 2, -verbose => 1, -output => \*STDERR});

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
$max_seq_size = int($max_seq_size); ## Convert to integer if user didn't do so
$max_N = int($max_N); ## Convert to integer if the user didn't do so
pod2usage( -msg  => "\n\n ERROR!  Required argument --fasta not found.\n\n", -exitval => 2, -verbose => 1)  if (! $fasta );
pod2usage( -msg  => "\n\n ERROR!  Required argument --max_seq_size must be an integer > 0 and <= 1 billion .\n\n", -exitval => 2, -verbose => 1)  if ( $max_seq_size <= 0  || $max_seq_size > 1000000000 );
pod2usage( -msg  => "\n\n ERROR!  Required argument --max_N must be an integer > 0 and <= 1 billion .\n\n", -exitval => 2, -verbose => 1)  if ( $max_N <= 0  || $max_N > 1000000000 );

### Global variables
my $header = "";
my $seq = "";
my ($no_seq,$dup_header,$bad_bases,$seqs_too_big,$too_many_N) = (0,0,0,0,0); ## counting the reasons sequences are tossed, only if we're spitting out a cleaned fasta file
my %Headers;
my ($fh,$ofh);

if ($cleaned_output) {
    open($ofh,'>',"$cleaned_output") || die "\n Cannot write to the cleaned_output file: $\n";
}
if ($fasta =~ m/\.gz$/) { ## if a gzip compressed infile
    open($fh, '-|', 'gunzip', '-c', $fasta) || die "\n\n Cannot open the input file: $fasta\n\n";
}
else { ## If not gzip comgressed
    open($fh,'<',"$fasta") || die "\n\n Cannot open the input file: $fasta\n\n";
}
while(<$fh>) {
    $_ =~ s/\r[\n]*/\n/gm;
    chomp;
    if ($_ =~ m/^>/) {
	unless ($header eq "") {
	    $seq = uc($seq);
	    evaluate_seq($header,$seq);
	}
	$seq = "";
	$header = $_;
    }
    else {
	$seq = $seq . $_;
    }
}
close($fh);
evaluate_seq($header,$seq);
if ($cleaned_output) { close($ofh); }

if ($cleaned_output && $no_seq > 0 || $cleaned_output && $dup_header > 0 || $cleaned_output && $bad_bases > 0 || $cleaned_output && $seqs_too_big > 0 || $cleaned_output && $too_many_N > 0) {
    print "\n WARNING: Cleaned output written to: $cleaned_output\n However, some sequences were thrown out for the following reasons:\n Duplicate headers = $dup_header\n Header with no sequence data = $no_seq\n Sequence contains invalid bases = $bad_bases\n Sequences were too big: $seqs_too_big\n Sequences with too many 'N's: $too_many_N\n\n";
}

exit 0;

sub evaluate_seq
{
    ## Input: A string we presume to be a nuclotide or Amino acid sequence
    ## Output: Nothing, die if this isn't a valid nucleotide or AA sequence in FASTA format
    my $h = $_[0];
    my $s = uc($_[1]);
    my $hold_s = uc($_[1]); ## Grab the 2nd element passed to the subroutine and uppercase it
    my $print=1;
    ## Seq length can't be zero
    if (length($hold_s) == 0) {
	if ($cleaned_output) {
	    $print = 0;
	}
	else {
	    die "\n FAIL: Your input FASTA file contains empty sequences (no base pairs): $header\n\n"; }
	$no_seq++;
    }
    if (length($hold_s) > $max_seq_size) { ## if a sequence is larger than the max size limit
	if ($cleaned_output) {
	    $print=0;
	}
	else {
	    die "\n FAIL: Your input FASTA file contains sequences larger than your cut-off of $max_seq_size : $header\n\n";}
	$seqs_too_big++;
    }
    if ($peptide) { ## if this is supposed to be a peptide sequence
	$s =~ tr/ARNDCQEGHILKMFPSTWYVX/                     /;
    }
    else {
	my $number_of_ns = $hold_s =~ tr/N/N/;
	if ($number_of_ns >= $max_N) {
	    if ($cleaned_output) {
		$print = 0;
	    }
	    else {
		die "\n FAIL: Your input FASTA file contains a sequence with too many 'N' bases.\n Limit = $max_N\n Actual = $number_of_ns\n Seq: $h\n\n";
	    }
	    $too_many_N++;
	}
	$s =~ tr/ATGCWSMKRYBDHVN/               /;
    }
    $s =~ s/ //g;
    if (length($s) > 0) {
	if ($cleaned_output) {
            $print = 0;
	}
	else {
	    die "\n FAIL: Your input FASTA file contains invalid bases (valid NT bases: A,T,G,C,W,S,M,K,R,Y,B,D,H,V,N; valid AA residues: A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X) the following sequence contains the characters '$s':\n\n$header\n\n";
	}
	$bad_bases++;
    }
    if (exists $Headers{$h}) {
	if ($cleaned_output) {
	    $print = 0;
	}
	else {
	    die "\n FAIL: Your input FASTA file does not contain unique identifiers. Offending sequence ID:\n\n$h\n\n";
	}
	$dup_header++;
    }
    $Headers{$h} = 1;
    if ($cleaned_output && $print == 1) {
	print $ofh $h . "\n" . $hold_s . "\n";
    }
}
