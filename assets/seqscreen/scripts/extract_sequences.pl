#!/usr/bin/env perl

# MANUAL FOR extract_sequences.pl

=pod

=head1 NAME

extract_sequences.pl -- extracts sequences from a FASTA file

=head1 SYNOPSIS

 extract_sequences.pl --lookup=/Path/to/lookup.txt --fasta=/Path/to/infile.fasta --out=/Path/to/output.fasta [--inverse]
                     [--help] [--manual]

=head1 DESCRIPTION

 Extracts sequences that are in the lookup file from a FASTA file.
 Does the inverse with the --inverse flag
 
=head1 OPTIONS

=over 3

=item B<-f, --fasta>=FILENAME

Input file in FASTA format. (Required) 

=item B<-o, --out>=FILENAME

Output file in FASTA format. (Required) 

=item B<-l, --lookup>=FILENAME

Input lookup file. One seq per line. (Required)

=item B<-v, --inverse>

Perform the inverse operation. (Optional)

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
my($fasta,$out,$lookup,$inverse,$help,$manual);

GetOptions (	
				"f|fasta=s"	=>	\$fasta,
                                "l|lookup=s"    =>      \$lookup,
				"o|out=s"	=>	\$out,
                                "v|inverse"     =>      \$inverse,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual) || pod2usage({-exitval => 2, -verbose => 1, -output => \*STDERR});

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --fasta not found.\n\n", -exitval => 2, -verbose => 1)  if (! $fasta );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $out );
pod2usage( -msg  => "\n\n ERROR!  Required argument --lookup not found.\n\n", -exitval => 2, -verbose => 1)  if (! $lookup );

my %Lookup;
my $print_flag = 0;

open(my $lfh,'<', "$lookup") || die "\n\n Cannot open the input file: $lookup\n\n";
if ($lookup =~ m/\.btab$/) { ## If this is a BTAB file we need to grab the query field
    while(<$lfh>) {
	chomp;
	my @a = split(/\t/, $_);
	my $query_basename = get_basename($a[0]);
	$Lookup{$query_basename} = 0;
    }
}
else { ## if this is a file with just one field (i.e. a list of queries we want)
    while(<$lfh>) {
        chomp;
	my @a = split(/\t/, $_);
	if (scalar(@a) > 1) { die "\n The lookup file needs to either be a .btab file or a one-column list file. This file has: " . scalar(@a) . " fields\n"; }
	my $query_basename = $_;
        $Lookup{$query_basename} = 0;
    }
}
close($lfh);

check_fasta($fasta);

open(my $ofh,'>',"$out") || die "\n Cannot write to the file: $out\n";
open(my $ffh,'<',"$fasta") || die "\n Cannot open the FASTA file: $fasta\n";
while(<$ffh>) {
    chomp;
    if ($_ =~ m/^>/) {
	$print_flag = 0; ## Reset when you see a header
	my $h = $_;
	$h =~ s/^>//;
	$h =~ s/ .*//;
	my $base_h = get_basename($h);
	if ($inverse) {
	    unless (exists $Lookup{$base_h}) {
		print $ofh $_ . "\n";
		$print_flag = 1; ## print the next lines too!
	    }
	}
	else {
	    if (exists $Lookup{$base_h}) {
		print $ofh $_ . "\n";
		$print_flag = 1; ## print the next lines too!
	    }
	}
    }
    elsif ($print_flag == 1) {
	print $ofh $_ . "\n";
    }
}
close($ffh);
close($ofh);

exit 0;

sub get_basename
{
    ## Input: An seqscreen unambigous sequence ID
    ## Output: A sequence ID without the unambig ID
    my $s = $_[0];
    $s =~ s/_unambig_\d+$//;
    return $s;
}

sub check_fasta
{
    ## Input: the input FASTA file
    ## Output: Nothing, just exit if two cases are met.
    my $file = $_[0];
    if ( -z $file ) { die "\n Error: The input FASTA file appears to be empty\n"; } ## If the file is empty then just exit 0
    my $seq_count=0;
    open(my $fh, '<', $file) || die "\n Error: Cannot open the fasta file $file\n";
    while(<$fh>) {
        chomp;
        if ($_ =~ m/^>/) { $seq_count++; }
    }
    close($fh);
    if ($seq_count == 0) { die "\n Error: the input FASTA file does not appear to be in FASTA format as there are no header lines.\n\n"; }
    return 0;
}
