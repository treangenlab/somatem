#!/usr/bin/env perl

# MANUAL FOR blast_formatter_parallel.pl

=pod

=head1 NAME

blast_formatter_parallel.pl -- converts the BLAST ASN file into other formats all at once

=head1 SYNOPSIS

 blast_formatter_parallel.pl --asn=/Path/to/infile.asn --out=/Path/to/output --format=5,6,0
                     [--help] [--manual]

=head1 DESCRIPTION

 Converts a BLAST ASN file to one or more BLAST formats in parallel
 
=head1 OPTIONS

=over 3

=item B<-a, --asn>=FILENAME

Input BLAST archive file in ASN.1 format. (Required) 

=item B<-o, --out>=FILENAME

Output file base name. (Required)

=item B<-f, --format>=OUTFMT

Output formats to write, multiple can be specified with CSV. (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.
threads


=head1 AUTHOR

Written by Daniel Nasko, 
Center for Bioinformatics and Computational Biology, University of Maryland.

=head1 REPORTING BUGS

Report bugs to dnasko@umiacs.umd.edu

=head1 COPYRIGHT

Copyright 2017 Daniel Nasko.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use threads;

#ARGUMENTS WITH NO DEFAULT
my($asn,$outroot,$format,$help,$manual);

GetOptions (	
				"a|asn=s"	=>	\$asn,
				"o|out=s"	=>	\$outroot,
                                "f|format=s"    =>      \$format,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual) || pod2usage({-exitval => 2, -verbose => 1, -output => \*STDERR});

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --asn not found.\n\n", -exitval => 2, -verbose => 1)  if (! $asn );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outroot);
pod2usage( -msg  => "\n\n ERROR!  Required argument --format not found.\n\n", -exitval => 2, -verbose => 1)  if (! $format);


my @THREADS;
my @Formats = split(/,/, $format); ## gather desired formats into the array
my %Extensions = ( ## Dictionary of comman suffixes for output formats
    0 => '.raw',
    5 => '.xml',
    6 => '.btab' );

foreach my $i (@Formats) {
    my $form_actual = $i;
    $form_actual =~ s/\"//g; $form_actual =~ s/ .*//; ## sometimes you can pass in formats that aren't just a number, e.g. "6 std salltitles".
                                                      ## this code gets us to what the number is
    my $blast_exe;
    my $output_file = "$outroot" . "$Extensions{$form_actual}";
    open(my $fh, '>', "$output_file") || die "\n Error: you cannot write to the output area you chose: $outroot\n\n";
    close($fh);
    if (exists $Extensions{$form_actual}) { ## validate that this is an acceptable blast format value
	$blast_exe = "diamond view -a $asn --top 5 --out $output_file --outfmt $i";
    }
    else { die "\n Error unkonwn format passed: $i\n"; }
    push (@THREADS, threads->create('task',"$blast_exe")); ## run the blast_formatter commands in parallel on the system
}

foreach my $thread (@THREADS) {
    $thread->join();
}

exit 0;

sub task
{
    system( @_ );
}
