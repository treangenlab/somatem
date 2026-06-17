#!/usr/bin/env perl

# MANUAL FOR remove_outliers.pl

=pod

=head1 NAME

remove_outliers.pl -- remove the outliers that outlier dectection found

=head1 SYNOPSIS

 remove_outliers.pl --outliers /Path/to/outliers.txt --btab /Path/to/blastn.btab --out=/Path/to/output.btab
                     [--help] [--manual]

=head1 DESCRIPTION

 Gets rid of all the outliers that outlier detection found and writes them to a new btab file.
 
=head1 OPTIONS

=over 3

=item B<-ol, --outliers>=FILENAME

Outlier detection file in tsv format (Required).

=item B<-b, --btab>=FILENAME

The blast output that went in to outlier detection (Required).

=item B<-o, --out>=FILENAME

Output file in blast tabular format. (Required) 

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

Copyright 2018 Daniel Nasko.  
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
my($outlier_file,$btab_file,$outfile,$help,$manual);

GetOptions (	
				"ol|outliers=s"	=>	\$outlier_file,
                                "b|btab=s"      =>      \$btab_file,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual) || pod2usage({-exitval => 0, -verbose => 1, -output => \*STDERR});

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --outliers not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outlier_file );
pod2usage( -msg  => "\n\n ERROR!  Required argument --btab not found.\n\n", -exitval => 2, -verbose => 1)  if (! $btab_file );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

## Global variables
my %Outliers; ## 2D dictionary to hold the outlier detection info

open(my $fh, '<', $outlier_file) || die "\n\n Cannot open the input file: $outlier_file\n\n";
while(<$fh>) {
    chomp;
    my @a = split(/\t/, $_);
    unless ($a[0] eq "#query_sequence") {
	unless ($a[1] eq "NA") {
	    my @b = split(/;/, $a[1]);
	    foreach my $i (@b) {
		$Outliers{$a[0]}{$i} = 1;
	    }
	}
    }
}
close($fh);

open(my $ofh, '>', $outfile) || die "\n\n Error: Cannot write to the output file: $outfile\n\n";
open(my $bfh, '<', $btab_file) || die "\n\n Cannot open the input file: $btab_file\n\n";
while(<$bfh>) {
    chomp;
    my @a = split(/\t/, $_);
    if (exists $Outliers{$a[0]}) {
	if (exists $Outliers{$a[0]}{$a[1]}) {
	    print $ofh $_ ."\n";
	}
    }
    else {
	print $ofh $_ . "\n";
    }
}
close($fh);
close($ofh);

exit 0;
