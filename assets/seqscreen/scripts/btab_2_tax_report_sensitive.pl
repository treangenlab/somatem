#!/usr/bin/env perl

# MANUAL FOR btab_2_tax_report.pl

=pod

=head1 NAME

btab_2_tax_report.pl -- parse btab files into seqscreen taxonomy report format

=head1 SYNOPSIS

 btab_2_tax_report.pl --blastx=/Path/to/blastx.btab --blastn=/Path/to/blastn.btab --out=/Path/to/output.txt [--cutoff=5]
                     [--help] [--manual]

=head1 DESCRIPTION

 Convert btab files to a seqscreen-style report. Requires BLASTn, BLASTx tabular results.
 
=head1 OPTIONS

=over 3

=item B<-x, --blastx>=FILENAME

Input BLASTx results file in tabular format. (Required) 

=item B<-n, --blastn>=FILENAME

Input BLASTn results file in tabular format. (Required)

=item B<-o, --out>=FILENAME

Output file in tabular format. (Required) 

=item B<-c, --cutoff>=INT

Report Tax ID's that are within --cutoff percent of the top-hit bit score.
--cutoff should be an integer from 0 to 50. E.g. 5 would mean report all Tax IDs
with hits that are within 5% of the top-hit bit score. (Default=0).

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
my($blastx,$blastn,$outfile,$help,$manual);

##ARGS WITH DEFAULT
my $cutoff=0;

GetOptions (	
				"x|blastx=s"	=>	\$blastx,
                                "n|blastn=s"    =>      \$blastn,
                                "c|cutoff=s"    =>      \$cutoff,
				"o|out=s"	=>	\$outfile,
                                "h|help"	=>	\$help,
                                "m|manual"	=>	\$manual) || pod2usage({-exitval => 2, -verbose => 1, -output => \*STDERR});

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --blastx not found.\n\n", -exitval => 2, -verbose => 1)  if (! $blastx );
pod2usage( -msg  => "\n\n ERROR!  Required argument --blastn not found.\n\n", -exitval => 2, -verbose => 1)  if (! $blastn );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);
pod2usage( -msg  => "\n\n ERROR!  Argument --cutoff between be between 0 and 50, not $cutoff\n\n", -exitval => 2, -verbose => 1)  if ($cutoff > 50 || $cutoff < 0);
$cutoff = (100 - int($cutoff)) / 100; ## transforming the cut off to a number we can just multiply against the top-hit. e.g. 5 becomes 0.95.

my %TopHit;  ## holds the top-hit bit score values
my %Results; ## holds the set of results based on the cut-off value
my %Source;  ## holds the source of the results (blastn or blastx)
my %Confidence;  ## holds the confidence value of the results (blastn or blastx)

## Loop through both BTAB files to identify the top hit E value for each
## query sequence
find_top_hits_in_btab_file($blastx); ## BLASTx last

## Now loop through the BTABs again and collect results after considering the cut-off
if (-e $blastn) {
    find_top_hits_in_btab_file($blastn); ## BLASTn next
    open(my $blastn_hand,'<', "$blastn") || die "\n Cannot open the file: $blastn\n";
    while(<$blastn_hand>) {
        chomp;
        my @a = split(/\t/, $_);
        my $query = base_query($a[0]);
        my $maxbitscore = $a[7];
        my $bitscore = $a[11]; 
        my $taxid = $a[14]; ## Tax ID's are in a different spot for BLASTn vs. BLASTx results!!!
        my $bitscore_cutoff = $cutoff * $TopHit{$query};
        if ($bitscore >= $bitscore_cutoff) {
        if (exists $Results{$query}{$taxid}) {
            if ($Results{$query}{$taxid} < $bitscore) {
            $Results{$query}{$taxid} = $bitscore;

            }
            $Source{$query}{$taxid} = "blastn";
                $Confidence{$query}{$taxid} = $a[2]/100.0;#$bitscore/$maxbitscore;
        }
        else {
            $Results{$query}{$taxid} = $bitscore;
            $Source{$query}{$taxid} = "blastn";
                $Confidence{$query}{$taxid} = $a[2]/100.0;#$bitscore/$maxbitscore;
        }
        }
    }
    close($blastn_hand);
}


open(my $blastx_hand,'<', "$blastx") || die "\n Cannot open the file: $blastx\n";
while(<$blastx_hand>) {
    chomp;
    my @a = split(/\t/, $_);
    my $query = base_query($a[0]);
    my $maxbitscore = $a[7];
    my $bitscore = $a[11];
    my $taxid = $a[15]; ## Tax ID's are in a different spot for BLASTn vs. BLASTx results!!!
    $taxid =~ s/.*TaxID=//;
    $taxid =~ s/ .*//;
    my $bitscore_cutoff = $cutoff * $TopHit{$query};
    if ($bitscore >= $bitscore_cutoff) {
	if (exists $Results{$query}{$taxid}) {
            if ($Results{$query}{$taxid} < $bitscore) {
                $Results{$query}{$taxid} = $bitscore;
            }
	    $Source{$query}{$taxid} = "all";
            $Confidence{$query}{$taxid} = $a[2]/100.0;#$bitscore/$maxbitscore;
	}
        else {
            $Results{$query}{$taxid} = $bitscore;
            $Source{$query}{$taxid} = "blastx";
            $Confidence{$query}{$taxid} = $a[2]/100.0;#$bitscore/$maxbitscore;
        }
    }
}
close($blastx_hand);

open(my $ofh,'>', "$outfile") || die "\n Cannot write to: $outfile\n";
print $ofh join("\t", "#query","taxid","source", "confidence") . "\n";
foreach my $qry (sort keys %Results) {
    print $ofh $qry . "\t";
    my @a; ## tmp array to hold all tax ids this query could be
    my @b; ## tmp array to hold all sources the tax ids came from (blastn / blastx)
    my @c; ## tmp array to hold all confidence values
    foreach my $tid ( sort { $Results{$qry}{$b} <=> $Results{$qry}{$a} } keys %{$Results{$qry}} ) { ## Sort taxid hits for a query from largest bit score to smallest
	push(@a, $tid);
	push(@b, $Source{$qry}{$tid});
	push(@c, $Confidence{$qry}{$tid});
    }
    print $ofh join(",", @a) . "\t" . join(",", @b) . "\t". join(",", @c) . "\n";
}
close($ofh);

exit 0;

sub find_top_hits_in_btab_file
{
    ## Input: a btab file to parse through for top hits
    ## Output: will fill up the global dictionary %TopHit with top-hit info
    my $file = $_[0];
    open(my $fh,'<', "$file") || die "\n Cannot open the file: $file\n";
    while(<$fh>) {
	chomp;
	my @a = split(/\t/, $_);
	my $query = base_query($a[0]);
	my $bitscore = $a[11];
	unless ($bitscore >= 0) { die "\nError: parsing this BTAB file didnt go as expected. Might be grabbing the wrong field for Bit Score\n\n"; }
	if (exists $TopHit{$query}) {
	    if ($TopHit{$query} < $bitscore) { $TopHit{$query} = $bitscore; }
	}
	else {
	    $TopHit{$query} = $bitscore;
	}
    }
    close($fh);
}

sub base_query
{
    ## Input: a query sequence from the seqscreen pipeline, will have _unambig_ identifiers that we added, but want to remove
    ## Output: the "base" name of the query sequence
    my $s = $_[0];
    $s =~ s/_unambig_\d+$//;
    return $s;
}
