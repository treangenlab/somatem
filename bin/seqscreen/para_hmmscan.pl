#!/usr/bin/env perl

# MANUAL FOR para_hmmscan.pl

=pod

=head1 NAME

para_hmmscan.pl -- embarasingly parallel HMMscan

=head1 SYNOPSIS

 para_hmmscan.pl --query=/Path/to/infile.fasta --db=/Path/to/db.hmm --out=/Path/to/output.tab [--evalue=1e-3] [--threads=4]
                     [--help] [--manual]

=head1 DESCRIPTION

=head1 OPTIONS

=over 3

=item B<-q, --query>=FILENAME

Input peptide query file in FASTA format. (Required) 

=item B<-d, --db>=FILENAME

Input peptide HMM DB. The .hmm file needs to have been run through hmmpress. (Required)

=item B<-o, --out>=FILENAME

Path to output tab file. (Required)

=item B<-e, --evalue>=INT

E-value. (Default = 10)

=item B<-t, --threads>=INT

Number of CPUs to use. (Default = 1)

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
use threads;
use FindBin;
use Cwd 'abs_path';
my $script_working_dir = $FindBin::Bin;

#ARGUMENTS WITH NO DEFAULT
my($query,$db,$out,$help,$manual);
my $threads = 1;
my $evalue = 10;
my @THREADS;

GetOptions (	
				"q|query=s"	=>	\$query,
                                "d|db=s"        =>      \$db,
                                "o|out=s"       =>      \$out,
                                "e|evalue=s"    =>      \$evalue,
                                "t|threads=i"   =>      \$threads,
             			"h|help"	=>	\$help,
				"m|manual"	=>	\$manual)  || pod2usage({-exitval => 2, -verbose => 1, -output => \*STDERR});

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
$threads = int($threads);
pod2usage( -msg  => "\n\n ERROR!  Required arguments --query not found.\n\n", -exitval => 2, -verbose => 1)  if (! $query );
pod2usage( -msg  => "\n\n ERROR!  Required arguments --db not found.\n\n", -exitval => 2, -verbose => 1)     if (! $db );
pod2usage( -msg  => "\n\n ERROR!  Required arguments --out not found.\n\n", -exitval => 2, -verbose => 1)    if (! $out );
pod2usage( -msg  => "\n\n ERROR!  Required arguments --threads must be an integer > 0.\n\n", -exitval => 2, -verbose => 1) if ( $threads < 1);
pod2usage( -msg  => "\n\n ERROR!  Required arguments --threads must be a number >= 0.\n\n", -exitval => 2, -verbose => 1)  if ( $evalue < 0 );
pod2usage( -msg  => "\n\n ERROR!  The input FASTA file cannot be compressed.\n\n", -exitval => 2, -verbose => 1)  if ( $query =~ m/\.gz$/ );

my $program = "hmmscan";
my $splitby = 2; ## How many threads should each split get?
my @chars = ("A".."Z", "a".."z");
my $rand_string;
$rand_string .= $chars[rand @chars] for 1..8;
my $outdir=dirname($out);
my $tmp_dir = $outdir . "/$program" . "_tmp_" . $rand_string; ## a temporary working directory with a unique name

## Check that phmmer in installed on this machine and in PATH
my $PROG = `which $program`; unless ($PROG =~ m/$program/) { die "\n\n ERROR: External dependency '$program' not installed in system PATH\n\n";}
my $date = `date`;
print STDERR " Using $threads threads\n";
print STDERR " Using this $program: $PROG Beginning: $date\n";

## If only 1 thread is selected, just run the program as-is...
if ($threads == 1) {
    my $returnCode = system( "$program --cpu 1 -o $out.txt --tblout $out.tab -E $evalue $db $query" );
    if ($returnCode != 0) {
	die "\n $program failed and exited with code: $returnCode\n\n";
    }
}
else {
    print `mkdir -p $tmp_dir`;
    print `chmod 700 $tmp_dir`;
    my %CoreDist = distribute_cores($threads, $splitby); ## split the input file up based on the number of CPUs the user allowed
    my $nfiles = keys %CoreDist;
    my $seqs = count_seqs($query);
    my $seqs_per_thread = seqs_per_thread($seqs, $nfiles);
    $nfiles = split_multifasta($query, $tmp_dir, "split", $seqs_per_thread, $threads);
    `mkdir -p $tmp_dir/result_splits`;
    for (my $i=1; $i<=$nfiles; $i++) {
	my $hmmer_exe = "$program --cpu $splitby -o $tmp_dir/result_splits/split.$i.txt --tblout $tmp_dir/result_splits/split.$i.tab -E $evalue $db $tmp_dir/split-$i.fsa";
	push (@THREADS, threads->create('task',"$hmmer_exe"));
    }
    foreach my $thread (@THREADS) {
	$thread->join();
    }
    `cat $tmp_dir/result_splits/*.txt > $out.txt`; ## Concatenate the txt output files into the final output file
    `cat $tmp_dir/result_splits/*.tab > $out.tab`; ## Same, but for the .tab files
    `rm -rf $tmp_dir`; ## remove the temporary directory
}
$date = `date`;
print STDERR "\n Parallel $program complete: $date\n";

exit 0;

sub task
{
    ## Input: a bunch of commands to send to the system
    ## Output: None.
    my $returnCode = system( @_ );
    if ($returnCode != 0) {
	die "\n Error: $program exitted with value: $returnCode\n\n";
    }
}

sub count_seqs
{
    ## Input: A FASTA file
    ## Output: the number of sequences in the FASTA file
    my $q = $_[0];
    my $s = 0;
    open(my $fh,'<',"$q") || die "\n Cannot open the file: $q\n";
    while(<$fh>) {
	chomp;
	if ($_ =~ m/^>/) { $s++; }
    }
    close($fh);
    return $s;
}

sub split_multifasta
{
    ## Inputs:
    my $q       = $_[0]; ## Input FASTA file
    my $working = $_[1]; ## Termporary working directory
    my $prefix  = $_[2]; ## Prefix you want to give the split FASTA files
    my $spt     = $_[3]; ## Number of sequences in each split FASTA
    my $nfiles  = $_[4]; ## Number of files that the FASTA should be split into
    ## Output: Will split the input file into nfiles with spt sequences in each. Return number of files.
    my $j=0;
    my $fileNumber=1;
    print `mkdir -p $working`;
    open(my $fh,'<',"$q") || die "\n Cannot open the file: $q\n";
    open(my $ofh,'>',"$working/$prefix-$fileNumber.fsa") or die "Error! Cannot create output file: $working/$prefix-$fileNumber.fsa\n";
    while(<$fh>) {
        chomp;
        if ($_ =~ /^>/) { $j++; }
        if ($j > $spt && $fileNumber < $nfiles) { #if time for new output file
            close($ofh);
            $fileNumber++;
            open($ofh, '>', "$working/$prefix-$fileNumber.fsa") or die "Error! Cannot create output file: $working/$prefix-$fileNumber.fsa\n";
            $j=1;
        }
        print $ofh $_ . "\n";
    }
    close($fh);
    close($ofh);
    return $fileNumber;
}

sub seqs_per_thread
{
    ## Input: Number of sequences in the FASTA file AND number of threads the user provided
    ## Output: The number of sequences each thread to analyze. I.e. the number of sequences each split should get.
    my $s = $_[0];
    my $t = $_[1];
    my $seqs_per_file = $s / $t;
    if ($seqs_per_file =~ m/\./) {
        $seqs_per_file =~ s/\..*//;
        $seqs_per_file++;
    }
    return $seqs_per_file;
}

sub distribute_cores
{
    ## Inputs:
    my $t = $_[0];  ## The number of CPUs the user is allowing the program to use
    my $by = $_[1]; ## The number of threads each split is allowed to use
    ## Output a hash telling us how many sequences should be in each of the split up files
    my %Hash;
    my $nsplits = calc_splits($t, $by); ## How many temporary files should we create based on number of CPUs/threads
    my $file=1;
    for (my $i=1; $i<=$t; $i++){
        $Hash{$file}++;
        if ($file==$nsplits) { $file = 0;}
        $file++;
    }
    return %Hash;
}

sub calc_splits
{
    ## Inputs: (1) Number of CPUs the user allowed (2) the number of threads for each CPU
    ## Output: The number of temporary files we need by dividing these numbers
    my $t = $_[0];
    my $by = $_[1];
    my $n = roundup($t/$by);
    return $n;
}

sub roundup {
    ## Input: a number
    ## Output: round that number up
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}
