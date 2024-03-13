#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Data::Dumper;


# Variables
my $fastq_file;
my $k_mers;
my %seqs = ();  #hash of hashes
my %edges = ();
my %eulerian_nodes = ();
my %seqs2 = ();
my %skip_sequences = ();

# Check if the correct number of command-line arguments is provided
if (@ARGV != 2) {
    die "Incorrect number of arguments passed";
}

# Check if the file is compressed (.gz)
if ($ARGV[0] =~ /\.gz$/i) {
    my $fastq_in = $ARGV[0];
    gunzip $fastq_in => "fastq_file.fq" or die "gunzip failed: $GunzipError\n";
    $fastq_file = "fastq_file.fq";
} elsif ($ARGV[0] =~ /\.fq$/i) {
    $fastq_file = $ARGV[0];
} else {
    die "Unrecognised file type";
}

# Check if the provided k-mer size is valid
if (($ARGV[1] > 0) && ($ARGV[1] < 200)) {
    $k_mers = $ARGV[1];
} elsif ($ARGV[1] >= 200) {
    die "K-mers size provided is too big";
} else {
    die "Invalid K-mer argument";
}

###################
###     Main    ###
###################

parse_fastq($fastq_file, \%seqs, \%seqs2, $k_mers, \%edges, \%eulerian_nodes);

skips_hash(\%skip_sequences, \%eulerian_nodes); 

simplify_graph(\%eulerian_nodes, \%seqs, \%skip_sequences);

my $final_path = '';
for my $key (keys %eulerian_nodes) {
    if ($eulerian_nodes{$key} > 0) {
        my $current = find_longest_path($key, \%eulerian_nodes);
        if (length($current) > length($final_path)) {
            $final_path = $current;
        }
    }
}

# Save final sequence as a file
my $file_path = "final_sequence.fasta";

open my $file_handle, '>', $file_path or die "Cannot open file '$file_path' for writing: $!";

print $file_handle ">" . $ARGV[0] . "\n";
print $file_handle $final_path;

close $file_handle;

print "Final sequence has been saved to $file_path\n";

###################
###  Functions  ###
###################

sub parse_fastq {
    my ($fastq_file, $seqs, $seqs2, $kmer_size, $edges, $nodes) = @_;

    open my $fh, '<', $fastq_file or die "Cannot open file: $!\n";

    my $nlines = 0;

    # Iterate through each sequence in the FASTQ file
    while (my $line = <$fh>) {
        chomp $line;

        # Check if the line is the start of a new sequence
        if ($line =~ /^@/) {
            $nlines++;
            my $seq_id = $line;

            # Read the next three lines to get the sequence, a plus sign, and quality
            my $seq = <$fh>;

            chomp($seq);
            # Store the sequence in the hash

            for (my $i = 0; $i <= length($seq) - $kmer_size; $i++) {
                my $kmer = substr($seq, $i, $kmer_size);
                my $left_k = substr($kmer, 0, $kmer_size - 1);
                my $right_k = substr($kmer, 1, $kmer_size);

                $seqs->{$left_k}{$right_k}++;
                $seqs2->{$left_k}{$right_k}++;
                $edges->{ $left_k . "::" . $right_k }++;
                $nodes->{$left_k}++;
                $nodes->{$right_k}--;
            }
        }
    }

    # Close the FASTQ file handle
    close $fh;

}   # parse_fastq

sub simplify_graph {
    my ($nodes, $graph, $skips) = @_;

    foreach my $node (keys %$nodes) {
        if ($nodes->{$node} > 1) {
            foreach my $rightk (keys %{ $graph->{$node} }) {
                $nodes->{$node} -= $graph->{$node}{$rightk};
                delete $graph->{$node}{$rightk};
            }

            $graph->{$node}{substr($skips->{$node}, ($k_mers - 1) * -1)} += $nodes->{$node} * -1;
        }
    }
} # simplify_graph

sub skips_hash {
    my ($skips, $nodes) = @_;

    foreach my $node (keys %$nodes) {
        if ($nodes->{$node} > 1) {
            $skips->{$node} = find_longest_path($node, $nodes);
        }
    }
} # skips_hash

sub find_longest_path {
    my ($start, $nodes) = @_;

    # 1. generate all possible paths across the bubble
    my @paths = ();
    while (@paths == 0 || length($paths[-1]) > $k_mers - 1) {
        push @paths, traverse_graph($start, \%seqs, $nodes, \%skip_sequences);
    }
    
    # 2. Iterate through the paths and find the longest path
    # Alternatively, select the path based on other criteria
    # e.g., path coverage
    my $longest_length = 0;
    my $longest_path;

    foreach my $path (@paths) {
        my $length = length($path);

        if ($length > $longest_length) {
            $longest_length = $length;
            $longest_path = $path;
        }
    }
    
    # return the selected path as a string
    return $longest_path;
} # find_longest_path

sub traverse_graph {
    my ($start, $graph, $nodes, $skips) = @_;

    # This function returns one available route across a bubble.
    my @stack = ($start);
    my @visited = ();
    my $stop = 0;

    while (!$stop) {
        my $current = pop @stack;

        if (defined $current) { # Check if $current is defined
            push @visited, $current;

            if (exists $graph->{$current} && $nodes->{$current} > -1) {
                # Still in the bubble
                my $subnode_found = 0;

                foreach my $subnode (keys %{$graph->{$current}}) {
                    my $count = $graph->{$current}{$subnode};
                    

                    if ($count > 0) {
                        push @stack, $subnode;
                        # if (@visited > 2) {
                        #     $graph->{$visited[-2]}{$visited[-1]}++;
                        # }
                        $graph->{$current}{$subnode}--;
                        $subnode_found = 1;
                        last;
                    }
                }

                unless ($subnode_found) {
                    $stop = 1;
                }
            } else {
                $stop = 1;
            }
        } else {
            $stop = 1; # Stop if $current is undefined
        }
    }

    my $skip = "";
    for my $j (0..$#visited){
        if ($j < $#visited) {

            if (exists($skips->{$visited[$j]})) {
                $skip = $skip . $skips->{$visited[$j]};
            } else {
                $skip = $skip . substr($visited[$j], 0, 1);
            }
            
        } else {
            $skip = $skip . $visited[$j];
        }
    }


    return $skip;
} # traverse_graph