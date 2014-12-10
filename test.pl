#!/usr/bin/perl -w

use strict;

use Cwd 'abs_path';
my $dir = abs_path($0);
$dir =~ s/\/[^\/]+$//;

# generic tests
print "\nRunning generic tests:\n";
system "$^X $dir/t/generic.pl $dir/t/generic.gff3";

# ncbi module tests
print "\nTesting module 'GFFTree::NCBI':\n";
system "$^X $dir/t/ncbi.pl $dir/t/ncbi.gff3";

print "\n";