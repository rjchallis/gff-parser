#!/usr/bin/perl -w

use strict;
use GFFTree::NCBI;
use Test::More tests => 1;


my $gff = GFFTree::NCBI->new({});

# Test 1
ok( defined $gff, 'new()' );
