#!/usr/bin/perl -w

=head1 NAME

GFFTree::NCBI

=head1 SYNOPSIS

  my $gff = GFFTree::NCBI->new({});
  $gff->name('GFF');
  $gff->parse_file();
  while (my $gene = $gff->next_feature('gene')){
	while (my $mrna = $gene->next_feature('mrna')){
		$mrna->validate;
	}
  }


=head1 DESCRIPTION

sets NCBI specific expectation for working with a gff file.

=cut

use strict;
package GFFTree::NCBI;
use base qw(GFFTree);

=head2 new
  Function : Creates a new GFFTree::NCBI
  Example  : $gff = GFFTree::NCBI->new({});
=cut

sub new {
    my $self = shift;
    my $options = shift;
    my $gff = GFFTree->new($options);
    $gff->add_expectation('gene','hasParent','region','force');
	$gff->add_expectation('match','hasParent','region','force');
	$gff->add_expectation('cDNA_match','hasParent','region','force');
	$gff->add_expectation('Genomic','hasParent','region','force');
	$gff->add_expectation('mrna','hasParent','gene','find');
	$gff->add_expectation('trna','hasParent','gene','find');
	$gff->add_expectation('transcript','hasParent','gene','find');
	$gff->add_expectation('exon','hasParent','mrna','find');
	$gff->add_expectation('cds','hasParent','mrna','find');
	$gff->multiline('CDS');
	$gff->multiline('tRNA');
	$gff->multiline('cDNA_match');
	$gff->multiline('five_prime_UTR');
	$gff->multiline('three_prime_UTR');
	$gff->multiline('match');
	$gff->add_expectation('cds|exon|mrna|trna|transcript|gene|region','<=[_start,_end]','SELF','warn');
	$gff->add_expectation('cds|exon|mrna|trna|transcript','>=[_start]','PARENT','warn');
	$gff->add_expectation('cds|exon|mrna|trna|transcript','<=[_end]','PARENT','warn');

    
    return $gff;
}

1;

