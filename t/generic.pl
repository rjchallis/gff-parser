#!/usr/bin/perl -w

use strict;
use Test::More tests => 20;
use GFFTree;


my $gff = GFFTree->new({});
# Test 1
ok( defined $gff, 'new()' );

is($gff->add_expectation('exon','<[_start,_end]','PARENT','warn'), 1, 'add_expectation()');
$gff->add_expectation('mrna','hasParent','gene','force');
$gff->add_expectation('trna','hasParent','gene','force');
is($gff->multiline('cds'),1,'multiline()');

is($gff->is_multiline('cds'),1,'is_multiline()');
$gff->multiline('cdna_match');
$gff->multiline('trna');

is($gff->name('GFF'), 'GFF', 'name()');

is($gff->map_types({'initial' => 'exon', 'internal' => 'exon', 'terminal' => 'exon'}),3,'map_types()');

$gff->lacks_id('all','make');
is($gff->lacks_id(),'make','lacks_id()');

$gff->undefined_parent('make');
is($gff->undefined_parent(),'make','undefined_parent()');

$gff->separator('\t');

my $comment_count = $gff->has_comments('#',['\[','\]']);
#my $comment_count = $gff->has_comments();
print "file has $comment_count comment types\n";

ok($gff->parse_file(),'parse_file()');

# Test 10
ok($gff->validate_all('exon'),"validate_all('exon')");

is($gff->order_features('gene'),19,"order_features('gene')");

is(my @genes = $gff->by_type('gene'),19,"by_type('gene') as array");

ok(my $gene = $gff->by_type('gene'),"by_type('gene') as scalar");

ok($gene = $gff->by_id('gene1'),"by_id()");

is($gene->next_feature('exon')->name,'id20',"next_feature('exon') - with multiple exons");

$gene = $gff->by_id('gene2');
is($gene->next_feature('exon')->name,'id31',"next_feature('exon') - reverse strand");

is(length $gene->as_string(), 136 ,"\$gene->as_string()");

my $cds = $gff->by_id('cds2');
is($cds->{attributes}->{'_start_array'}[0],2227,"\$cds->{attributes}->{'_start_array'}");

is(length $cds->as_string(), 609 ,"\$cds->as_string()");

$gene = $gff->by_id('gene9');
$gene->fill_gaps('cds','5utr','before');

# Test 20
is($gene->next_feature('5utr')->_length(),109,"fill_gaps('exon','5utr','external')");
$gff->validate_all('trna');
$gene = $gff->by_id('gene1317');
print $gene->as_string();
while (my $trna = $gene->next_feature('trna')){
	print $trna->as_string();

}



my $mrna = $gff->by_id('mRNA00002');
$mrna->validate();
print $mrna->as_string();
my @exons = $mrna->order_features('exon');
print scalar @exons,"\n";
while (my $exon = $mrna->next_feature('exon')){
	print $exon->as_string();

}
print "\n";
print $mrna->mother->as_string();
while (my $exon = $mrna->mother->next_feature('exon')){
	print $exon->as_string(1);

}

my @match = ($gff->by_type('match'),$gff->by_type('cdna_match'));
print $match[0]->{attributes}->{Target},"\n";
print $match[1]->{attributes}->{Target},"\n";

__END__
	# next_feature is broken for exon so really need to include seq_name in any ordering - maybe make seq_regions become parents
$gene = $gff->by_name('gene16');
print $gene->as_string();
$gene->validate;
print $gene->as_string();
print $gene->mother()->as_string();
$gene->unlink_from_mother();
$gff->add_daughter($gene);
$gene->validate;
print $gene->as_string();
print $gene->mother()->as_string();

