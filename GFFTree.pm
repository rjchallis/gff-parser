#!/usr/bin/perl -w

=head1 NAME

GFFTree

=head1 SYNOPSIS

  my $gff = GFFTree->new({});
  $gff->name('GFF');
  $gff->parse_file();
  $gff->add_expectation('mrna','hasParent','gene','find');
  $gff->add_expectation('mrna','<[_start,_end]','SELF','warn');
  $gff->add_expectation('mrna','>=[_start]','PARENT','warn');
  while (my $gene = $gff->next_feature('gene')){
	while (my $mrna = $gene->next_feature('mrna')){
		$mrna->validate;
	}
  }


=head1 DESCRIPTION

This module is intended to be as close as possible to a universal gff3 file parser. 

Calling parse_file causes the gff3 file to be read into a graph structure, allowing 
relationships to be queried and modified during an optional validation step, which can be 
applied to individual features independently.

The parser makes no assumptions about the feature types expected in column 3 nor the 
content of column 9 (other than that it is a set of key-value pairs).  Expected properties 
of, and relationships among and between, features may be defined as a set of expectations, 
which can be tested during an optional validation step. The behaviour of the parser on 
finding an unmet assumption (to ignore, warn, die, link to an existing feature or create a 
new feature) can be defined independently for each assumption through the use of flags.

=cut

use strict;
package GFFTree;
use Tree::DAG_Node;
our @ISA=qw(Tree::DAG_Node);
use Encode::Escape::ASCII;

=head2 new
  Function : Creates a new GFFTree
  Example  : $gff = GFFTree->new({});
=cut

sub new {
    my $class = shift;
    my $options = shift;
    my $self = bless $class->SUPER::new();
    $self->attributes($options);
    return $self;
}


# use nesting to allow subs to share and retain access to private variables
{
	my %ids;
	my %starts;
	my %type_map;

=head2 map_types
  Function : Loads a mapping of types to allow treating of features with different types 
             as one
  Example  : $gff->map_types({'internal' => 'exon'});
=cut

	sub map_types {
		my $mapping = pop;
	    foreach my $type (keys %{$mapping}){
	    	$type_map{$type} = $mapping->{$type};
	    }
	    return scalar keys %type_map;
	}


=head2 parse_file
  Function : Reads a gff3 file into a tree structure
  Example  : $gff->parse_file();
=cut

	sub parse_file {
		my $node = shift;
	    while (<>){
			next if is_comment($_);
			my $parent = $node;
			my ($data,$attribs) = parse_gff_line($_);
			my %attributes;
			$attributes{'_seq_name'} = $data->[0];
			$attributes{'_source'} = $data->[1];
			$attributes{'_type'} = $type_map{$data->[2]} ? $type_map{$data->[2]} : $data->[2];
			$attributes{'_start'} = $data->[3];
			$attributes{'_end'} = $data->[4];
			$attributes{'_score'} = $data->[5];
			$attributes{'_strand'} = $data->[6];
			$attributes{'_phase'} = $data->[7];
			if ($attribs->{'Parent'}){
				$parent = $ids{$attribs->{'Parent'}};
			}
			else {
				# check type to decide what to do with features without parents
			}
			$ids{$attribs->{'ID'}} = $parent->new_daughter({%attributes,%$attribs});
			$ids{$attribs->{'ID'}}->name($attribs->{'ID'});
			push @{$starts{$data->[0]}{$data->[2]}{$data->[3]}},$ids{$attribs->{'ID'}};
		}
	}

=head2 by_id
  Function : Fetch a node by unique id
  Example  : $node = by_id('id');
=cut

	sub by_id  {
		my $id = pop;
		return $ids{$id};
	}

=head2 make_id
  Function : generate a unique id for a node
  Example  : $id = $node->make_id('prefix');
=cut

	sub make_id  {
		my $node = shift;
		my $prefix = shift;
		my $suffix = shift;
		$suffix = 0 unless $suffix;
		while (by_id($prefix.$suffix)){
			$suffix++;
		}
		my $id = $prefix.$suffix;
		$node->name($id);
		$node->{attributes}->{'ID'} = $id;
		$ids{$id} = $node;
		return $id;
	}

	
=head2 by_start
  Function : Fetch an arrayref of nodes start position
  Example  : $node_arrayref = by_start('scaf1','exon',432);
=cut

	sub by_start  {
		my $start = pop;
		my $type = pop;
		my $seq_name = pop;
		return $starts{$seq_name}{$type}{$start};
	}

}


# use nesting to allow sub to retain access to private variables

{
	my %features;
	
=head2 next_feature
  Function : Sequentially fetch daughter features of a node by type
  Example  : $gene->next_feature('exon');
=cut

	sub next_feature {
		my ($self, $type) = @_;
		unless ($features{$type} && @{$features{$type}}){
			$self->order_features($type);
		}
		return shift @{$features{$type}};
	}
	
=head2 order_features
  Function : order daughter features of a given  type sequentially
  Example  : $gene->order_feature('exon');
=cut

	sub order_features {
		my ($self, $type) = @_;
		my @unsorted = by_type($self,$type);
		@{$features{$type}} = sort { $a->{attributes}->{_start} <=> $b->{attributes}->{_start} } @unsorted;
		push @{$features{$type}},0;
	}
	
=head2 fill_gaps
  Function : fill in gaps between features, eg make introns based on exons, should check 
  whether such features exist before running this (maybe need a fill_gaps_unless sub)
  Example  : $gene->fill_gaps('exon','intron','internal');
             $gene->fill_gaps('exon','utr','external');
             $gene->fill_gaps('exon','5utr','before');
             $gene->fill_gaps('exon','3utr','after');
=cut

	sub fill_gaps {
		my ($self, $type, $new_type, $location) = @_;
		$self->order_features($type);
		my %attributes;
		$attributes{'_seq_name'} = $self->{attributes}->{_seq_name};
		$attributes{'_source'} = 'GFFTree';
		$attributes{'_type'} = $new_type;
		$attributes{'_score'} = '.';
		$attributes{'_strand'} = $self->{attributes}->{_strand};
		$attributes{'_phase'} = '.';
		$attributes{'Parent'} = $self->name;
		if ($location eq 'internal'){
			if (@{$features{$type}} > 2){
				for (my $i = 1; $i < @{$features{$type}} - 1; $i++){
					$attributes{'_start'} = $features{$type}[($i-1)]->{attributes}->{_end} + 1;
					$attributes{'_end'} = $features{$type}[$i]->{attributes}->{_start} - 1;
					my $node = $self->new_daughter(\%attributes);
					$node->make_id($new_type);
				}
			}
		}
		if ($location eq 'before' || $location eq 'external'){
			if (@{$features{$type}} > 1 && $features{$type}[0]->{attributes}->{_start} > $self->{attributes}->{_start}){
				$attributes{'_start'} = $self->{attributes}->{_start};
				$attributes{'_end'} = $features{$type}[0]->{attributes}->{_start} - 1;
				my $node = $self->new_daughter(\%attributes);
				$node->make_id($new_type);
			}
		}
		if ($location eq 'after' || $location eq 'external'){
			if (@{$features{$type}} > 1 && $features{$type}[-2]->{attributes}->{_end} < $self->{attributes}->{_end}){
				$attributes{'_end'} = $self->{attributes}->{_end};
				$attributes{'_start'} = $features{$type}[-2]->{attributes}->{_end} + 1;
				my $node = $self->new_daughter(\%attributes);
				$node->make_id($new_type);
			}
		}
		# TODO - check for existing features with validation_find()
	}

}


# use nesting to allow subs to share and retain access to private variables
{
	my %expectations;
	# lookup table for flags
	my %actions = ( 'ignore' => \&validation_ignore,
					'warn' => \&validation_warning,
					'die' => \&validation_die,
					'find' => \&validation_find,
					'make' => \&validation_make,
					'force' => \&validation_force,
              );
	
=head2 add_expectation
  Function : Specify conditions that feature-types should meet in order to pass validation 
             expectations can be applied to multiple feature types using the pipe symbol in 
             $Arg[0]. More documentation to be written...
  Example  : $gff->add_expectation('mrna','hasParent','gene','find');
=cut

	sub add_expectation {
		# define expectations for gff validation
		# feature_type,relation,alt_type,flag
		# flags: ignore, warn, die, find, make, force
		# relations: hasParent, hasChild, >, gt, <, lt, ==, eq, >=, <=, !=, ne 
		# mrna hasParent gene 
		# mrna|exon <[start,end] SELF
		# mrna <=[end] PARENT warn
		# exon hasParent mrna|transcript|gene
		# cds hasParent exon
		
		my ($self,$type,$relation,$alt_type,$flag) = @_;
		$type =~ tr/[A-Z]/[a-z]/;
		my @type = split /\|/,$type;
		for (my $t = 0; $t < @type; $t++){
			push @{$expectations{$type[$t]}},{'relation' => $relation, 'alt_type' => $alt_type, 'flag' => $flag} || return;
		}
		return scalar keys %expectations;
	}

=head2 validate
  Function : Test whether a feature meets the conditions defined in any relevant added 
             expectations and respond according to the specified flag
  Example  : $mrna->validate();
=cut
	
	sub validate {
		my $self = shift;
		my $type = $self->{attributes}->{_type};
		$type =~ tr/[A-Z]/[a-z]/;
		if ($expectations{$type}){
			for (my $i = 0; $i < @{$expectations{$type}}; $i++){
				my $hashref = $expectations{$type}[$i];
				if ($hashref->{'relation'} eq 'hasParent'){
					my $message = $type." ".$self->name.' does not have a parent of type '.$hashref->{'alt_type'};
					$actions{$hashref->{'flag'}}->($self,$message,$expectations{$type}[$i]) unless $self->mother->{attributes}->{_type} =~ m/$hashref->{'alt_type'}/i;
				}
				elsif ($hashref->{'relation'} eq 'hasChild'){
					#;
				}
				else {
					my @relation = split /[\[\]]/,$hashref->{'relation'};
					my @attrib = split /,/,$relation[1];
					my $message = $type.' '.$self->name.'->('.$attrib[0].') is not '.$relation[0].' '.$hashref->{'alt_type'}.'->('.$attrib[-1].')';
					my $first = $self->{attributes}->{$attrib[0]};
					my $second = $hashref->{'alt_type'} =~ m/self/i ? $self->{attributes}->{$attrib[-1]} : $self->mother->{attributes}->{$attrib[-1]};
					$actions{$hashref->{'flag'}}->($message) unless compare($first,$second,$relation[0]);
				}
			}
		}
	}
	
=head2 compare
  Function : compare two values based on a specified operator
  Example  : compare(10,20,'<'); # returns true
=cut
	
	sub compare {
		my ($first,$second,$operator) = @_;
		return $first > $second if $operator eq '>';
		return $first gt $second if $operator eq 'gt';
		return $first < $second if $operator eq '<';
		return $first lt $second if $operator eq 'lt';
		return $first == $second if $operator eq '==';
		return $first eq $second if $operator eq 'eq';
		return $first >= $second if $operator eq '>=';
		return $first <= $second if $operator eq '<=';
		return $first != $second if $operator eq '!=';
		return $first ne $second if $operator eq 'ne';
	}
	
	sub validation_ignore {
		# nothing happens
	}
	
	sub validation_warning {
		my $message = pop;
		warn "WARNING: $message\n";
	}
	
	sub validation_die {
		my $message = pop;
		die "ERROR: $message\n";
	}

=head2 validation_find
  Function : find a feature to satisfy an expectation - limited functionality at present
  Example  : validation_find($expectation_hashref);
=cut
	
	sub validation_find {
		# TODO 	- add parent/child relationships
		# 		- handle relationships other than parent
		my $self = shift;
		my $expectation = pop;
		if ($expectation->{'relation'} eq 'hasParent'){
			my @possibles = by_start($self->{attributes}->{'_seq_name'},$expectation->{'alt_type'},$self->{attributes}->{'_start'});
			while (my $parent = shift @possibles){
				if ($self->{attributes}->{'_end'} == $parent->[0]->{attributes}->{'_end'}){
					warn 'found '.$self->name.' a parent '.$expectation->{'alt_type'}.' '.$parent->[0]->name."\n";
				}
			}
		}
	}

=head2 validation_make
  Function : make a feature to satisfy an expectation - to be implemented
  Example  : validation_make($expectation_hashref);
=cut
	
	sub validation_make {
		# TODO 	- everything
		# 		- need to generate a new feature with a new id (use make_id to check generated id is not already in use)
	}
	
=head2 validation_force
  Function : find a feature to satisfy an expectation if possible, otherwise make one - 
             requires other functions to work first...
  Example  : validation_force($expectation_hashref);
=cut

	sub validation_force {
		my $find = validation_find;
		return $find if $find;
		return validation_make;
	}

}

=head2 is_comment
  Function : returns true if a line is a comment
  Example  : is_comment($line);
=cut

sub is_comment {
	return 1 if $_ =~ m/^#/;
}

=head2 parse_gff_line
  Function : splits a line of gff into 8 fields and a key-value hash, escaping encoded 
             characters
  Example  : parse_gff_line($line);
=cut

sub parse_gff_line {
	my @data = split /\t/,$_;
	chomp $data[8];
	my %attribs = split /[=;]/,$data[8];
	pop @data;
	foreach my $key (keys %attribs){
		$attribs{$key} =~ s/\%/\\x/g;
		$attribs{$key} = decode 'ascii-escape', $attribs{$key};
	}
	return \@data,\%attribs;
}

=head2 _seq_name
  Function : get/set the sequence name for a feature
  Example  : $name = $node->_seq_name();
             $node->_seq_name('new_name');
=cut

sub _seq_name {
    my $node = shift;
    my $val = shift;
    $node->attributes->{_seq_name} = $val if $val;
    return $node->attributes->{_seq_name};
}

=head2 _length
  Function : get the length of a feature
  Example  : $length = $node->_length();
=cut

sub _length {
    my $node = shift;
    my $length = $node->attributes->{_end} - $node->attributes->{_start} + 1;
    return $length;
}

=head2 _phase
  Function : get/set the phase name for a feature
  Example  : $phase = $node->_phase();
             $node->_phase(2);
=cut

sub _phase {
    my $node = shift;
    my $val = shift;
    $node->attributes->{_phase} = $val if defined $val && length $val > 0;
    return $node->attributes->{_phase};
}

=head2 _end_phase
  Function : get/set the end_phase name for a feature - ensembl specific at the moment
  Example  : $phase = $node->_end_phase();
             $node->_end_phase(2);
=cut

sub _end_phase {
    my $node = shift;
    my $val = shift;
    $node->attributes->{_end_phase} = $val if defined $val && length $val > 0;
    unless (defined $node->attributes->{_end_phase}){
    	if ($node->_phase eq '-1'){
    		$node->attributes->{_end_phase} = '-1';
    	}
    	else {
	    	$node->attributes->{_end_phase} = $node->_phase + $node->_length % 3;
    		$node->attributes->{_end_phase} -= 3 if $node->attributes->{_end_phase} > 2;
    	}
    }
    return $node->attributes->{_end_phase};
}

=head2 by_name
  Function : walk through the tree to find a feature by name, works in scalar or array 
             context
  Example  : $node = $gff->by_name('id2');
=cut


sub by_name {
    my ($self, $name) = @_;
    my @found =();
    my $retvalue = wantarray ? 1 : 0;
    $self->walk_down({callback=>sub{
        if ($_[0]->name eq $name) {
            push @found, $_[0];
            return $retvalue;
        }
        1}});
    return wantarray? @found : @found ? $found[0] : undef;
}

=head2 by_type
  Function : calls by_attribute to find an array of features of a given type
  Example  : @nodes = $gff->by_type('exon');
=cut


sub by_type {
    my ($self, $type) = @_;
    return by_attribute($self,'_type',$type);
}

=head2 by_attribute
  Function : returns a scalar or array of features with a given attribute or where an 
             attribute matches a specific value
  Example  : @nodes = $gff->by_attribute('anything');
  			 $node = $gff->by_attribute('anything','something');
=cut

sub by_attribute {
    my ($self, $attrib, $value) = @_;
    my @found =();
    my $retvalue = wantarray ? 1 : 0;
    $self->walk_down({callback=>sub{
        if ($_[0]->attributes->{$attrib} && (!defined $value || $_[0]->attributes->{$attrib} =~ m/^$value$/i)) {
            push @found, $_[0];
            return $retvalue;
        }
        1}});
    return wantarray? @found : @found ? $found[0] : undef;
}

=head2 by_not_attribute
  Function : returns a scalar or array of features without a given attribute or where an 
             attribute does not match a specific value
  Example  : @nodes = $gff->by_attribute('anything');
  			 @nodes = $gff->by_attribute('anything','something');
=cut

sub by_not_attribute {
    my ($self, $attrib, $value) = @_;
    my @found =();
    my $retvalue = wantarray ? 1 : 0;
    $self->walk_down({callback=>sub{
    	if (defined $value){
	        if ($_[0]->attributes->{$attrib} && $_[0]->attributes->{$attrib} !~ m/^$value$/i) {
    	        push @found, $_[0];
       	     	return $retvalue;
       	 	}
       	}
       	else {
       		if (!$_[0]->attributes->{$attrib}) {
            	push @found, $_[0];
            	return $retvalue;
        	}
        }
        1}});
    return wantarray? @found : @found ? $found[0] : undef;
}


1;