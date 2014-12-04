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
			$attributes{'_type'} = $data->[2];
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
  Example  : my $node = by_id('id');
=cut

	sub by_id  {
		my $id = pop;
		return $ids{$id};
	}
	
=head2 by_start
  Function : Fetch an arrayref of nodes start position
  Example  : my $node_arrayref = by_start('scaf1','exon',432);
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
			@{$features{$type}} = by_type($self,$type);
			# order features by start_position (and seq_name);
			
			push @{$features{$type}},0;
		}
		return shift @{$features{$type}};
	}

}

=head2 make_introns
  Function : Experimental feature to fill in gaps between exons, needs to become more 
             generic and to actually create new features
  Example  : $gene->make_introns();
=cut

sub make_introns {
	my $self = shift;
	my @exons = by_attribute($self,'_type','exon');
	return unless $exons[0];
	print $exons[0]->{name},"\n";
	my $parent = $exons[0]->mother();
	# TODO: sort exons by start value
	# for now, assume exons are sequential
	my @introns;
	if ($exons[0]->{attributes}->{_start} <= $exons[-1]->{attributes}->{_start}){
    	for (my $i = 0; $i < @exons; $i++){
    		if ($i == 0){
    			if ($exons[$i]->{attributes}->{_start} > $parent->{attributes}->{_start}){
    				push @introns,$i;
    			}
    		}
    		else {
    			if ($exons[$i]->{attributes}->{_start} > $exons[($i-1)]->{attributes}->{_end} + 1){
    				push @introns,$i;
    			}
    		}
    		if ($i == @exons - 1) {
    			if ($exons[$i]->{attributes}->{_end} < $parent->{attributes}->{_end}){
    				push @introns,$i;
    			}
    		}
    	}
    }
    print "$self->{name}\t@introns\t",scalar @exons,"\n" if @introns;
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
			push @{$expectations{$type[$t]}},{'relation' => $relation, 'alt_type' => $alt_type, 'flag' => $flag};
		}
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
		# 		- need to generate a new feature with a new id (use by_id to check generated id is not already in use)
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


sub is_comment {
	return 1 if $_ =~ m/^#/;
}

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

sub _seq_name {
    my $node = shift;
    my $val = shift;
    $node->attributes->{_seq_name} = $val if $val;
    return $node->attributes->{_seq_name};
}

sub _length {
    my $node = shift;
    my $length = $node->attributes->{_end} - $node->attributes->{_start} + 1;
    return $length;
}

sub _phase {
    my $node = shift;
    my $val = shift;
    $node->attributes->{_phase} = $val if defined $val && length $val > 0;
    return $node->attributes->{_phase};
}

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


sub by_type {
    my ($self, $type) = @_;
    return by_attribute($self,'_type',$type);
}


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