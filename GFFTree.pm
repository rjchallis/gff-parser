#!/usr/bin/perl -w
use strict;
package GFFTree;
use Tree::DAG_Node;
our @ISA=qw(Tree::DAG_Node);
use Encode::Escape::ASCII;


sub new {
    my $class = shift;
    my $options = shift;
    my $self = bless $class->SUPER::new();
    $self->attributes($options);
    return $self;
}

sub parse_file {
	my $node = shift;
	my %ids;
    while (<>){
		next if is_comment($_);
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
			$node = $ids{$attribs->{'Parent'}};
		}
		else {
			# check type to decide what to do with features without parents
		}
		$ids{$attribs->{'ID'}} = $node->new_daughter({%attributes,%$attribs});
		$ids{$attribs->{'ID'}}->name($attribs->{'ID'});
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


1;