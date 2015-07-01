#!/usr/bin/perl

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
use warnings;
package GFFTree;
use Tree::DAG_Node;
use Encode::Escape::ASCII;
use IO::Unread;
our @ISA=qw(Tree::DAG_Node);

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
	my %suffices;
	my %by_start;
	my %type_map;
	my %multiline;
	my %parents;
	my $separator = '\t';
	my $has_comments = undef;
	my $lastline;


=head2 separator
  Function : Sets/returns the separator to be used when parsing a gff file
  Example  : $gff->separator('\t');
=cut

	sub separator {
		my $sep = pop;
		if (ref $sep ne 'HASH'){
	    	$separator = $sep;
	    }
	    return $separator;
	}


=head2 has_comments
  Function : treat anything after (optionally between) specified characters as a comment
			 each call will overwrite any previous values
  Example  : $gff->has_comments('#'); # anything after the hash character will be ignored
  			 $gff->has_comments(['[',']']); # anything in square brackets will be ignored
  			 $gff->has_comments('#', ['[',']']); # multiple comments can be specified at once
  			 $gff->has_comments(); # the file has no comments
=cut

	sub has_comments {
		my ($node,@arr) = @_;
		$has_comments = \@arr;
	    return scalar @arr;
	}


=head2 delete_comments
  Function : delete comments from lines (as specified by has_comments)
  Example  : delete_comments($line);
=cut

	sub delete_comments {
		my $line = pop;
		my @array = @$has_comments;
		while (my $char = shift @array){
			if (ref $char && ref $char eq 'ARRAY'){
				my $pattern = $char->[0].'.+?'.$char->[1];
				$line =~ s/$pattern//g;
			}
			else {
				my $pattern = $char.'.+';
				$line =~ s/$pattern//;
			}
		}
		return $line;
	}



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


=head2 multiline
  Function : Specify feature types that can be split over multiple lines use 'all' to
  allow any feature be multiline
  Example  : $gff->multiline('cds');
=cut

	sub multiline {
		my $types = pop;
		$types =~ tr/[A-Z]/[a-z]/;
		my @types = split /\|/,$types;
	    while (my $type = shift @types){
	    	$multiline{$type} = 1;
	    }
	    return scalar keys %multiline;
	}


=head2 is_multiline
  Function : Test whether feature type can be split ove multiple lines
  Example  : $gff->is_multiline('cds');
=cut

	sub is_multiline {
		return 1 if $multiline{'all'};
		my $type = pop;
		$type =~ tr/[A-Z]/[a-z]/;
		return 1 if $multiline{$type};
		return;
	}


=head2 parse_file
  Function : Convenience method, calls parse_chunk with no separator to read an entire
             gff3 file into a tree structure.
  Example  : $gff->parse_file();
=cut

	sub parse_file {
		my $node = shift;
		$node->parse_chunk();
		return 1;
	}


=head2 parse_chunk
  Function : Reads a chunk of a gff3 file into a tree structure, handling multiline 
             features and ordering features on the same region.  if called with no 
             parameters the whole file will be parsed.
  Example  : $gff->parse_chunk('separator','###');
  Example  : $gff->parse_chunk('change','region');
  Example  : $gff->parse_chunk('separator','###');
=cut

	sub parse_chunk {
		my ($node,$split_by,$break_on) = @_;
		foreach my $it ($node->clear_daughters) { $it->delete_tree }
		if ($lastline){
			$lastline = undef;
			return;
		}
		#$ids{'root'} = $node;
		if ($split_by){
			return unless $break_on;
			if ($split_by eq 'change'){
				return unless $break_on eq 'region';
				# TODO: consider adding support for other changes
			}
		}
		my $fasta;
		my $region;
		my $seq = '';
		my $ctr = 0;
		while (<>){
		    $lastline = $_ if eof();
	    	chomp;
			if ($split_by && $split_by eq 'separator'){
				last if $_ =~ m/^$break_on/;
			}
			if (my $ret = is_comment($_)){
				if ($ret == -9){
					$fasta = substr $_,1;
				}
				else {
					$fasta = undef;
					$seq = '';
					$region = undef;
				}
				next;
			}
			elsif ($fasta){
				$seq .= $_;
				$region = $node->by_attributes(['_type','region'],['_seq_name',$fasta]) unless $region;
				unless ($region){
					$region = $node->make_region($fasta,'region');
				}
				$region->{attributes}->{_seq} = $seq;
				$region->{attributes}->{_end} = length $seq;
				next;
			}
			my $parent = $node;
			if ($has_comments){
				$_ = delete_comments($_);
			}
			my ($data,$attribs) = parse_gff_line($_,$separator);
			next unless $data;
			my %attributes;
			$attributes{'_seq_name'} = $data->[0];
			$attributes{'_source'} = $data->[1];
			$attributes{'_type'} = $type_map{$data->[2]} ? $type_map{$data->[2]} : $data->[2];
			$attributes{'_start'} = $data->[3];
			$attributes{'_end'} = $data->[4];
			$attributes{'_score'} = $data->[5];
			$attributes{'_strand'} = $data->[6];
			$attributes{'_phase'} = $data->[7];
			if ($attribs->{'Parent'} && $ids{$attribs->{'Parent'}}){
				$parent = $ids{$attribs->{'Parent'}};
				$ctr++;
			}
			else {
				# use the root node as an orphanage
				$parent = $node;
				#$attribs->{'Parent'} = 'root';
			}

			if (!$attribs->{'ID'}){ # need to do something about features that lack IDs
				my $behaviour = $node->lacks_id($attributes{'_type'});
				next if $behaviour eq 'ignore';
				if ($behaviour eq 'warn'){
					warn "WARNING: the feature on line $. does not have an ID, skipping feature\n";
					next;
				}
				elsif ($behaviour eq 'die'){
					die "ERROR: the feature on line $. does not have an ID, cannot continue\n";
				}
				elsif ($behaviour ne 'make'){
					if ($attribs->{$behaviour}){
						$attribs->{'ID'} = $attribs->{$behaviour};
					}
					else {
						$behaviour = 'make';
					}
				}
				if ($behaviour eq 'make'){
					if (is_multiline($attributes{'_type'})){
						# test whether parent has a child of this type with an ID already
						if ($attribs->{'Parent'} && $parents{$attribs->{'Parent'}}{$attributes{'_type'}}){
							$attribs->{'ID'} = $parents{$attribs->{'Parent'}}{$attributes{'_type'}};
						}
					}
					if (!$attribs->{'ID'}){
						my $prefix = $attributes{'_type'}.'___';
						my $suffix = 0;
						if ($suffices{$prefix}){
							$suffix = $suffices{$prefix};
							$suffices{$prefix}++;
						}
						while (by_id($prefix.$suffix)){
							$suffix++;
							$suffices{$prefix} = $suffix + 1;
						}
						$attribs->{'ID'} = $prefix.$suffix;
						if (is_multiline($attributes{'_type'}) && $attribs->{'Parent'}){
							$parents{$attribs->{'Parent'}}{$attributes{'_type'}} = $attribs->{'ID'};
						}
					}
				}
			}

			$attribs->{'ID'} =~ s/'//g;

			if (ref $attribs->{'Parent'} eq 'ARRAY'){ # multiparent feature
				my $base_id = $attribs->{'ID'};
				delete $attribs->{'ID'};
				my @parents = @{$attribs->{'Parent'}};
				delete $attribs->{'Parent'};
				for (my $p = 0; $p < @parents; $p++){
					$attributes{'Parent'} = $parents[$p];
					my $id;
					if ($p == 0){
						$attributes{'_duplicate'} = 0;
						$id = $base_id;
					}
					else {
						$attributes{'_duplicate'} = 1;
						$id = $base_id.'._'.$p;
					}
					$id = $p == 0 ? $base_id : $base_id.'._'.$p;
					$attributes{'ID'} = $id;
					if ($ids{$parents[$p]}){
						$ids{$id} = $ids{$parents[$p]}->new_daughter({%attributes,%$attribs});
					}
					else {
						$ids{$id} = $node->new_daughter({%attributes,%$attribs});
					}
					$ids{$id}->id($id);
					if ($attribs->{'Name'}){
						$ids{$id}->name($attribs->{'Name'});
					}
					else {
						$ids{$id}->name($id);
					}
					push @{$by_start{$attributes{'_seq_name'}}{$attributes{'_type'}}{$attributes{'_start'}}},$id;
				}
			}
			elsif (my $existing = $ids{$attribs->{'ID'}}){ # test if ID already exists, then treat as clash or multline
				if (is_multiline($attributes{'_type'}) &&
					$existing->{attributes}->{'_seq_name'} eq $attributes{'_seq_name'} &&
					$existing->{attributes}->{'_type'} eq $attributes{'_type'} &&
					$existing->{attributes}->{'_strand'} eq $attributes{'_strand'} &&
					(!$attribs->{'Parent'} ||
					($attribs->{'Parent'} &&
					$existing->{attributes}->{'Parent'} eq $attribs->{'Parent'}))
					){
					unless ($existing->{attributes}->{'_attributes'}){
						foreach my $attr (keys %{$existing->{attributes}}){
							$existing->{attributes}->{'_attributes'}->{$attr} = 1;
							push @{$existing->{attributes}->{$attr.'_array'}},$existing->{attributes}->{$attr};
						}
					}
					# change _start and _end, storing individual component values in an array
					if ($attributes{'_start'} > $existing->{attributes}->{'_start'}){
						# don't use unshift but choose between push and splice according to the existing values
						if ($attributes{'_end'} > $existing->{attributes}->{'_end'}){
							foreach my $attr (keys %{$existing->{attributes}->{'_attributes'}}){
								push @{$existing->{attributes}->{$attr.'_array'}},defined $attributes{$attr} ? $attributes{$attr} : defined $attribs->{$attr} ? $attribs->{$attr} : undef;
							}
							foreach my $attr (keys %$attribs){
								unless ($existing->{attributes}->{'_attributes'}->{$attr}){
									for (my $i = 0; $i - 1 < @{$existing->{attributes}->{'_start_array'}}; $i++){
										push @{$existing->{attributes}->{$attr.'_array'}},undef;
									}
									push @{$existing->{attributes}->{$attr.'_array'}},$attribs -> {$attr};
									$existing->{attributes}->{'_attributes'}->{$attr} = 1;
								}
							}
							$existing->{attributes}->{'_end'} = $attributes{'_end'};
						}
						else {
							# find the position in the array to insert the current values
							my $index = 0;
							while ($index < @{$existing->{attributes}->{'_start_array'}} && $attributes{'_start'} > $existing->{attributes}->{'_start_array'}[$index]){
								$index++;
							}
							foreach my $attr (keys %{$existing->{attributes}->{'_attributes'}}){
								splice @{$existing->{attributes}->{$attr.'_array'}},$index,0,defined $attributes{$attr} ? $attributes{$attr} : defined $attribs->{$attr} ? $attribs->{$attr} : undef;
							}
							foreach my $attr (keys %$attribs){
								unless ($existing->{attributes}->{'_attributes'}->{$attr}){
									for (my $i = 0; $i - 1 < @{$existing->{attributes}->{'_start_array'}}; $i++){
										push @{$existing->{attributes}->{$attr.'_array'}},undef;
									}
									splice @{$existing->{attributes}->{$attr.'_array'}},$index,0,$attribs -> {$attr};
									$existing->{attributes}->{'_attributes'}->{$attr} = 1;
								}
							}
						}

					}
					else {
						foreach my $attr (keys %{$existing->{attributes}->{'_attributes'}}){
							unshift @{$existing->{attributes}->{$attr.'_array'}},defined $attributes{$attr} ? $attributes{$attr} : defined $attribs->{$attr} ? $attribs->{$attr} : undef;
						}
						for (my $s = 0; $s < @{$by_start{$attributes{'_seq_name'}}{$attributes{'_type'}}{$existing->{attributes}->{'_start'}}}; $s++){
							my $possible = $by_start{$attributes{'_seq_name'}}{$attributes{'_type'}}{$existing->{attributes}->{'_start'}}[$s];
							if ($possible eq $existing){
								# remove and replace
								splice(@{$by_start{$attributes{'_seq_name'}}{$attributes{'_type'}}{$existing->{attributes}->{'_start'}}}, $s, 1);
								$existing->{attributes}->{'_start'} = $attributes{'_start'};
								push @{$by_start{$attributes{'_seq_name'}}{$attributes{'_type'}}{$attributes{'_start'}}},$existing;
								last;
							}
						}
						foreach my $attr (keys %$attribs){
							unless ($existing->{attributes}->{'_attributes'}->{$attr}){
								for (my $i = 1; $i < @{$existing->{attributes}->{'_start_array'}}; $i++){
									unshift @{$existing->{attributes}->{$attr.'_array'}},undef;
								}
								unshift @{$existing->{attributes}->{$attr.'_array'}},$attribs -> {$attr};
								$existing->{attributes}->{'_attributes'}->{$attr} = 1;
							}
						}
						$existing->{attributes}->{'_start'} = $attributes{'_start'};
					}
				}
				else {
					# ID clash
					#print $parent->id(),"\n";
					#print $attributes{'_type'},"\n";
					#print is_multiline($attributes{'_type'}),"\n";
					die "ERROR: feature ID $attribs->{'ID'} has already been used (line $.)\nERROR: should you call multiline('".$attributes{'_type'}."') before parse_file()?\n";
				}
			}
			else {
				#print "$attribs->{'ID'}\n";
				$ids{$attribs->{'ID'}} = $parent->new_daughter({%attributes,%$attribs});
				$ids{$attribs->{'ID'}}->id($attribs->{'ID'});
				if ($attribs->{'Name'}){
					$ids{$attribs->{'ID'}}->name($attribs->{'Name'});
				}
				else {
					$ids{$attribs->{'ID'}}->name($attribs->{'ID'});
				}
				push @{$by_start{$attributes{'_seq_name'}}{$attributes{'_type'}}{$attributes{'_start'}}},$ids{$attribs->{'ID'}};
			}
			if ($split_by && $split_by eq 'change'){
				my $nextline = <>;
				if (defined($nextline)){	
					next if is_comment($nextline);
					IO::Unread::unread ARGV, $nextline;
					next if $fasta;
					if ($has_comments){
						$nextline = delete_comments($nextline);
					}
					my ($nextdata,undef) = parse_gff_line($nextline,$separator);
					last if $nextdata && $nextdata->[0] ne $data->[0];
				}
				else {
					last;
				}
			}
			
			
		}

		# loop through orphanage to see if anything can be done with unparented features
		my $orphans = 0;
		my @orphans = $node->daughters();
		while (scalar @orphans != $orphans){
			$orphans = scalar @orphans;
			for (my $o = 0; $o < @orphans; $o++){
				if ($orphans[$o]->{attributes}->{'Parent'} && $ids{$orphans[$o]->{attributes}->{'Parent'}}){
					# move the orphan node to a new parent
					$orphans[$o]->unlink_from_mother();
					$ids{$orphans[$o]->{attributes}->{'Parent'}}->add_daughter($orphans[$o]);
				}
			}
			@orphans = $node->daughters();
		}
		for (my $o = 0; $o < @orphans; $o++){
			if ($orphans[$o]->{attributes}->{'Parent'}){
				my $behaviour = $node->undefined_parent();
				if ($behaviour eq 'die'){
					die "ERROR: no parent feature exists for ",$orphans[$o]->{attributes}->{_type}," ",$orphans[$o]->id()," with the ID $orphans[$o]->{attributes}->{'Parent'}, cannot continue\n";
				}
				else {
					# wait for validation to take care of things - needs improving
				}
			}
		}
		return 1;
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
		my $suffix = 0;
		$suffix = 0 unless $suffix;
		if ($suffices{$prefix}){
			$suffix = $suffices{$prefix};
			$suffices{$prefix}++;
		}
		while (by_id($prefix.$suffix)){
			$suffix++;
			$suffices{$prefix} = $suffix + 1;
		}
		my $id = $prefix.$suffix;
		$node->id($id);
		$node->{attributes}->{'ID'} = $id;
		if ($node->{attributes}->{'Name'}){
			$node->name($node->{attributes}->{'Name'});
		}
		else {
			$node->name($id);
		}
		$ids{$id} = $node;
		push @{$by_start{$node->{attributes}->{_seq_name}}{$node->{attributes}->{_type}}{$node->{attributes}->{_start}}},$node;
		return $id;
	}


=head2 by_start
  Function : Fetch an arrayref of nodes start position. If multiple types are given, and
             more than one type has a match, the first specified type will be returned
             preferentially.
  Example  : $node_arrayref = by_start('scaf1','exon',432);
             $node_arrayref = by_start('scaf1','mrna|exon',432);
=cut

	sub by_start  {
		my $start = pop;
		my $type = pop;
		my $seq_name = pop;
		if ($type =~ m/\|/){
			my @types = split /\|/,$type;
			while ($type = shift @types){
				return $by_start{$seq_name}{$type}{$start} if $by_start{$seq_name}{$type};
			}
		}
		return $by_start{$seq_name}{$type}{$start};
	}


=head2 nearest_start
  Function : Fetch an arrayref of nodes as close as possible to the start position. If
             multiple types are given, and more than one type has a closest match, the
             first specified type will be returned preferentially.
  Example  : $node_arrayref = nearest_start('scaf1','exon',432);
=cut

	sub nearest_start  {
		my $start = pop;
		my $type = pop;
		my $seq_name = pop;
		my $prev_begin = 0;
		my @types;
		if ($type =~ m/\|/){
			@types = split /\|/,$type;
		}
		else {
			push @types, $type;
		}
		while (my $mtype = shift @types){
			next unless $by_start{$seq_name}{$mtype};
			foreach my $begin (sort { $a <=> $b } keys %{$by_start{$seq_name}{$mtype}}){
				next if $begin < $prev_begin;
				last if $begin > $start;
				$prev_begin = $begin;
				$type = $mtype;
			}
		}
		return $by_start{$seq_name}{$type}{$prev_begin};
	}


}

# use nesting to allow sub to retain access to private variables

{
	my %lacks;

=head2 lacks_id
  Function : get/set behaviour for features that lack an id
  Example  : $gff->lacks_id('region','ignore'); # default
  Example  : $gff->lacks_id('exon','warn');
  Example  : $gff->lacks_id('gene','die');
  Example  : $gff->lacks_id('all','make');
  Example  : $behaviour = $gff->lacks_id('gene');
=cut

	sub lacks_id  {
		my $node = shift;
		my $type = shift || 'all';
		my $behaviour = shift;
		$type =~ tr/[a-z]/[A-Z]/;
		$lacks{$node}{$type} = $behaviour if $behaviour;
		return $lacks{$node}{$type} ? $lacks{$node}{$type} : $lacks{$node}{'all'} ? $lacks{$node}{'all'} : 'ignore';
	}
}

=head2 undefined_parent
  Function : get/set behaviour for features that should have a parent but don't
  Example  : $gff->undefined_parent('die');
  Example  : $gff->undefined_parent('make'); # default
  Example  : $behaviour = $gff->undefined_parent();
=cut

sub undefined_parent  {
	my $node = shift;
	my $behaviour = shift;
	$node->{_attributes}->{'_undefined_parent'} = $behaviour if $behaviour;
	return $node->{_attributes}->{'_undefined_parent'} ? $node->{_attributes}->{'_undefined_parent'} : 'make';
}


# use nesting to allow sub to retain access to private variables

{
	my %features;
	my %parents;

=head2 by_type
  Function : returns an ordered array of features of a given type
  Example  : @nodes = $gff->by_type('exon');
=cut

	sub by_type {
    	my ($self, $type) = @_;
    	$self->order_features($type);
    	return @{$features{$type}}[ -@{$features{$type}} .. -2 ];;
	}

=head2 next_feature
  Function : Sequentially fetch daughter features of a node by type
  Example  : $gene->next_feature('exon');
=cut

	sub next_feature {
		my ($self, $type) = @_;
		unless ($features{$type} && @{$features{$type}} && $parents{$type} && $parents{$type} eq $self->id()){
			$self->order_features($type);
		}
		return shift @{$features{$type}};
	}

=head2 order_features
  Function : order daughter features of a given type sequentially
  Example  : $gene->order_feature('exon');
=cut

	sub order_features {
		my ($self, $type, $strand) = @_;
		my @unsorted = by_attribute($self,'_type',$type);
		@{$features{$type}} = ($strand && $strand eq '-') ? sort { $b->{attributes}->{_start} <=> $a->{attributes}->{_start} } @unsorted : sort { $a->{attributes}->{_start} <=> $b->{attributes}->{_start} } @unsorted;
		push @{$features{$type}},0;
		$parents{$type} = $self->id;
		return (scalar(@{$features{$type}})-1);
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
		$attributes{'Parent'} = $self->id;
		if ($location eq 'internal'){
			if (@{$features{$type}} > 2){
				for (my $i = 1; $i < @{$features{$type}} - 1; $i++){
					$attributes{'_start'} = $features{$type}[($i-1)]->{attributes}->{_end} + 1;
					$attributes{'_end'} = $features{$type}[$i]->{attributes}->{_start} - 1;
					next if $attributes{'_end'} <= $attributes{'_start'};
					if (my $feature = $self->find_daughter(\%attributes)){
						$self->add_daughter($feature);
					}
					else {
						my $node = $self->new_daughter(\%attributes);
						$node->make_id($new_type);
					}
				}
			}
		}
		if ($features{$type} > 1){
			if ($location eq 'before' || $location eq 'external'){
				if ($features{$type}[0]->{attributes}->{_end} < $self->{attributes}->{_end}){
					if ($features{$type}[0]->{attributes}->{_strand} eq '-'){
						$attributes{'_start'} = $features{$type}[0]->{attributes}->{_end} + 1;
						$attributes{'_end'} = $self->{attributes}->{_end};
					}
					else {
						$attributes{'_start'} = $self->{attributes}->{_start};
						$attributes{'_end'} = $features{$type}[0]->{attributes}->{_start} - 1;
					}
					return if $attributes{'_end'} <= $attributes{'_start'};
					if (my $feature = $self->find_daughter(\%attributes)){
						$self->add_daughter($feature);
					}
					else {
						my $node = $self->new_daughter(\%attributes);
						$node->make_id($new_type);
					}
				}
			}
			if ($location eq 'after' || $location eq 'external'){
				if ($features{$type}[-2]->{attributes}->{_end} < $self->{attributes}->{_end}){
					if ($features{$type}[0]->{attributes}->{_strand} eq '-'){
						$attributes{'_end'} = $features{$type}[-2]->{attributes}->{_start} - 1;
						$attributes{'_start'} = $self->{attributes}->{_start};
					}
					else {
						$attributes{'_end'} = $self->{attributes}->{_end};
						$attributes{'_start'} = $features{$type}[-2]->{attributes}->{_end} + 1;
					}
					return if $attributes{'_end'} <= $attributes{'_start'};
					if (my $feature = $self->find_daughter(\%attributes)){
						$self->add_daughter($feature);
					}
					else {
						my $node = $self->new_daughter(\%attributes);
						$node->make_id($new_type);
					}
				}
			}
		}
		return;
	}

=head2 find_daughter
  Function : Find out whether an element already has a daughter with a given set of
             attributes
  Example  : $gene->find_daughter(\%attributes);
=cut

	sub find_daughter {
		my $self = shift;
		my $attributes = shift;
		my @possibles = by_start($attributes->{'_seq_name'},$attributes->{'_type'},$attributes->{'_start'});
		while (my $feature = shift @possibles){
			if ($attributes->{'_end'} == $feature->[0]->{attributes}->{'_end'}){
				return $feature->[0];
			}
		}
		return;
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
					'skip' => \&validation_skip,
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
		# relations: hasParent, hasChild, hasSister, >, gt, <, lt, ==, eq, >=, <=, !=, ne
		# mrna hasParent gene
		# mrna|exon <[start,end] SELF
		# mrna <=[end] PARENT warn
		# exon hasParent mrna|transcript|gene
		# cds hasSister exon

		my ($self,$type,$relation,$alt_type,$flag) = @_;
		$type =~ tr/[A-Z]/[a-z]/;
		my @type = split /\|/,$type;
		for (my $t = 0; $t < @type; $t++){
			push @{$expectations{$type[$t]}},{'relation' => $relation, 'alt_type' => $alt_type, 'flag' => $flag} || return;
		}
		return scalar keys %expectations;
	}


=head2 find_sister
  Function : find sister features accounting for multiline, nested matching and encompassing
             features
  Example  : $sister = $cds->find_sister('exon');
           : $sister = $exon->find_sister('cds');
=cut

	sub find_sister {
		my $self = shift;
		my $alt_type = shift;
		my $scale = shift;
		my $parent = $self->mother();
		my ($sister,$size);
		if (!is_multiline($self->{attributes}->{_type}) && !is_multiline($alt_type) ||
			is_multiline($self->{attributes}->{_type}) && is_multiline($alt_type)){

			while (my $feature = $parent->next_feature($alt_type)){
				if ($self->{attributes}->{_start} == $feature->{attributes}->{_start} && $self->{attributes}->{_end} == $feature->{attributes}->{_end}){
					# twinSister
					$sister = $feature;
					$size = 'twin';
					last; # match found, don't check any more features
				}
				elsif ($self->{attributes}->{_start} >= $feature->{attributes}->{_start} && $self->{attributes}->{_end} >= $feature->{attributes}->{_end}){
					# littleSister
					$sister = $feature;
					$size = 'little';
					# continue the loop to find a better match (i.e. twin sister)
				}

				elsif ($self->{attributes}->{_start} >= $feature->{attributes}->{_start} && $self->{attributes}->{_end} <= $feature->{attributes}->{_end}){
					# bigSister
					$sister = $feature;
					$size = 'big';
					# continue the loop to find a better match (i.e. twin sister)
				}
			}
		}
		elsif (is_multiline($self->{attributes}->{_type}) && !is_multiline($alt_type)){
			my @starts = $self->{attributes}->{_start_array} ? @{$self->{attributes}->{_start_array}} : ($self->{attributes}->{_start});
			my @ends = $self->{attributes}->{_end_array} ? @{$self->{attributes}->{_end_array}} : ($self->{attributes}->{_end});
			for (my $i = 0; $i < @starts; $i++){
				my @features = $parent->by_type($alt_type);
				while (my $feature = shift @features){
					$sister = undef;
					if ($starts[$i] == $feature->{attributes}->{_start} && $ends[$i] == $feature->{attributes}->{_end}){
						# twinSister
						$sister = $feature;
						$size = 'twin';
					}
					elsif ($starts[$i] <= $feature->{attributes}->{_start} && $ends[$i] >= $feature->{attributes}->{_end}){
						# littleSister
						$sister = $feature;
						$size = 'little';
					}

					elsif ($starts[$i] >= $feature->{attributes}->{_start} && $ends[$i] <= $feature->{attributes}->{_end}){
						# bigSister
						$sister = $feature;
						$size = 'big';
					}
					last if $sister; # # match found for this part of the multiline feature, don't check any more features
				}
				last unless $sister; # all parts of the multiline feature must have a match
			}
		}
		else { # !is_multiline($self->{attributes}->{_type}) and is_multiline($alt_type)
			my @features = $parent->by_type($alt_type);
			while (my $feature = shift @features){
				my @starts = $feature->{attributes}->{_start_array} ? @{$feature->{attributes}->{_start_array}} : ($feature->{attributes}->{_start});
				my @ends = $feature->{attributes}->{_end_array} ? @{$feature->{attributes}->{_end_array}} : ($feature->{attributes}->{_end});
				for (my $i = 0; $i < @starts; $i++){
					$sister = undef;
					if ($starts[$i] == $self->{attributes}->{_start} && $ends[$i] == $self->{attributes}->{_end}){
						# twinSister
						$sister = $feature;
						$size = 'twin';
					}
					elsif ($starts[$i] <= $self->{attributes}->{_start} && $ends[$i] >= $self->{attributes}->{_end}){
						# littleSister
						$sister = $feature;
						$size = 'big';
					}

					elsif ($starts[$i] >= $self->{attributes}->{_start} && $ends[$i] <= $self->{attributes}->{_end}){
						# bigSister
						$sister = $feature;
						$size = 'little';
					}
					last if $sister; # match found for this part of the multiline feature, don't check any more features
				}
				last unless $sister; # all parts of the multiline feature must have a match
			}
		}
		return $sister;
	}


=head2 make_sister
  Function : make sister features accounting for multiline and single line interactions
  Example  : $sister = $cds->make_sister('exon');
           : $sister = $exon->make_sister('cds');
=cut

	sub make_sister {
		my $self = shift;
		my $alt_type = shift;
		my $parent = $self->mother();
		my $sister;
		my @attributes = ('_seq_name','_source','_start','_end','_score','_strand','_phase','Parent');
		if (is_multiline($self->{attributes}->{_type}) && is_multiline($alt_type) ||
			!is_multiline($self->{attributes}->{_type}) && !is_multiline($alt_type)){
			$sister = $self->copy({no_attribute_copy => 1});
			foreach my $attribute (@attributes){
				$sister->{attributes}->{$attribute} = $self->{attributes}->{$attribute};
				if ($self->{attributes}->{$attribute.'_array'}){
					$sister->{attributes}->{$attribute.'_array'} = $self->{attributes}->{$attribute.'_array'};
				}
			}
			$parent->add_daughter($sister);
			$sister->{attributes}->{_type} = $alt_type;
			$sister->make_id($alt_type);
			$sister->{attributes}->{Name} = $sister->name();
		}
		elsif (is_multiline($self->{attributes}->{_type}) && !is_multiline($alt_type)){
			my @starts = $self->{attributes}->{_start_array} ? @{$self->{attributes}->{_start_array}} : ($self->{attributes}->{_start});
			my @ends = $self->{attributes}->{_end_array} ? @{$self->{attributes}->{_end_array}} : ($self->{attributes}->{_end});
			for (my $i = 0; $i < @starts; $i++){
				my @features = $parent->by_type($alt_type);
				while (my $feature = shift @features){
					$sister = undef;
					if ($starts[$i] == $feature->{attributes}->{_start} && $ends[$i] == $feature->{attributes}->{_end}){
						# twinSister
						$sister = $feature;
					}
					elsif ($starts[$i] <= $feature->{attributes}->{_start} && $ends[$i] >= $feature->{attributes}->{_end}){
						# littleSister
						$sister = $feature;
					}

					elsif ($starts[$i] >= $feature->{attributes}->{_start} && $ends[$i] <= $feature->{attributes}->{_end}){
						# bigSister
						$sister = $feature;
					}
					last if $sister; # # match found for this part of the multiline feature, don't check any more features
				}
				if (!$sister){ # make a new sister as all parts of the multiline feature must have a match
					$sister = $self->copy({no_attribute_copy => 1});
					foreach my $attribute (@attributes){
						$sister->{attributes}->{$attribute} = $self->{attributes}->{$attribute.'_array'} ? $self->{attributes}->{$attribute.'_array'}[$i] : $self->{attributes}->{$attribute};
					}
					$parent->add_daughter($sister);
					$sister->{attributes}->{_type} = $alt_type;
					$sister->make_id($alt_type);
					$sister->{attributes}->{Name} = $sister->name();
				}
			}
		}
		else { # !is_multiline($self->{attributes}->{_type}) and is_multiline($alt_type)
			my @features = $parent->by_type($alt_type);
			while (my $feature = shift @features){
				my @starts = $feature->{attributes}->{_start_array} ? @{$feature->{attributes}->{_start_array}} : ($feature->{attributes}->{_start});
				my @ends = $feature->{attributes}->{_end_array} ? @{$feature->{attributes}->{_end_array}} : ($feature->{attributes}->{_end});
				for (my $i = 0; $i < @starts; $i++){
					$sister = undef;
					if ($starts[$i] == $self->{attributes}->{_start} && $ends[$i] == $self->{attributes}->{_end}){
						# twinSister
						$sister = $feature;
					}
					elsif ($starts[$i] <= $self->{attributes}->{_start} && $ends[$i] >= $self->{attributes}->{_end}){
						# littleSister
						$sister = $feature;
					}

					elsif ($starts[$i] >= $self->{attributes}->{_start} && $ends[$i] <= $self->{attributes}->{_end}){
						# bigSister
						$sister = $feature;
					}
					last if $sister; # match found for this part of the multiline feature, don't check any more features
				}
				if (!$sister){ # make a new sister as all parts of the multiline feature must have a match
					####################################################################################
					# TODO: make a new/append to an existing multiline feature.                        #
					#       make new is relatively straightforward but appending gets more complicated #
					####################################################################################
					die "ERROR: Making a new multiline sister feature from a non-multiline feature is not yet supported\n       Cannot run \$".$self->{attributes}->{_type}."->make_sister with multiline('$alt_type') and not multiline('".$self->{attributes}->{_type}."')\n       Consider removing add_expectation('$alt_type','hasSister','".$self->{attributes}->{_type}."','make')\n";
				}
			}
		}
		return $sister;
	}


=head2 make_child
  Function : make child features
  Example  : $child = $mrna->make_child('exon');
=cut

	sub make_child {
		my $self = shift;
		my $alt_type = shift;
		my $child;
		my @attributes = ('_seq_name','_source','_start','_end','_score','_strand','_phase');
		$child = $self->copy({no_attribute_copy => 1});
		foreach my $attribute (@attributes){
			$child->{attributes}->{$attribute} = $self->{attributes}->{$attribute};
			if ($self->{attributes}->{$attribute.'_array'}){
				$child->{attributes}->{$attribute.'_array'} = $self->{attributes}->{$attribute.'_array'};
			}
		}
		$self->add_daughter($child);
		$child->{attributes}->{_type} = $alt_type;
		$child->{attributes}->{Parent} = $self->id();
		$child->make_id($alt_type);
		$child->{attributes}->{Name} = $child->name();
		return $child;
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
					my $message = $type." ".$self->id.' does not have a parent of type '.$hashref->{'alt_type'};
					$actions{$hashref->{'flag'}}->($self,$message,$hashref) unless $self->mother->{attributes}->{_type} && $self->mother->{attributes}->{_type} =~ m/$hashref->{'alt_type'}/i;
				}
				elsif ($hashref->{'relation'} eq 'hasChild'){
					my $message = $type." ".$self->id.' does not have a child of type '.$hashref->{'alt_type'};
					$actions{$hashref->{'flag'}}->($self,$message,$hashref) unless $self->by_type($hashref->{'alt_type'});
				}
				elsif ($hashref->{'relation'} eq 'hasSister'){
					my $message = $type." ".$self->id.' does not have a sister of type '.$hashref->{'alt_type'};
					$actions{$hashref->{'flag'}}->($self,$message,$hashref) unless $self->find_sister($hashref->{'alt_type'});
				}
				else {
					my @relation = split /[\[\]]/,$hashref->{'relation'};
					my @attrib = split /,/,$relation[1];
					my $message = $type.' '.$self->id.'->('.$attrib[0].') is not '.$relation[0].' '.$hashref->{'alt_type'}.'->('.$attrib[-1].')';
					$message .=  '('.$self->mother->id.')' if $self->mother->id;
					my $first = $self->{attributes}->{$attrib[0]};
					my $second = $hashref->{'alt_type'} =~ m/self/i ? $self->{attributes}->{$attrib[-1]} : $self->mother->{attributes}->{$attrib[-1]};
					$actions{$hashref->{'flag'}}->($self,$message,$hashref) unless compare($first,$second,$relation[0]);
				}
			}
		}
		return;
	}

=head2 validate_all
  Function : run validate on all features, optionally of type specified in parentheses
  Example  : $gff->validate_all();
  Example  : $gff->validate_all('gene');
  Example  : $gff->validate_all('exon');
=cut

	sub validate_all {
		my $self = shift;
		my $type = shift;
		my @features;
		if ($type){
			@features = $self->by_type($type);
		}
		else {
			@features = $self->descendants();
		}
		while (my $feature = shift @features){
			$feature->validate();
		}
		return 1;
	}


	sub validation_ignore {
		# nothing happens
		return;
	}

	sub validation_warning {
		my $self = shift;
		my $message = shift;
		warn "WARNING: $message\n";
		return;
	}

	sub validation_die {
		my $self = shift;
		my $message = shift;
		die "ERROR: $message\n";
	}

	sub validation_skip {
		my $self = shift;
		my $message = shift;
		warn "WARNING: $message\n";
		$self->{attributes}->{_skip} = 'true';
	}


=head2 validation_find
  Function : find a feature to satisfy an expectation - limited functionality at present
  Example  : validation_find($expectation_hashref);
=cut

	sub validation_find {
		# TODO 	- handle relationships other than parent
		my $self = shift;
		my $message = shift;
		my $expectation = shift;
		my %attributes;
		$attributes{'_seq_name'} = $self->{attributes}->{_seq_name};
		$attributes{'_source'} = 'GFFTree';
		$attributes{'_type'} = $expectation->{'alt_type'};
		$attributes{'_score'} = '.';
		$attributes{'_strand'} = $self->{attributes}->{_strand};
		$attributes{'_phase'} = '.';

		my $relative;
		if ($expectation->{'relation'} eq 'hasParent'){
			my @possibles = by_start($self->{attributes}->{'_seq_name'},$expectation->{'alt_type'},$self->{attributes}->{'_start'});
			while ($relative = shift @possibles){
				if ($self->{attributes}->{'_end'} == $relative->[0]->{attributes}->{'_end'}){
					last;
				}
			}
			unless ($relative){
				@possibles = nearest_start($self->{attributes}->{'_seq_name'},$expectation->{'alt_type'},$self->{attributes}->{'_start'}) unless @possibles;
				while ($relative = shift @possibles){
					if ($self->{attributes}->{'_end'} <= $relative->[0]->{attributes}->{'_end'}){
						last;
					}
				}
			}
			if ($relative){
				my $parent = $relative->[0];
				$attributes{'_start'} = $parent->{attributes}->{_start};
				$attributes{'_end'} = $parent->{attributes}->{_end};
				$self->{attributes}->{Parent} = $parent->id();
				$self->unlink_from_mother();
				$parent->add_daughter($self);
				return $relative->[0];
			}

		}
		return;
	}

=head2 validation_make
  Function : make a feature to satisfy an expectation - limited to parents at the moment
  Example  : validation_make($expectation_hashref);
=cut

	sub validation_make {
		my $self = shift;
		my $expectation = pop;
		my %attributes;
		$attributes{'_seq_name'} = $self->{attributes}->{_seq_name};
		$attributes{'_source'} = 'GFFTree';
		$attributes{'_type'} = $expectation->{'alt_type'};
		$attributes{'_score'} = '.';
		$attributes{'_strand'} = $self->{attributes}->{_strand};
		$attributes{'_phase'} = '.';
		if ($expectation->{'relation'} eq 'hasParent'){
			if ($expectation->{'alt_type'} eq 'region'){
				# find limits of the region
				my @features = by_attribute($self,'_seq_name',$self->{attributes}->{_seq_name});
				# assuming regions always start at 1, comment out the code below if not true
				# my @sorted = sort { $a->{attributes}->{_start} <=> $b->{attributes}->{_start} } @features;
				$attributes{'_start'} = 1; # $sorted[0]->{attributes}->{_start};
				my @sorted = sort { $b->{attributes}->{_end} <=> $a->{attributes}->{_end} } @features;
				# TODO - if a sequence is available, use that to find the length of a region
				$attributes{'_end'} = $sorted[0]->{attributes}->{_end};
				$attributes{'_strand'} = '+';
			}
			else {
				$attributes{'_start'} = $self->{attributes}->{_start};
				$attributes{'_end'} = $self->{attributes}->{_end};
			}
			$attributes{'Parent'} = $self->{attributes}->{Parent} if $self->{attributes}->{Parent};
			my $node = $self->mother->new_daughter(\%attributes);
			$node->make_id($expectation->{'alt_type'});
			$self->{attributes}->{Parent} = $node->id();
			$self->unlink_from_mother();
			$node->add_daughter($self);
			return $node;
		}
		elsif ($expectation->{'relation'} eq 'hasSister'){
			my $sister = $self->make_sister($expectation->{'alt_type'});
			return $sister;
		}
		elsif ($expectation->{'relation'} eq 'hasChild'){
			my $child = $self->make_child($expectation->{'alt_type'});
			return $child;
		}
		return;
	}

=head2 validation_force
  Function : find a feature to satisfy an expectation if possible, otherwise make one
  Example  : validation_force($expectation_hashref);
=cut

	sub validation_force {
		my $self = shift;
		my $expectation = pop;
		my $relative = $self->validation_find('',$expectation);
		return $relative if $relative;
		$relative = $self->validation_make('',$expectation);
		return $relative;
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
	return;
}


=head2 make_region
  Function : make a feature required during parsing, maybe combine with validation_make?
  Example  : parse_make($expectation_hashref);
=cut

sub make_region {
	my $self = shift;
	my $name = shift;
	my $type = shift;
	my %attributes;
	$attributes{'_seq_name'} = $name;
	$attributes{'_source'} = 'GFFTree';
	$attributes{'_type'} = $type;
	$attributes{'_score'} = '.';
	$attributes{'_strand'} = '+';
	$attributes{'_phase'} = '.';
	$attributes{'_start'} = 1;
	my $node = $self->new_daughter(\%attributes);
	$node->make_id($type);
	return $node;
}


=head2 as_string
  Function : returns a gff representation of a single feature
  Example  : $gene->as_string();
  Example  : $gene->as_string(1); # skip features labeled as duplicates
=cut

sub as_string {
	my $self = shift;
	my $skip_dups = shift;
	my $line = '';
	my @col_nine;
	if (is_multiline($self->{attributes}->{_type}) && $self->{attributes}->{_start_array}){
		for (my $s = 0; $s < @{$self->{attributes}->{_start_array}}; $s++){
			foreach my $key (sort keys %{$self->{attributes}}){
				next if $key =~ m/^_/;
				next if $key =~ m/_array$/;
				my $attr = $self->{attributes}->{$key};
				if ($self->{attributes}->{$key.'_array'}){
					$attr = $self->{attributes}->{$key.'_array'}[$s];
				}
				if (ref $attr eq 'ARRAY') {
					my $value = join(',',@{$attr});
					$value =~ s/=/\%3D/g;
  					$value =~ s/;/\%3B/g;
  					$col_nine[$s] .= $key.'='.$value.';';
				}
				else {
					my $value = $attr;
					$value =~ s/\._\d+$//;
					$value =~ s/=/\%3D/g;
  					$value =~ s/;/\%3B/g;
					$col_nine[$s] .= $key.'='.$value.';';
				}
			}
			chop $col_nine[$s];
		}
	}
	else {
		foreach my $key (sort keys %{$self->{attributes}}){
			next if $key =~ m/^_/;
			if (ref $self->{attributes}->{$key} eq 'ARRAY') {
				my $value = join(',',@{$self->{attributes}->{$key}});
				$value =~ s/=/\%3D/g;
  				$value =~ s/;/\%3B/g;
	  			$col_nine[0] .= $key.'='.$value.';';
			}
			else {
				my $value = $self->{attributes}->{$key};
				$value =~ s/\._\d+$//;
				$value =~ s/=/\%3D/g;
  				$value =~ s/;/\%3B/g;
				$col_nine[0] .= $key.'='.$value.';';
			}
		}
		chop $col_nine[0];
	}
	my @start_array = ($self->{attributes}->{_start});
	my @end_array = ($self->{attributes}->{_end});
	my @phase_array = ($self->{attributes}->{_phase});
	my @score_array = ($self->{attributes}->{_score});
	if (is_multiline($self->{attributes}->{_type}) && $self->{attributes}->{_start_array}){
		for (my $s = 0; $s < @{$self->{attributes}->{_start_array}}; $s++){
			$start_array[$s] = $self->{attributes}->{_start_array}[$s];
			$end_array[$s] = $self->{attributes}->{_end_array}[$s];
			$phase_array[$s] = $self->{attributes}->{_phase_array}[$s];
			$score_array[$s] = $self->{attributes}->{_score_array}[$s];
		}
	}
	for (my $s = 0; $s < @start_array; $s++){
		next if $skip_dups && $self->{attributes}->{_duplicate};
		$line .= $self->{attributes}->{_seq_name}."\t";
		$line .= $self->{attributes}->{_source}."\t";
		$line .= $self->{attributes}->{_type}."\t";
		$line .= $start_array[$s]."\t";
		$line .= $end_array[$s]."\t";
		$line .= $score_array[$s]."\t";
		$line .= $self->{attributes}->{_strand}."\t";
		$line .= $phase_array[$s]."\t";
		$line .= $col_nine[$s]."\n";
	}
	return $line;
}

=head2 structured_output
  Function : returns a gff representation of a feature and all descendants
  Example  : $gene->structured_output();
=cut

sub structured_output {
	my $self = shift;
	return if $self->{attributes}->{_skip};
	my $output;
	$output .= $self->as_string(1);
	my @daughters = $self->daughters();
	while (my $daughter = shift @daughters){
		my $out = $daughter->structured_output();
		$output .= $out if $out;
		return if $daughter->{attributes}->{_skip};
	}
	return $output;
}



=head2 is_comment
  Function : returns true if a line is a comment
  Example  : is_comment($line);
=cut

sub is_comment {
	return -1 if $_ =~ m/^$/;
	return -9 if $_ =~ m/^>/;
	return length($1) if $_ =~ m/^(#+)/;
	return;
}

=head2 parse_gff_line
  Function : splits a line of gff into 8 fields and a key-value hash, escaping encoded
             characters
  Example  : parse_gff_line($line);
=cut

{
	my $col_count;
	my $col_count_flag = 'ignore';

=head2 expect_columns
  Function : expected number of columns in the input file with flag to ignore warn, die 
             or skip if the number of columns found does not match the expected number.  
  Example  : expect_columns(9,'skip');
=cut

	sub expect_columns {
		my $node = shift;
		$col_count = shift;
		my $flag = shift;
		$col_count_flag = $flag if $flag;
		return $col_count;
	}

=head2 parse_gff_line
  Function : splits a line of gff into 8 fields and a key-value hash, escaping encoded
             characters and building arrays of comma-separated values
  Example  : parse_gff_line($line);
=cut

	sub parse_gff_line {
		my ($line,$sep) = @_;
		my @data = split /$sep/,$line;
		if ($col_count){
			if ($col_count != @data){
				die "ERROR: Expected $col_count columns but found ".scalar(@data) if $col_count_flag eq 'die';
				warn "WARNING: Expected $col_count columns but found ".scalar(@data) if $col_count_flag eq 'warn';
				if ($col_count_flag eq 'skip'){
					#warn "WARNING: Expected $col_count columns but found ".scalar(@data).", skipping line";
					return;
				}
			}
		}
		chomp $data[8];
		my %attribs = split /[=;]/,$data[8];
		pop @data;
		foreach my $key (keys %attribs){
			# treat differently if commas are present
			if (!$attribs{$key}){
				delete $attribs{$key};
				next;
			}
			my @parts = split /,/,$attribs{$key};
			if ($parts[1]){
				$attribs{$key} = [];
				while (my $part = shift @parts){
					$part =~ s/\%/\\x/g;
					$part = decode 'ascii-escape', $part;
					push @{$attribs{$key}},$part;
				}
			}
			else {
				$attribs{$key} =~ s/\%/\\x/g;
				$attribs{$key} = decode 'ascii-escape', $attribs{$key};
			}
		}
		return \@data,\%attribs;
	}

}

=head2 id
  Function : get/set the id for a feature
  Example  : $name = $node->id();
             $node->id('new_ID');
=cut

sub id {
    my $node = shift;
    my $val = shift;
    $node->attributes->{ID} = $val if $val;
    return $node->attributes->{ID};
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

# use nesting to allow sub to retain access to private variable
{
	my %pf_conv = ( 0 => 0, 1 => 2, 2 => 1 );

=head2 _frame
  Function : get/set the frame name for a feature
  Example  : $frame = $node->_frame();
             $node->_frame(1);
=cut

	sub _frame {
		my $node = shift;
		my $val = shift;
		if (defined $val && length $val > 0){
			$node->attributes->{_phase} = $pf_conv{$val} ;
		}
		return $pf_conv{$node->attributes->{_phase}};
	}

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
    	if ($_[0]->name() && $_[0]->name eq $name) {
            push @found, $_[0];
            return $retvalue;
        }
        1}});
    return wantarray? @found : @found ? $found[0] : undef;
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
        if ($_[0]->attributes->{$attrib}){
        	my $match  = 0;
        	if (!defined $value) {
        		$match++;
            }
            elsif (ref $_[0]->attributes->{$attrib} eq 'ARRAY'){
            	for (my $i = 0; $i < @{$_[0]->attributes->{$attrib}}; $i++){
            		if ($_[0]->attributes->{$attrib}->[$i] =~ m/^$value$/i){
        				$match++;
            		}
            	}
            }
            elsif ($_[0]->attributes->{$attrib} =~ m/^$value$/i){
            	$match++;
            }
            if ($match > 0){
            	push @found, $_[0];
            	return $retvalue;
            }
        }
        1}});
    return wantarray? @found : @found ? $found[0] : undef;
}


=head2 by_attributes
  Function : returns a scalar or array of features with each of a given set of attributes
             or where each of a set of attributes match specific values
  Example  : @nodes = $gff->by_attributes(['anything','anything_else]);
  			 $node = $gff->by_attributes(['anything','something'],['anythingelse','somethingelse']);
=cut

sub by_attributes {
    my $self = shift;
    my @matches;
    while (my $arrayref = shift){
    	my @nodes;
    	if ($arrayref->[1]){
    		@nodes = $self->by_attribute($arrayref->[0],$arrayref->[1]);
    	}
    	else {
    		@nodes = $self->by_attribute($arrayref->[0]);
    	}
    	if (@matches){
    		my @tmp;
    		for (my $i = 0; $i < @matches; $i++){
    			for (my $j = 0; $j < @nodes; $j++){
    				push @tmp,$matches[$i] if $matches[$i] == $nodes[$j];
    			}
    		}
    		@matches = @tmp;
    	}
    	else {
    		@matches = @nodes;
    	}
    	return unless @matches;
    }
    return wantarray? @matches : $matches[0];
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
    	my $match  = 0;
    	if ($_[0]->attributes->{$attrib}){
    		if (defined $value){
    			if (ref $_[0]->attributes->{$attrib} eq 'ARRAY'){
    				for (my $i = 0; $i < @{$_[0]->attributes->{$attrib}}; $i++){
            			if ($_[0]->attributes->{$attrib}->[$i] =~ m/^$value$/i){
        					$match++;
            			}
            		}
    			}
    			elsif ($_[0]->attributes->{$attrib} =~ m/^$value$/i) {
    				$match++;
       	 		}
       	 	}
       	 	else {
       	 		$match++;
       	 	}
       	}
        if ($match == 0){
            push @found, $_[0];
       	    return $retvalue;
        }
        1}});
    return wantarray? @found : @found ? $found[0] : undef;
}


1;