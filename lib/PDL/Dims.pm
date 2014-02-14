#!/usr/bin/perl

package PDL::Dims;

=head1 NAME

PDL::Dims - work on named dimensions and meaningful ranges

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.001';

use strict;


use parent 'Exporter','PDL';
#use base 'PDL';
#@PDL::Dims::ISA=qw/PDL Exporter/;
use PDL::NiceSlice;
use Storable qw(dclone);
use PDL;
use 5.012;
our @EXPORT=qw(i2pos pos2i initdim idx didx dimsize dimname sln rmdim davg dinc dmin dmax nreduce nop ncop copy_dim );



# returns a list of all dims for one parameter -- intended as an auxillary for internal use only

sub _dimpar {
	my $self=shift;
	my $p=shift; # par
	#say "All of ".%{$self->hdr}.", $p";
	return unless $p;
	my @s;
	for my $i (@{dimname($self)}) {
		say "name $i $p";
		barf ("Unknown dim $i") unless (ref ($self->hdr->{$i}) eq  'HASH');
		push @s,$self->hdr->{$i}->{$p};
	#	say "name $i: ",$self->hdr->{$i}->{$p};
	}
	#say "Dims ".@{dimname($self)}."Par $p: @s";
	return @s;
}
sub dimsize { # sets size for a dim
	my $self=shift;
	my $d=shift; # dimname
	my $v=shift; # value
	return #([values %{$self->hdr->{dimsize}}]
#		,[keys %{$self->hdr->{dimsize}}],
		#)
		[_dimpar($self,'dimsize')]
	unless ($d) ;
	if (defined $v) {
		barf ("Unknown dim $d") unless (ref ($self->hdr->{$d}) eq  'HASH');
		int $v;
		$self->hdr->{$d}->{dimsize}=$v; 
		idx($self,$d,$v-1) if (idx($self,$d)>=$v-1);
	}
	return $self->hdr->{$d}->{dimsize}; #, $r{$d};
}# get/set name of a dim by number
sub dimname {
	my $self=shift;
	my $d=shift; # dim number
	my $n=shift; # new name
	#die unless $self->isa('PDL');
	#say keys %{$self->hdr} unless defined $d;
	return $self->hdr->{dimnames} unless defined $d;
	#barf "Unknown dim $d" unless (ref $self->hdr->{dimnames} eq  'HASH');
	$self->hdr->{dimnames}->[$d]=$n if defined $n;
	return $self->hdr->{dimnames}->[$d];
}

sub davg {
	my $self=shift;
	my $d=shift; # dimension
	my $r=(shift || [0,-1,1]); # range [a,b,step] 
	#say "dataset $self dim $d range ",@{$r};
	return ($$self->mv(didx($self,$d),0)->($$r[0]:$$r[1]:$$r[2])->average);
}
# boundaries of dims in appropriate units (mm, s, ..)

sub dinc {
	my $self=shift;
	my $d=shift; # dimname
	my $v=shift; # max 
	return [_dimpar($self,'diminc')] unless $d;
	if (defined $v) {
		barf ("Unknown dim $d") unless (ref $self->hdr->{$d} eq  'HASH');
		$self->hdr->{$d}->{diminc}=$v; 
		$self->hdr->{$d}->{spacing}=1;
	}
	return $self->hdr->{$d}->{diminc}; #, $r{$d};
}

sub dmax {
	my $self=shift;
	my $d=shift; # dimname
	my $v=shift; # max 
	return [_dimpar($self,'dimmax')] unless $d;
	#say "$d ".$self->hdr;
	barf ("Unknown dim $d") if (defined $v and ref $self->hdr->{$d} ne  'HASH');
	$self->hdr->{$d}->{dimmax}=$v if defined ($v);
	return $self->hdr->{$d}->{dimmax}; #, $r{$d};
}

sub dmin {
	my $self=shift;
	my $d=shift; # dimname
	my $v=shift; # min 
	return [_dimpar($self,'dimmin')] unless $d;
	barf ("Unknown dim $d") if (defined $v and ref $self->hdr->{$d} ne  'HASH');
	#barf "Unknown dim $d" unless (ref $self->hdr->{$d} eq  'HASH');
	$self->hdr->{$d}->{dimmin}=$v if defined ($v);
	return $self->hdr->{$d}->{dimmin}; #, $r{$d};
}


sub didx { # set/get index - complementary to dimname
	my $self=shift;
	my $d=shift; # dimname
	return [_dimpar($self,'pos')] unless $d;
	if (ref $d eq 'ARRAY') {

		@$d;
	} else {
		my $n=shift; # new position
		barf ("Unknown dim $d") if (defined $n and ref $self->hdr->{$d} ne  'HASH');
		#barf "Unknown dim $d" unless (ref $self->hdr->{$d} eq  'HASH');
		$self->hdr->{$d}->{pos}=$n if defined $n;
		#say "type $self idx $d ".$self->hdr->{$d}->{pos};
		return $self->hdr->{$d}->{pos};
	}
}

#transformations between piddle and world coords. 
sub i2pos{
	my $self=shift; # dataset
	my $d=shift || return; # dimname
	my $i=shift ; #value
	return $i*(dmax($self,$d)-dmin($self,$d))/dimsize($self,$d)+dmin($self,$d);
}

sub pos2i{
	my $self=shift;
	my $d=shift || return;
	my $i=shift ;
	return rint (($i-dmin($self,$d))/(dmax($self,$d)-dmin($self,$d))*dimsize($self,$d)-0);
}

sub nreduce { # analogue to PDL::Reduce, perform operation on named rather than numbered dims
	my $self=shift;
	my $op=shift; # operation - passed on to reduce
	require PDL::Reduce;
	my @d=map {didx($self,$_)} @_;
	say "reduce $op @d (@_)";
	return $self->reduce($op,@d);
}

sub nop { # wrapper to functions like sin, exp, rotate operating on one named dimension
	my $self=shift;
	my $op=shift;
	my $dim=shift;
	#my $arg=shift; # suitably shaped second argument 
	$self=$self->mv(didx($self,$dim),0);
	my $res;
	if ($op eq 'rotate'){
		my $s=shift;
		say "schifing $s";
		$res=$self->rotate($s);
	} else {
		$res=$self->$op (@_);
	}
	return ($res->mv(0,didx($self,$dim)));
}

sub ncop { # named combine operation -- +*-/ ... 
	my $self=shift;
	my $other=shift;
	my $op=shift;
	my @d=@{dimname($self)};
	my @e=@{dimname($other)};
	my $m=0;
	my @nother; # new order in other
	my @nself; # new order in self
	my @tself; # thread dims
	my @tother; 
	my @add;
	my @co=();
	my $i=0;
	my $j=0;
	for my $ds (@d) {
		my $n=didx($self,$ds);
		if (defined ($m=didx($other,$ds))) {
			push @nself,$n;
			push @nother,$m;
	#		say "$ds $m $n";
	#		say "self $ds $m $n i $i";
			push @co,$i+0;
			$i++;
		} else {
			push @tself,$n;
	#		say "$ds $m $n i: $i ", ($other->ndims+$i);
			push @co,($other->ndims+$j);
		}
	}
	for my $ds (@e) {
		my $n=didx($other,$ds);
		if (eval {$m=didx($self,$ds)}) {
			1;
		} else {
			push @tother,$n;
			push @add,$ds;
			push @co,$i+0;
	#		say "other $ds $m $n";
			$i++;
		}
	}
	#say "@nself, @tself; other @nother, @tother.";
	#say "Co: @co";
	my $ns=$self->sever->reorder(@nself,@tself); 
	my $no=$other->sever->reorder(@nother,@tother);
	for my $n (0..$#tother) { # fill in dummy dims
		$ns=$ns->dummy($#nself+1,1);
	}
	say $ns->info,$no->info;
	say $self->info,$other->info;
	my $res=$ns->$op($no,@_);
	say $res->info;
	$res->sethdr($self->hdr_copy);
	#say @{dimname($res)};
	#say @{dimname($other)};
	#say @{dimname($self)};
	$other->hdr;
	my $i=$self->ndims;
	for my $ds (0..$#add) {
		say $ds;
		copy_dim($other,$res,$add[$ds],$i);
		$i++;
	}
	#@co=@{didx($res)};
	#say "@co";
	$res=$res->reorder(@co);
	say $res->info;
	return $res;
}

sub sln { # returns a slice by dimnames and patterns
	my $self=shift;
	my %s=@_; # x=>'0:54:-2', t=>47, ... #
	my $str;
	#say "%s";
	#say ("dimnames @{$self->dimname}");
	my (@n)=@{$s{names}||dimname($self)};
	for my $i (0.. $#n) {
		#say "$i $n[$i] $idx[$i]";

		$str.=$s{$n[$i]} if ($s{$n[$i]});
		$str.=',' unless ($i==$#n);
	}
	#say $str;
	return $self->slice($str);
}

sub initdim {
	my $self=shift;
	my $d=shift || return ;
	#say "Init dim $d ...";
	$self->hdr->{ndims}=0 unless ($self->hdr->{ndims});
	warn "$d is defined!" if (ref $self->hdr->{$d} eq  'HASH');
	my %p=@_;
	#say "pars: ",%p;
	$self->hdr->{$d}=\%p;
	#say "Creating dim $d at pos. $p{pos}; Ndims ".$self->hdr->{ndims};
	if ((not $p{pos}) or ($p{pos}>$self->hdr->{ndims})) {
		$p{pos}=$self->hdr->{ndims};
	}
	if ($p{pos}<$self->hdr->{ndims}) {
		for (my $i=$self->hdr->{ndims}-1;$i>=$p{pos};$i--) {
			dimname($self,$i,dimname($self,$i-1)); # shift the position up!
			didx($self,dimname($self,$i-1),$i);
		}		
	}
	$p{size}=$self->dim($p{pos}) unless ($p{size});
	#say "Pos $p{pos}";
	didx ($self,$d,$p{pos});
	dimname ($self,$p{pos},$d);
	$self->hdr->{$d}=\%p;
	dimsize ($self,$d,$p{size}||1);
	dmin ($self,$d,$p{min}||0);
	if ($p{inc} and $p{max}) {
		barf ("Increment and maximum don't fit! ($self $d $p{min} $p{max} $p{inc} ".dimsize($self,$d))
			#." ".$p{max}-$p{min}-(dimsize($self,$d)-1)*$p{inc)
#			),mb::Error) 
			unless ($p{max}-$p{min} - (dimsize($self,$d)-1)*$p{inc} < 1e-8);
	} elsif ($p{inc}) {
		$p{max}=$p{inc}*dimsize($self,$d)+$p{min};
	} elsif ($p{max}) {
		$p{inc}=($p{max}-$p{min})/((dimsize($self,$d)-1)||1);
	} else {
		$p{max}=dimsize($self,$d)-1;
		$p{inc}=1;
	}
	#$p{inc}=($p{inc}||1);
	dinc ($self,$d,$p{inc});
	dmax ($self,$d,$p{max}); #||(dimsize($self,$d)-1)*$p{inc};
	$self->hdr->{ndims}++; 
	idx($self,$d,($p{index}||dmin($self,$d)));
	if (ref $p{vals}) {
		vals ($self,$d,$p{vals}); 
	} else {
		#say [list (sequence (dimsize($self,$d))*$p{inc}+$p{min})];
		vals ($self,$d,[list (sequence (dimsize($self,$d))*$p{inc}+$p{min})]);
		$self->hdr->{$d}->{spacing}=1;
	}
}

sub copy_dim {
	my $old=shift;
	my $new=shift;
	my $dim=shift;
	#say "old: $old; new: $new; dim: $dim";
	my $d=dclone($old->hdr->{$dim});
	#say "old $old new $new dim %$d";
	#say @{dimname $new};
	$$d{pos}=shift;
	initdim($new,$dim,%$d);
}

sub rmdim {
	my $self=shift;
	my $d=shift;
	#say "removing $self, $d";
	splice $self->hdr->{dimnames},didx($self,$d),1; # cut out the array
	for my $i (didx($self,$d)..$self->hdr->{ndims}-1) { 
		didx($self,dimname($self,$i),$i);	#update position of following dims
	}
	delete $self->hdr->{$d};
	die "This should be undefined! ",$self->hdr->{$d} if defined ($self->hdr->{$d});
	$self->hdr->{ndims}--;
}

sub idx {
	#my $self=shift;
	my $self=shift;
	die "I don't have data to work on (idx)" unless defined $self->hdr;
	#say "$self type array ";
	my $index=shift;
	return [_dimpar($self,'index')] unless $index;
	if (defined (my $v=shift)) {
		$v<0? 0: $v;
		$v>=dimsize($self,$index)? dimsize($self,$index)-1 : $v;
		$self->hdr->{$index}->{index}=$v ;
	}
	return $self->hdr->{$index}->{index}; 
}

sub vals { #individual values of dims -- like the position along an axis in world coords., echo times
	my $self=shift;
	my $d=shift; #dim
	return unless $d;
	$self->hdr->{$d}->{vals}=[] unless defined $self->hdr->{$d}->{vals};
	#say "Vals: $d ",@{$self->hdr->{$d}->{vals}};
	my $i=shift; #index
	if (defined $i) {
		if (ref $i eq 'ARRAY' and dimsize($self,$d) == $#{$i}+1) {
			$self->hdr->{$d}->{vals}=$i ;
			$self->hdr->{$d}->{spacing}=0;
		} else {
			my $v=shift; #value
			if ( defined $v) { 
				$self->hdr->{$d}->{vals}->[$i]=$v ;
				$self->hdr->{$d}->{spacing}=0;
			}
		}
		return $self->hdr->{$d}->{vals}->[$i];
	} else {
		return $self->hdr->{$d}->{vals};
	}
}

1;
#BEGIN {
#        if ($_[0] eq q/-d/) { 
#		require Carp; 
#		$SIG{__DIE__} = sub {print Carp::longmess(@_); die;}; 
#	} 
#}


=head1 SYNOPSIS

    use PDL::Dims;

	
    initdim ($piddle,'x',pos=>0);


=head1 DESCRIPTION

This module is an extension to PDL by giving the dimensions a meaningful name,
values or ranges. The module also provides wrappers to perform most PDL
operations based on named dimensions. Each dimension is supposed to have a unique name.



Each dim has its own hash with several important keys; names can be set (x,y,z,t, ...)

=over 

=item * dimpos - hash of index names and their position (i.e. x=>0, z=>2, t=>3, y=>1 ...)

=item * dimnames - the opposite of 'dimpos'; an array naming the indices

=item * dims - hash of index sizes (by name)

=item * index - hash of current values of each index (by name)

=back




=head1 SUBROUTINES/METHODS

The module has two different types of functions. 

First, there are several utitility functions to manipulate and retrieve the dimension info. 
Currently, they only work on the metadata, no operations on the piddle are performed.
It is your responsibility to call the appropriate function whenever piddle dims change.

Then there are functions to perform most PDL operations based on named dimensions instead
of numbered dims. Wrappers to slice (sln), reduce (nreduce) and other functions 
(nop, ncop) are implemented. 

The following parameters are currently used by PDL::Dims. It is *strongly* discouraged
to access them directly except during calls to initdim. Use the functions listed below
to manipulate and retrieve values.

They are usually called like that:

	func ($piddle , [$dim, [value]]);

if no dim is given, an array reference to all values are returned.

Otherwise, it returns the given value for a particular $dim, setting it to $value if defined.

=over 

=item * dimsize - set/get dimension size of piddle by name

=item * dimname - get/set dimension name of piddle by index

=item * idx - set/get current index values by name

=item * didx - get/set dimension names by index

=item * dmin/dmax/dinc - get/set min/max/increment values in appropriate units (mm/s/ms...)

=item * unit - NOT YET IMPLEMENTED. In the future, each assignment or operation
should be units-avare, e.g. converting autmatically from mm to km or feet,
preventing you from adding apples to peas ...

=back


Use the following for creating/deleting dimensions (meta data).


=head2 initdim - 

initializes a dimenson

usage: 
	initdim($piddle,$dimname,%args);

Arguments are the fields listed above.

=head2 rmdim 

removes a dminenso 

	rmdim($piddle,$dim);

=head2 * copy_dim 

copies the diminfo from one piddle to the next.

	copy_dim($a,$b,$dim,$pos);

calls initdim on $b with parameters from $a. You can supply a new position.

=head2 ncop

operates on two piddles, combining them by names. Equally named dimensions have
to obey the usual threding rules. For opertators like '+' use the named version
and an additional argument 0.

usage:
	$c=ncop($a,$b,'operation',@args);

=head2 nreduce

wrapper around reduce calling with names instead.

=head2 sln

	sln ($piddle,%slices);

perform slicing based on names.

Example:
	$sl=sln($a,x=>'0:10:2',t=>'4:',);

You can omit dims; it parses each key-value pair and constructs a slice string,
so what typically works for slice, works here.

=head2 nop

perform non-aggregate function on a dimension.

usage:
	$res=nop($a,$operation,@args);

=head1 AUTHOR

Ingo Schmid




=head1 LICENSE AND COPYRIGHT

Copyright 2014 Ingo Schmid.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

1; # End of PDL::Dims
