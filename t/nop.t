#!perl 
use 5.012;
use strict;
use warnings FATAL => 'all';
#use PDL;
use PDL;
#use lib 'lib';
use PDL::Dims ; #qw/initdim/;
use PDL::NiceSlice;
use Test::Simple tests => 7;

#sub ok{}
$a=sequence(6)->reshape(2,3);
initdim ($a,'x',size=>2);
initdim ($a,'y',);
$b=zeroes(3,4,2,3);
initdim( $b,'y',);
initdim( $b,'a',vals =>['a','b','c','d']);
initdim( $b,'u',);
initdim( $b,'b',);
ok((my @x=pos2i($b,'a','c')==42),'init');
ok(pos2i($b,'a','c'),'pos2i');
ok(sclr(sln($a,x=>0,y=>1))**2==4 ,'sln'); #$a(0,1;-)**2,'sln');
ok(sclr (nagg($a,'sumover','y')->(1))==9,'nagg');
#undef $a,$b;
ok(sclr(nop($a,'rotate','y',1)->(1,0))==5,'rotate');
ok(max (ncop($a,$b,'plus',0))==5,'ncop');
ok((nreduce ($a,'add','y','x'))==15,'nreduce');
#diag( "Testing PDL::Dims $PDL::Dims::VERSION, Perl $], $^X" );
