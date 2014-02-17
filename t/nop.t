#!perl 
use 5.012;
use strict;
use warnings FATAL => 'all';
use Test::Simple tests => 3;
#use PDL;
use PDL;
use PDL::Dims ; #qw/initdim/;
use PDL::NiceSlice;

$a=sequence(6)->reshape(2,3);
initdim ($a,'x',2);
initdim ($a,'y',3);
$b=zeroes(3,4,2,3);
initdim( $b,'y',3);
initdim( $b,'a',4);
initdim( $b,'u',2);
initdim( $b,'b',3);
#ok((sln($a,x=>0,y=>2)**2)==$a(0,2;-)**2,'sln');
ok((nop($a,'rotate','y',1)->(1,0))==5,'rotate');
ok(max (ncop($a,$b,'plus',0))==5,'ncop');
ok((nreduce ($a,'add','y','x'))==15,'nreduce');
#diag( "Testing PDL::Dims $PDL::Dims::VERSION, Perl $], $^X" );
