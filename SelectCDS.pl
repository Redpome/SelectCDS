#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Bio::Translator;
use Bio::Util::DNA 'cleanDNA';

if(scalar @ARGV !=2){
	print "\nUsage: perl SelectCDS.pl input.fa output.fa\n\n";
	exit;
}
unless ( open(E, "<$ARGV[0]") ) {
	print "Cannot open file \n\n";
	exit;
}
open O,">$ARGV[ 1 ]";
open T,">$ARGV[ 0 ].err";
my $t = Bio::Translator->new(1);
my ($id,$seq);
$/=">";
my $Scod ="ATG";
my $Ecod1;
my $Ecod2;
my $Ecod3;
my $Site;
my $seqnum;
my $cod3;
my $seqnumT;
<E>;
while(<E>){
	chomp;
	my ($id,$seq)=split(/[\n\r]+/,$_,2);
	$seq=~s/[\r\n\s]//g;
	$seq=uc($seq);
	$seqnum=length($seq);
	$seqnumT=$seqnum-3;
    $Site = index ($seq,$Scod);
	$cod3 = $seqnum % 3;
	$Ecod1 =rindex ($seq,"TAG");
    $Ecod2 =rindex ($seq,"TGA");
    $Ecod3 =rindex ($seq,"TAA");
    if ($Site eq 0) {
        if ($cod3 eq 0) {
			if ($Ecod1 eq $seqnumT || $Ecod2 eq $seqnumT || $Ecod3 eq $seqnumT) {
				my $pep_ref = $t->translate( cleanDNA( \$seq) );
				$$pep_ref =~ s/(.{1,60})/$1\n/g;
				my$aa=$$pep_ref;
				$aa=~s/[\r\n\s]//g;
				my$al=length($aa);
				$al=$al*3;
				my$end="*";
				my @Pro=split(/\*/,$aa);
				my$an=grep/^*/,@Pro;
				if (($an eq 1) and ($al eq $seqnum)) {
                     $seq=~s/...$//g;
					 print O ">$id\n$seq\n";
                }
                else {
					 print T "$id\tInsert TAG|TAA|TGA\n";
				}
            }
			else {
				print T "$id\tNot TAG|TAA|TGA end\n";
			}
        }
        else {
			print T "$id\tLength 3 times error\n";
		}
    }
    else {
		print T "$id\tATG error\n";
	}
}
close E;
close O;
close T;
