#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl table.extendLines.pl [options] table.txt indexes key1 ... > table.extended.txt

Options: -h       display this help message

EOF
}
my ($tableFile, $valueIndexes, @keysList) = @ARGV;
my %valueIndexHash = map {$_ => 1} (my @valueIndexList = eval($valueIndexes));
if(my @keyList = map {split(/,/, $_)} @keysList) {
	open(my $reader, $tableFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my @constantTokenList = @tokenList[grep {!defined($valueIndexHash{$_})} 0 .. $#tokenList];
		print join("\t", @constantTokenList, $_), "\n" foreach(map {join("\t", $keyList[$_], $tokenList[$valueIndexList[$_]])} 0 .. $#valueIndexList);
	}
	close($reader);
} else {
	open(my $reader, $tableFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my @constantTokenList = @tokenList[grep {!defined($valueIndexHash{$_})} 0 .. $#tokenList];
		print join("\t", @constantTokenList, $_), "\n" foreach(@tokenList[@valueIndexList]);
	}
	close($reader);
}
