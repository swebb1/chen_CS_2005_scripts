#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %opt = ();
            
GetOptions(\%opt, "in=s", "chr=s","remove");


my %chrmap;
open (CHR, "$opt{chr}") || die "cannot open $opt{chr}";
while (<CHR>){
    chomp;
    my @l = split (/\t/,$_);
    $chrmap{$l[0]} = $l[1];
}
close(CHR);

open (IN, "$opt{in}") || die "cannot open $opt{in}";

while (<IN>){
    if(m/^#/){
	print $_;
    }
    else{
	my @line = split (/\t/,$_);
	my $chromo = $line[0];
	if(defined $chrmap{$chromo}){
	    $line[0] = $chrmap{$chromo};
	    print join("\t", @line);
	}
	else{
	    unless(defined $opt{remove}){
		print join("\t", @line);
	    }
	}
    }
}
close(IN);
