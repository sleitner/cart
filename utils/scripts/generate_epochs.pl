#!/usr/bin/perl

open EPOCH, "> HC/epochs.dat";

foreach $file ( <HC/hlist_*.dat> ) {
	($aexp) = ( $file =~ /hlist_(.*)\.dat/ );

	print EPOCH "$aexp\n";
}

close EPOCH;
