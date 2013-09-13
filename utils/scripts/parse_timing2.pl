#!/usr/bin/perl

$directory = $ARGV[0];

if ( defined($ARGV[1]) ) {
	$print_step = $ARGV[1];
} else {
	$print_step = -1;
}

my @selected_timers = ( "hydro", "gravity", "density", "refinement", "communication", "work", "RT_cooling", "RT_level_update", "level_total" );

($file) = glob("$directory/timing.0*.log");
open FILE, $file or die "Unable to load timing file!\n";

$line = <FILE>;
($num_levels,$num_timers) = ($line =~ /^\# (\d+) levels (\d+) timers/ );

$num_real_levels = $num_levels -1;
print "$num_real_levels levels, $num_timers timers\n";

$line = <FILE>;
chomp $line;
($c, $step, $t, $aexp, $total, @timer_labels) = split /\s+/, $line;

$old_total = 0.0;
@old_timers = ();

while ( $line = <FILE> ) {
	next if $line =~ /^\#/;
	chomp $line;
	($step, $t, $aexp, $total, @cols) = split /\s+/, $line;

	if ( $total > $old_total ) {
		$current_step = $total - $old_total;
		for $i ( 0 .. @cols-1 ) {
			$current_timers[$i] = $cols[$i] - $old_timers[$i];
		}
	} else {
		$current_step = $total;
		for $i ( 0 .. @cols-1 ) {
			$current_timers[$i] = $cols[$i];
		}
	}

	$old_total = $total;
	@old_timers = @cols;

	# subtract lower level timers from higher level
	for $level ( 0 .. $num_levels-2 ) {
		$current_timers[($level+1)*$num_timers-1] -=
			$current_timers[($level+2)*$num_timers-1];
	}

	next if $print_step != -1 and $step != $print_step;

	print "Step: $step total time: $current_step\n";

	$headerline = sprintf "%30s", "Timer";
	$headerline = sprintf "%s  %8s %5s\n", $headerline, "Total", "%";

	print $headerline;	
	print (("-" x length($headerline)) . "\n");

	for $timer ( 0 .. $num_timers-1 ) {

	    my $sel = 0;
	    for $i ( 0 .. $#selected_timers )
	    {
		if ( $timer_labels[$timer] eq $selected_timers[$i] )
		{ 
		    $sel = 1;
		}
	    }

	    if ( $sel ) 
	    {

		printf "%30s", $timer_labels[$timer];
		$sum = 0;
		for $level ( 0 .. $num_levels-1 ) {
			$sum += $current_timers[$level*$num_timers + $timer];
		}
		printf "%10.2f %6.2f\n", $sum, ($current_step > 0) ? 100.0*$sum/$current_step : 0.0;
	    }
	}
}

close FILE;
