#!/usr/bin/perl

$directory = $ARGV[0];

if ( defined($ARGV[1]) ) {
	$print_step = $ARGV[1];
} else {
	$print_step = -1;
}

open FILE, "$directory/timing.000.log" or die "Unable to load timing file!\n";

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

	printf "Step: %u total time: %.2f\n", $step, $current_step;

	$headerline = sprintf "%30s %10s", "Levels:", "Global";
	for $level ( 1 .. $num_levels-1 ) {
		$headerline = sprintf "%s %8u  ", $headerline, $level-1;
	}
	$headerline = sprintf "%s  %8s %5s\n", $headerline, "Total", "%";

	print $headerline;	
	print (("-" x length($headerline)) . "\n");

	for $timer ( 0 .. $num_timers-1 ) {
		printf "%30s", $timer_labels[$timer];
		$sum = 0;
		for $level ( 0 .. $num_levels-1 ) {
			printf " %10.3f", $current_timers[$level*$num_timers + $timer];
			$sum += $current_timers[$level*$num_timers + $timer];
		}
		printf "%10.2f %6.2f\n", $sum, ($current_step > 0) ? 100.0*$sum/$current_step : 0.0;
	}
}

close FILE;
