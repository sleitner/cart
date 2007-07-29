#!/usr/bin/perl

$num_procs = 384;
$num_grid = 512;
$max_level = 7;

$min_local_root_cells = $num_grid**3;
$min_root_cells = $num_grid**3;
$min_particles = $num_grid**3;

$max_local_root_cells = 0;
$max_root_cells = 0;
$max_root_buffer_fraction = 0;
$max_buffer_fraction = 0;
$max_particles = 0;
$max_octs = 0;
$min_octs = $num_grid**3;

@total_level_cells = ();

for ( $proc = 0; $proc < $num_procs; $proc++ ) {
	$filename = sprintf "/gpfs/scratch/uam37/uam37895/CLS1/workload.%03u.dat", $proc;

	open FILE, $filename or die "Unable to open $filename\n";
	while ( $line = <FILE> ) {
		$last = $line;
	}
	close FILE;

	chomp $last;

	@cols = split /\s+/, $last;

	$max_particles = ( $cols[3] > $max_particles ) ? $cols[3] : $max_particles;
	$max_local_root_cells = ($cols[5] > $max_local_root_cells ) ? $cols[5] : $max_local_root_cells;
	$max_root_cells = ( ($cols[5]+$cols[6]) > $max_root_cells ) ? ($cols[5]+$cols[6]) : $max_root_cells;
	$max_root_buffer_fraction = ( ($cols[6]/($cols[5]+$cols[6])) > $max_root_buffer_fraction ) ?
			($cols[6]/($cols[5]+$cols[6])) : $max_root_buffer_fraction;

	$min_particles = ( $cols[3] < $min_particles ) ? $cols[3] : $min_particles;
	$min_local_root_cells = ( $cols[5] < $min_local_root_cells ) ? $cols[5] : $min_local_root_cells;
	$min_root_cells = ( $cols[5]+$cols[6] < $min_root_cells ) ? ($cols[5]+$cols[6]) : $min_root_cells;

	$cells = 0;
	$buffer = 0;
	for ( $level = 0; $level <= $max_level; $level++ ) {
		$cells += $cols[5+2*$level] + $cols[5+2*$level+1];
		$buffer += $cols[6+2*$level];
		$total_cells[$level] += $cols[5+2*$level];
	}
	
	$buffer_octs = $buffer / 8;
	$octs = $cells/8;

	$total_cells += $cells;
	$total_octs += $octs;
	$total_particles += $cols[3];	

	$max_octs = ( $octs > $max_octs ) ? $octs : $max_octs;
	$min_octs = ( $octs < $min_octs ) ? $octs : $min_octs;

	if ( $octs > 0 ) {
		$max_buffer_fraction = ( $buffer_octs/$octs > $max_buffer_fraction ) ? $buffer_octs/$octs : $max_buffer_fraction;
	}
}

print "Total cells = $total_cells\n";
print "Total octs = $total_octs\n";
print "Total particles = $total_particles\n";
print "Num particles: $min_particles -> $max_particles\n";
print "Num root cells: $min_root_cells -> $max_root_cells\n";
print "Num local root cells: $min_local_root_cells -> $max_local_root_cells\n";
print "Num octs: $min_octs -> $max_octs\n";
printf "Max root buffer fraction = %0.3f\n", $max_root_buffer_fraction;
printf "Max buffer fraction = %0.3f\n", $max_buffer_fraction;

$total_cells[0] = $num_grid**3;

for ( $level = 0; $level <= $max_level; $level++ ) {
	print "$level $total_cells[$level]";

	if ( $level < 8 ) {
		$num_leaves = $total_cells[$level] - ($total_cells[$level+1]/8);
		print " $num_leaves\n";
	} else {
		print " $total_cells[$level]\n";
	}
}
