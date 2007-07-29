#!/usr/bin/perl

###########################################################
#  Changelog:
#   5/18/07 - Added current physics and updated calculation
#		of static memory usage.
###########################################################

if ( $ARGV[0] ) {
	open DEFS, $ARGV[0] or die "Unable to open $ARGV[0]";
} else {
	open DEFS, "defs.h" or die "Unable to open defs.h";
}

@defines = ( "num_root_grid_refinements", "num_octs", 
	"num_particles", "num_star_particles", "num_tracers", "HYDRO_TRACERS",
	"PARTICLES", "GRAVITY","HYDRO", "STARFORM" );

$regex = "\#define\\s*(" . (join "|", @defines) . ")\\s*(.*)\\s*";

while ( $line = <DEFS> ) {
	chomp $line;
	$line =~ s/\/\*.*\\*//g;
	($def, $value) = ( $line =~ /$regex/g );
	next if $def eq "";

	$defs{$def} = $value;
}

$num_grid = 2**$defs{"num_root_grid_refinements"};
$num_root_cells = $num_grid**3;
$nDim = 3;
$num_children = 2**$nDim;

# add constants
foreach $k ( keys %defs ) {
	$defs{$k} =~ s/num_grid/$num_grid/g;
	$defs{$k} =~ s/num_root_cells/$num_root_cells/g;
	$defs{$k} =~ s/nDim/$nDim/g;
	$defs{$k} =~ s/num_children/$num_children/g;
}

# allow defines to contain other defines as values
foreach $k ( keys %defs ) {
	foreach $k2 ( keys %defs ) {
		next if $k2 eq $k;

		if ( $defs{$k} =~ /[\s+-\/\*(]$k2[\s+-\/\*)]/ ) {
			$defs{$k} =~ s/$k2/$defs{$k2}/;
		}
	}

	if ( $defs{$k} ne "" ) {
		$defs{$k} = eval $defs{$k};
	} else {
		$defs{$k} = 1;
	}
}

# get values into local variables
$num_octs = $defs{"num_octs"} || 0;
$num_cells = $num_children*$num_octs || 0;

print "num_octs = $num_octs\n";
print "num_cells = $num_cells\n";

if ( $defs{"PARTICLES"} ) {
	$num_particles = $defs{"num_particles"} || 0;
	print "num_particles = $num_particles\n";
} else {
	$num_particles = 0;
}

if ( $defs{"STARFORM"} ) {
	$num_star_particles = $defs{"num_star_particles"} || 0;
	print "num_star_particles = $num_star_particles\n";
} else {
	$num_star_particles = 0;
}

if ( $defs{"HYDRO_TRACERS"} ) {
	$num_tracers = $defs{"num_tracers"} || 0;
	print "num_tracers = $num_tracers\n";
} else {
	$num_tracers = 0;
}

print "\n";

$int_size = 4;
$float_size = 4;
$double_size = 8;
$long_size = 8;

if ( $defs{"GRAVITY"} ) {
	if ( $defs{"HYDRO"} ) {
		if ( $defs{"PARTICLES"} ) {
			$num_grav_vars = 4+$nDim;
		} else {
			$num_grav_vars = 3+$nDim;
		}
	} else {
		if ( $defs{"PARTICLES"} ) {
			$num_grav_vars = 3+$nDim;
		} else {
			$num_grav_vars = 2+$nDim;
		}
	}

	$num_refinement_vars = 0; # refinement shares space with GRAVITY
} else {
	$num_grav_vars = 0;
	$num_refinement_vars = 2;
}

if ( $defs{"PARTICLES"} ) {
	$particle_list = 1;
} else {
	$particle_list = 0;
}

if ( $defs{"HYDRO"} ) { 
	$num_hydro_vars = 5+$nDim;
} else {
	$num_hydro_vars = 0;
}

if ( not $defs{"STARFORM"} ) {
	$num_star_particles = 0;
}

if ( not $defs{"HYDRO_TRACERS"} ) {
	$num_tracers = 0;
}

$num_vars = $num_grav_vars+$num_hydro_vars+$num_refinement_vars;
print "num_vars per cell = $num_vars\n";

$bytes_per_oct = 1 * $int_size +  	# parent_cell
		 1 * $int_size +  	# level
                 2*$nDim*$int_size +	# neighbors
                 1 * $int_size +	# parent_root_sfc
                 2 * $int_size +	# next,prev
		 $nDim * $float_size 	# pos
		;

$bytes_per_cell =	1 * $int_size	+		# oct_child 
			$num_vars * $float_size +	# variables
			($num_hydro_vars-2)*$float_size+# backup_hvars
			1 * $float_size +		# ref
			$particle_list * $int_size;	# cell_particle_list

$bytes_per_particle =	2*$nDim*$double_size +	# x and v
		     	1 * $int_size + 	# id
			1 * $int_size +		# level
			1 * $float_size +	# mass
			1 * $double_size +	# particle_pot
			2 * $float_size +	# t, dt
			2 * $int_size		# particle_list
			;

$bytes_per_star = 	1 * $float_size +	# tbirth
			1 * $float_size + 	# initial mass
			1 * $float_size +	# ZII
			1 * $float_size		# ZIa
			;

$bytes_per_tracer = 	$nDim*$double_size +	# x
			1 * $int_size + 	# id
			2 * $int_size		# list
			;

$master_node_bytes_per_root_cell =	2 * $float_size;	# density/potential

print "bytes per oct = $bytes_per_oct\n";
print "bytes per cell = $bytes_per_cell\n";

if ( $defs{"PARTICLES"} ) {
print "bytes per particle = $bytes_per_particle\n";
}

if ( $defs{"STARFORM"} ) {
print "bytes per star particle = $bytes_per_star\n";
}

if ( $defs{"HYDRO_TRACERS"} ) {
print "bytes per tracer particle = $bytes_per_tracer\n";
}

print "\n";

$total_static =	$bytes_per_oct * $num_octs + 
		$bytes_per_cell * $num_cells +
		$bytes_per_particle * $num_particles +
		$bytes_per_star * $num_star_particles +
		$bytes_per_tracer * $num_tracer_particles;

$total_bytes = $total_static+$total_dynamic;

$print_total_static = print_format($total_static);
$print_total_dynamic = print_format($total_dynamic);
$print_total_bytes = print_format($total_bytes);

print "Total static memory: $print_total_static\n";
print "Total dynamic memory: $print_total_dynamic\n";
print "Total memory required: $print_total_bytes\n";
 
sub print_format {
	($bytes) = @_;

	if ( $total_bytes >= (1024)**2 and $total_bytes < (1024)**3 ) {
		return sprintf("%5.2fM", ($bytes/1024**2) );
	} elsif ( $total_bytes >= (1024)**3 ) {
		return sprintf("%5.2fG", ($bytes/1024**3) );
	} else {
		return sprintf("%5.2f", $bytes );
	}
}
