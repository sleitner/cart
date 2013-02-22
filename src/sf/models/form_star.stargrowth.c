#include 'form_star.stargrowth.h'
//snl add includes


void get_com_pos( int icell, double com_pos[nDim], int level ){
    int i, idir, isign;
    int neighbors[num_neighbors];
    double dcom_pos[nDim];
    double sum[nDim];
    
    cell_center_position(icell,com_pos); /* start at cell position */

    cell_all_neighbors(icell,neighbors);
    for(idir=0;idir<nDim; idir++){ dcom_pos[idir]=0; sum[idir]=cell_gas_density(icell); } 
    for(i=0;i<num_neighbors; i++){
	idir = i/2;
	isign = 1 - (i % 2)*2; // even is +
	dcom_pos[idir] += isign*cell_size[level]*cell_gas_density(neighbors[i]); 
	sum[idir] += cell_gas_density(neighbors[i]);
    }
    /* now assign COM position but stay within cell*/
    for(idir=0;idir<nDim; idir++){
	com_pos[idir] += (
	    dcom_pos[idir] > 0 ?
	    MIN( dcom_pos[idir]/sum[idir], cell_size[level]/2.01) :
	    MAX( dcom_pos[idir]/sum[idir],-cell_size[level]/2.01) );
    }
}
void grow_star_particle( int ipart, float delta_mass, int icell, int level) {
    int i;
    double add_mass;
    double new_density;
    double density_fraction, thermal_pressure;
    double sum_mass, pmass_orig;
    double com_pos[nDim];
	
    cart_assert( ipart < num_star_particles );
    add_mass = MIN( delta_mass, 0.667*cell_volume[cell_level(icell)]*cell_gas_density(icell) );
/*     get_com_pos(icell, com_pos, level ); */
/*     cart_assert( cell_contains_position(icell,com_pos) );  */
/*     for ( i = 0; i < nDim; i++ ) { */
/*      particle_x[ipart][i] = (com_pos[i]*add_mass + particle_x[ipart][i]*pmass_orig)/sum_mass; */
/*     } */
    /* add mass [at COM of neighbors] with cell momentum; */
    pmass_orig = particle_mass[ipart];
    sum_mass = pmass_orig + add_mass;
    cart_assert( cell_contains_position(icell,particle_x[ipart]) ); 
    for ( i = 0; i < nDim; i++ ) {
	particle_v[ipart][i] = 
	    ( cell_momentum(icell,i) / cell_gas_density(icell) * add_mass +
	      particle_v[ipart][i]*pmass_orig ) /sum_mass;
    }
    particle_mass[ipart] += add_mass;
    star_initial_mass[ipart] += add_mass; /* this is used for feedback */
    
#ifdef ENRICHMENT
    if(sf_metallicity_floor>0.0 && cell_gas_metal_density_II(icell)<sf_metallicity_floor*constants->Zsun*cell_gas_density(icell)){
	cell_gas_metal_density_II(icell) =  sf_metallicity_floor*constants->Zsun*cell_gas_density(icell);
    }
    star_metallicity_II[ipart] = 
	( cell_gas_metal_density_II(icell) / cell_gas_density(icell) * add_mass +
	  star_metallicity_II[ipart] * pmass_orig) / sum_mass ;
#ifdef ENRICHMENT_SNIa
    star_metallicity_Ia[ipart] = 
	( cell_gas_metal_density_Ia(icell) / cell_gas_density(icell) * add_mass +
	  star_metallicity_Ia[ipart] * pmass_orig) / sum_mass;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
    
    /* adjust cell values */
    new_density = cell_gas_density(icell) - add_mass * cell_volume_inverse[level];
    density_fraction = new_density / cell_gas_density(icell);
    
    /*
    // NG: this is to allow non-thermal pressure contribution
    */
    thermal_pressure = MAX((cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell),0.0);
    cell_gas_pressure(icell) = MAX(0.0,cell_gas_pressure(icell)-thermal_pressure);
    
    cell_gas_density(icell) = new_density;
    cell_gas_energy(icell) *= density_fraction;
    cell_gas_internal_energy(icell) *= density_fraction;
    cell_momentum(icell,0) *= density_fraction;
    cell_momentum(icell,1) *= density_fraction;
    cell_momentum(icell,2) *= density_fraction;
    
    cell_gas_pressure(icell) += thermal_pressure*density_fraction;
    
    for ( i = 0; i < num_chem_species; i++ ) {
	cell_advected_variable(icell,i) *= density_fraction;
    }
}
