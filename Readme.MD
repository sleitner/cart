## Edit Notes:

# Patterns
*   particle_id
*   particle_species_num
*   particle_species_indices
*   num_particles_total
*   max_first_species_id
*   total_new_stars
*   current_id
*   first_star
*   first_star_id
*   particle_species(
*   particle_id_is_star
*   particle_alloc(

#Files

*   cosmology/cosmics/cosmics_init.c [dhr - converted to particleid_t to avoid overflow]
*   cosmology/gic/gic_init.c [already 64-bit capable if sizeof(long) == sizeof(int64_t)]
*   src/core/agn.c [fixed_comparison function, moved to particle.c]
*   src/core/cache.c [fixed comparison function]
*   src/core/density.c [done]
*   src/core/halo_finder.c [done]
*   src/core/hydro_tracer.c [have questions]
*   src/core/io_artio.c [l.312?, more questions?]
*   src/core/io_cart.c [have questions]
*   src/core/io_cart1.def [have questions]
*   src/core/io_cart2.def [have questions]
*   src/core/particle.c [done]
*   src/core/particle.h [done]
*   src/core/particle_buffer.c [nothing found, dhr agrees]
*   src/core/tree_debug.c [MPI_INT?]
*   src/drivers/start_analysis.c [done]
*   src/extra/gic_tools.c [no changes found]
*   src/extra/halo_tools.c [done]
*   src/extra/ifrit.c [no changes found]
*   src/extra/sps.c [no changes necessary]
*   src/run/agn_step.c [changed qsort comparison function]
*   src/run/logging.c [no changes found]
*   src/run/particle_step.c [no changes necessary]
*   src/run/starformation_feedback_step.c [no changes found]
*   src/run/starformation_step.c [fixed star id remapping code]
*   src/run/step.c [formatting changes only]
*   tests/disk/disk.c [fixed particle id variables when calling particle_alloc]
*   tests/feedback/start_radp.c [no changes necessary]
*   tests/io/start_io.c [changed particle id type, left other 64-bit unsafe code since test is serial]
*   tests/orbit/start_orbit.c [no changes since intended for 2 particles]
*   tests/rt/et/start_rt.c [changed particle id type, left other 64-bit unsafe code]
*   tests/rt/iliev1a/start_rt.c [changed particle id type, left other 64-bit unsafe code]
*   tests/zeldovich/start_zeldovich.c [changed particle id calculation, num_particles_total made 64-bit safe]
*   utils/cart_file_utilities/convert_particle_endianness.c [questions have pages]
*   utils/cart_file_utilities/print_particle_aexp.c
*   utils/cart_file_utilities/print_particle_info.c

