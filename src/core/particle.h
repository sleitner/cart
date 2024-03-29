#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef PARTICLES

/*
//  C++ compatibility: INT64_MAX is defined only if __STDC_LIMIT_MACROS set explicitly (ISO C99 standard)
*/
#ifdef __cplusplus
#define __STDC_LIMIT_MACROS
#endif

#include <stdint.h>
#include <limits.h>

#ifdef OLDSTYLE_32BIT_PARTICLEID
#define particleid_t        int
#define MPI_PARTICLEID_T    MPI_INT
#define NULL_PARTICLE       (-1)
#define PARTICLEID_MAX      INT_MAX
#else
#define particleid_t        int64_t
#define MPI_PARTICLEID_T    MPI_LONG
#define NULL_PARTICLE       (-1L)
#define PARTICLEID_MAX      INT64_MAX
#endif

#define FREE_PARTICLE_LEVEL	(-1)

#define SAVE_PARTICLE_LISTS     0
#define FREE_PARTICLE_LISTS     1

#define MAX_PARTICLE_SPECIES	10

extern double particle_t[/* num_particles */];
extern double particle_dt[/* num_particles */];
extern double particle_x[/* num_particles */][nDim];
extern double particle_v[/* num_particles */][nDim];

#ifdef GRAVITY
extern float particle_pot[/* num_particles */];
#endif /* GRAVITY */

/* variables for monitoring energy */
extern double tintg;
extern double ekin;
extern double ekin1;
extern double au0;
extern double aeu0;
extern double ap0;
extern double ap1;

/* particle species */
extern int num_particle_species;
extern float particle_species_mass[MAX_PARTICLE_SPECIES];
extern particleid_t particle_species_num[MAX_PARTICLE_SPECIES];
extern particleid_t particle_species_indices[MAX_PARTICLE_SPECIES+1];

extern int particle_level[/* num_particles */];
extern float particle_mass[/* num_particles */];
extern particleid_t particle_id[/* num_particles */];
extern int particle_list_next[/* num_particles */];
extern int particle_list_prev[/* num_particles */];

extern int cell_particle_list[num_cells];

extern int num_local_particles;
extern particleid_t num_particles_total;
extern int next_free_particle;
extern int free_particle_list;
extern int particle_list_enabled;

#ifdef STAR_FORMATION
extern int next_free_star_particle;
extern int free_star_particle_list;
#endif /* STAR_FORMATION */

int particle_alloc( particleid_t id );
void particle_free( int ipart );
void particle_list_free( int ihead );
void particle_move( int ipart_old, int ipart_new );

/* qsort comparison functions */
int compare_particle_species_id( const void *a, const void *b );
int compare_particle_ids( const void *a, const void *b );
int compare_particle_mass( const void *a, const void *b );

void init_particles();
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
void build_mesh();
#ifdef REFINEMENT
void get_refinement_region();
void build_refinement_region(int do_load_balance);
#endif /* REFINEMENT */
#endif /* GRAVITY || RADIATIVE_TRANSFER */

void update_particle_list( int level );
void trade_particle_lists( int num_parts_to_send[MAX_PROCS], int *particle_list_to_send[MAX_PROCS], int trade_level, int free_particle_flag );
void build_particle_list();
void split_particle_list( int cell );
void join_particle_list( int cell );
void insert_particle( int cell, int part );
void delete_particle( int cell, int part );
void rebuild_particle_list();
int particle_species( particleid_t id );

#ifdef STAR_FORMATION
#define particle_id_is_star(id)		(id >= particle_species_indices[num_particle_species-1])
#define particle_is_star(index)		(particle_id_is_star( particle_id[index] ))
#endif /* STAR_FORMATION */

#endif /* PARTICLES */

#endif
