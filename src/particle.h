#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef PARTICLES 

#define NULL_PARTICLE   (-1)
#define FREE_PARTICLE_LEVEL	(-1)

#define SAVE_PARTICLE_LISTS     0
#define FREE_PARTICLE_LISTS     1

#define MAX_PARTICLE_SPECIES	10

extern int num_row;

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
extern double ekin2;
extern double au0;
extern double aeu0;
extern double ap0;
extern double ap1;

/* particle species */
extern int num_particle_species;
extern float particle_species_mass[MAX_PARTICLE_SPECIES];
extern int particle_species_num[MAX_PARTICLE_SPECIES];
extern int particle_species_indices[MAX_PARTICLE_SPECIES+1];

extern int particle_level[/* num_particles */];
extern float particle_mass[/* num_particles */];
extern int particle_id[/* num_particles */];
extern int particle_list_next[/* num_particles */];
extern int particle_list_prev[/* num_particles */];

/* extern int cell_particle_list[num_cells]; */
extern int CELL_ARRAY(cell_particle_list);

extern int num_local_particles;
extern long num_particles_total;
extern int next_free_particle;
extern int free_particle_list;
extern int particle_list_enabled;

#ifdef STARFORM
extern int next_free_star_particle;
extern int free_star_particle_list;
#endif /* STARFORM */

int particle_alloc( int id );
void particle_free( int ipart );
void particle_list_free( int ihead );
void particle_move( int ipart_old, int ipart_new );

void init_particles();
void build_mesh();

void update_particle_list( int level );
void trade_particle_lists( int num_parts_to_send[MAX_PROCS], int *particle_list_to_send[MAX_PROCS], int trade_level, int free_particle_flag );
void build_particle_list();
void split_particle_list( int cell );
void join_particle_list( int cell );
void insert_particle( int cell, int part );
void delete_particle( int cell, int part );
void rebuild_particle_list();
int particle_specie( int id );

#ifdef STARFORM
#define particle_id_is_star(id)		(id >= particle_species_indices[num_particle_species-1])
#define particle_is_star(index)		(particle_id_is_star( particle_id[index] ))
#endif /* STARFORM */

#endif /* PARTICLES */

#endif
