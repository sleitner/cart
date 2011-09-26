#define ONE_CELL_IN_THE_CENTER


#define refine_radius	(8.0)
int const slice_axis_z=2;
#define slice_hsize_pc  50

#ifdef COSMOLOGY
#define omm0 1.0
#define oml0 0.0
#define omb0 1.0
#define hubble 1.0
#define deltadc 0.0
#define a0 0.9
#define boxh (1.0e-3/a0*hubble) //1pkpc //for level4
#define rho0            (1e7)
#define blast_radius    (cell_size[max_level]) 

#ifdef STARFORM

#define E_ambient       (1.0e10)
const double mstar_one_msun = 10;
extern int last_star_id;
int last_star_id=-1;



#else /* !STARFORM */
#define E_ambient       (1.0)
#define E_det           (1.0e7*E_ambient)
#define P_ambient	(1.0) 
#endif /* STARFORM */


#else /* !COSMOLOGY */

#define blast_radius    (cell_size[max_level]) 
#define rho0            (1.0e7) 
#define P_ambient	(1.0) 
#define refine_radius	(4.0)

#endif  /* COSMOLOGY */
