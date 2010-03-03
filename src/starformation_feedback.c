#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"


typedef double(*fimf)(double);

double f_IMF_Salpeter(double m);
double f_IMF_MillerScalo(double m);
double f_IMF_Chabrier(double m);


struct InitialMassFunction
{
  char* name;
  fimf  f;
}
const IMF_fun[] = { { "Salpeter", f_IMF_Salpeter }, { "Miller-Scalo", f_IMF_MillerScalo }, { "Chabrier", f_IMF_Chabrier } };


struct
{
  int type;              /* type of IMF function */
  double slope;          /* used to be called a_IMF */
  double min_mass;       /* used to be called aM_stl */
  double max_mass;       /* used to be called aM_stu */
  double min_SNII_mass;  /* used to be called aM_SNII */
  double min_SNIa_mass;  /* used to be called aM_SNIa1 */
  double max_SNIa_mass;  /* used to be called aM_SNIa2 */
}
imf = { 1, 2.35, 0.1, 100.0, 8.0, 3.0, 8.0 };


void control_parameter_set_imf(const char *value, void *ptr, int ind)
{
  const int n = sizeof(IMF_fun)/sizeof(struct InitialMassFunction);

  int i;
  char *str = (char *)ptr;

  imf.type = -1;

  for(i=0; i<n; i++)
    {
      if(strcmp(str,IMF_fun[i].name) == 0)
	{
	  imf.type = i;
	}
    }

  if(imf.type < 0)
    {
      cart_debug("String '%s' is not a valid name of IMF. Valid names are:",value);
      for(i=0; i<n; i++) cart_debug("%s",IMF_fun[i].name);
      cart_error("ART is terminating.");
    }
}


void control_parameter_list_imf(FILE *stream, const void *ptr)
{
  control_parameter_list_string(stream,IMF_fun[imf.type].name);
}


#ifdef FEEDBACK
struct
{
  double energy_per_explosion;     /* used to be called E_51 */
  double time_delay;               /* used to be called t_fb */
}
snII = { 2.0, 1.0e3 };


#ifdef FEEDBACK_SNIa
struct
{
  double energy_per_explosion;          /* used to be called E_51 */
  double time_delay;                    /* used to be called t_SNIa */
  double exploding_fraction;            /* used to be called C_SNIa */
  double mass_in_metals_per_supernova;  /* used to be called ejM_SNIa */
}
snIa = { 2.0, 2.0e8, 1.5e-2, 1.3 };

void control_parameter_set_tSNIa(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  /*
  //  Backward compatibility
  */
  if(ind == 2) snIa.time_delay *= 1.0e9;
}
#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */


#ifdef STELLARMASSLOSS
struct
{
  double lost_fraction;   /* used to be called c0_ml */
  double time_interval;   /* used to be called T0_ml */
}
ml = { 0.05, 5.0e6 };

void control_parameter_set_t0ml(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  /*
  //  Backward compatibility
  */
  if(ind == 2) ml.time_interval *= 1.0e6;
}
#endif /* STELLARMASSLOSS */


double feedback_temperature_ceiling = 1.0e8;  /* Used to be called T_max_feedback; also, was a define in HART */


void config_init_star_formation_feedback()
{
  ControlParameterOps control_parameter_imf = { control_parameter_set_imf, control_parameter_list_imf };
#ifdef FEEDBACK
#ifdef FEEDBACK_SNIa
  ControlParameterOps control_parameter_tSNIa = { control_parameter_set_tSNIa, control_parameter_list_double };
#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */
#ifdef STELLARMASSLOSS
  ControlParameterOps control_parameter_t0ml = { control_parameter_set_t0ml, control_parameter_list_double };
#endif /* STELLARMASSLOSS */


  /*
  //  IMF
  */
  control_parameter_add2(control_parameter_imf,&imf.type,"imf","imf.type","the name of the IMF. Valid names are 'Salpeter', 'Miller-Scalo', and 'Chabrier'.");

  control_parameter_add3(control_parameter_double,&imf.slope,"imf:slope","imf.slope","a_imf","the slope of the power-law IMF; set it to zero for Miller-Scalo IMF.");

  control_parameter_add3(control_parameter_double,&imf.min_mass,"imf:min-mass","imf.min_mass","am_stl","the minimum stellar mass in the IMF model.");

  control_parameter_add3(control_parameter_double,&imf.max_mass,"imf:max-mass","imf.max_mass","am_stu","the maximum stellar mass in the IMF model.");

  control_parameter_add3(control_parameter_double,&imf.min_SNII_mass,"imf:min-SNII-mass","imf.min_SNII_mass","am_snii","the minimum mass of stars that explode as type II supernovae.");

  control_parameter_add3(control_parameter_double,&imf.min_SNIa_mass,"imf:min-SNIa-mass","imf.min_SNIa_mass","am_snia1","the minimum mass of stars that explode as type Ia supernovae.");

  control_parameter_add3(control_parameter_double,&imf.max_SNIa_mass,"imf:max-SNIa-mass","imf.max_SNIa_mass","am_snia2","the maximum mass of stars that explode as type Ia supernovae.");

#ifdef FEEDBACK
  /*
  //  type II supernova feedback
  */
  control_parameter_add3(control_parameter_double,&snII.energy_per_explosion,"snII:energy-per-explosion","snII.energy_per_explosion","e_51","average energy per type II supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_double,&snII.time_delay,"snII:time_delay","snII.time_delay","t_fb","time delay (in yrs) between the formation of the stellar particle and type II supernova explosions.");

#ifdef FEEDBACK_SNIa
  /*
  //  type Ia supernova feedback
  */
  control_parameter_add3(control_parameter_double,&snIa.energy_per_explosion,"snIa:energy-per-explosion","snIa.energy_per_explosion","e_51","average energy per type Ia supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_tSNIa,&snIa.time_delay,"snIa:time-delay","snIa.time_delay","t_snia","average time delay (in yrs) between the formation of the stellar particle and type Ia supernova explosions.");

  control_parameter_add3(control_parameter_double,&snIa.exploding_fraction,"snIa:exploding-fraction","snIa.exploding_fraction","c_snia","fraction of stars exploding as type Ia supernovae.");

  control_parameter_add3(control_parameter_double,&snIa.mass_in_metals_per_supernova,"snIa:mass-in-metals-per-supernova","snIa.mass_in_metals_per_supernova","ejm_snia","average mass (in solar masses) in metals ejected per type Ia supernova explosion.");
#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */

#ifdef STELLARMASSLOSS
  /*
  //  mass loss
  */
  control_parameter_add3(control_parameter_double,&ml.lost_fraction,"ml:lost-fraction","ml.lost_fraction","c0_ml","fraction of stellar mass ejected back into the ISM by stellar mass loss.");

  control_parameter_add3(control_parameter_t0ml,&ml.time_interval,"ml:time-interval","ml.time_interval","t0_ml","time interval (in yrs) over which the stellar mass loss occurs.");
#endif /* STELLARMASSLOSS */

  /*
  //  other
  */
  control_parameter_add3(control_parameter_double,&feedback_temperature_ceiling,"fb:temperature-ceiling","feedback_temperature_ceiling","T_max_feedback","maximum gas temperature for the feedback to operate. No feedback is allowed in the gas with the temperature above this limit.");
}


void config_verify_star_formation_feedback()
{
  /*
  //  IMF
  */
  cart_assert(imf.type>=0 && imf.type<2);

  cart_assert(imf.slope > 0.0);

  cart_assert(imf.min_mass > 0.0);

  cart_assert(imf.max_mass > 0.0);

  cart_assert(imf.min_SNII_mass > 1.0);

  cart_assert(imf.min_SNIa_mass > 1.0);

  cart_assert(imf.max_SNIa_mass > imf.min_SNIa_mass);

#ifdef FEEDBACK
  /*
  //  type II supernova feedback
  */
  cart_assert(!(snII.energy_per_explosion < 0.0));

  cart_assert(snII.time_delay > 0.0);

#ifdef FEEDBACK_SNIa
  /*
  //  type Ia supernova feedback
  */
  cart_assert(!(snIa.energy_per_explosion < 0.0));

  cart_assert(snIa.time_delay > 0.0);

  cart_assert(snIa.exploding_fraction>0.0 && snIa.exploding_fraction<1.0);

  cart_assert(snIa.mass_in_metals_per_supernova>0.0 && snIa.mass_in_metals_per_supernova<imf.max_SNIa_mass);
#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */

#ifdef STELLARMASSLOSS
  /*
  //  mass loss
  */
  cart_assert(ml.lost_fraction>0.0 && ml.lost_fraction<1.0);

  cart_assert(ml.time_interval > 0.0);
#endif /* STELLARMASSLOSS */

  /*
  //  other
  */
  cart_assert(feedback_temperature_ceiling > 1.0e6);
}


typedef struct
{
  double energy;
  double metals;
  double dt;
}
fb_pars;

#ifdef FEEDBACK
fb_pars fbp_snII_phys, fbp_snII_code;
#ifdef FEEDBACK_SNIa
fb_pars fbp_snIa_phys, fbp_snIa_code;
#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */


double f_IMF_Salpeter( double m )
{
  return pow( m, -imf.slope );
}


double f_IMF_MillerScalo( double m )
{
  /* Miller-Scalo (1979, ApJS 41, 513, eq. 30, Table 7) IMF */
  const double C_1 =  1.09;
  const double C_2 = -1.02;
  return exp( -C_1 * ( (log10(m) - C_2)*(log10(m) - C_2) ) ) / m;
}


double f_IMF_Chabrier( double m )
{
  /* Chabrier, G. (2001, ApJ 554, 1274) */
  const double m0_Ch = 716.4;
  const double beta_Ch = 0.25;
  const double alpha_Ch = -3.3;
  return exp( -pow( m0_Ch/m, beta_Ch) ) * pow(m,alpha_Ch);
}


/* Normalized */
double f_SNIa( double xd )
{
  return exp(-xd*xd) * sqrt(xd*xd*xd) / 1.812804954;
}


double f_IMF( double amstar )
{
  return IMF_fun[imf.type].f(amstar);
}


double fm_IMF( double amstar )
{
  return amstar * f_IMF(amstar);
}


double fmet_ej( double amstar )
{
  return min( 0.2, max( 0.01*amstar - 0.06, 1e-20 ) );
}


double fej_IMF( double amstar )
{
  return amstar * f_IMF(amstar) * fmet_ej(amstar);
}


void init_star_formation_feedback()
{
  const double erg = cgs->g*pow(cgs->cm/cgs->s,2.0);
  double total_mass;
  double number_SNII, number_SNIa;

#ifdef FEEDBACK
  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( fm_IMF, imf.min_mass, imf.max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  number_SNII = integrate( f_IMF, imf.min_SNII_mass, imf.max_mass, 1e-6, 1e-9 );
  cart_assert(number_SNII > 0.0);

#ifdef ENRICH 
  fbp_snII_phys.metals = integrate( fej_IMF, imf.min_SNII_mass, imf.max_mass, 1e-6, 1e-9 )/total_mass;
#endif /* ENRICH */

#ifdef FEEDBACK_SNIa

  number_SNIa = snIa.exploding_fraction*integrate( f_IMF, imf.min_SNIa_mass, imf.max_SNIa_mass, 1e-6, 1e-9 );
  cart_assert(number_SNIa > 0.0);

  fbp_snII_phys.dt = snII.time_delay;
  fbp_snII_phys.energy = 1e51*erg*snII.energy_per_explosion*number_SNII/(constants->Msun*total_mass);

  fbp_snIa_phys.dt = snIa.time_delay;
  fbp_snIa_phys.energy = 1e51*erg*snIa.energy_per_explosion*number_SNIa/(constants->Msun*total_mass);

#ifdef ENRICH_SNIa 
  fbp_snIa_phys.metals = snIa.mass_in_metals_per_supernova*number_SNIa/total_mass;
#endif /* ENRICH_SNIa */

#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */

  /*
  // The rest is for diagnostic only
  */
  if(local_proc_id == MASTER_NODE)
    {
#ifdef FEEDBACK

      cart_debug("Number of SNII explosions per unit mass: %le per Msun",number_SNII/total_mass);
      cart_debug("SNII specific energy: %le erg/g = (%le km/s)^2",fbp_snII_phys.energy,sqrt(fbp_snII_phys.energy)/constants->kms);

#ifdef ENRICH 
	cart_debug("SNII metal fraction : %le",fbp_snII_phys.metals);
#endif /* ENRICH */

#ifdef FEEDBACK_SNIa

      cart_debug("Number of SNIa explosions per unit mass: %le per Msun",number_SNIa/total_mass);
      cart_debug("SNIa specific energy: %le erg/g = (%le km/s)^2",fbp_snIa_phys.energy,sqrt(fbp_snIa_phys.energy)/constants->kms);

#ifdef ENRICH_SNIa 
	cart_debug("SNIa metal fraction : %le",fbp_snIa_phys.metals);
#endif /* ENRICH_SNIa */

#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */
    }
}


#if defined(HYDRO) && defined(PARTICLES)

double dUfact;  /* must be here to simplify OpenMP directives */

#ifdef STELLARMASSLOSS
double dt_ml_code;   /* used to be called T0_ml_code */
#endif /* STELLARMASSLOSS */


void stellar_feedback(int level, int cell, int ipart, double delta_t, double t_next, double vx, double vy, double vz)
{
  double dteff, phi, dU;
  double dmloss, rhor, e_old, rhofact;

  /* do feedback, enrichment, etc */
#ifdef FEEDBACK
  dteff = particle_t[ipart] - star_tbirth[ipart];
  if(dteff < fbp_snII_code.dt)
    {
      phi = min(delta_t,fbp_snII_code.dt-dteff)/fbp_snII_code.dt;

#ifdef ENRICH
      cell_gas_metal_density_II(cell) += phi*fbp_snII_code.metals*star_initial_mass[ipart];
#endif /* ENRICH */

      dU = min(phi*fbp_snII_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));

      /* limit energy release and don't allow to explode in hot bubble */
      if ( units->temperature*cell_gas_temperature(cell) < feedback_temperature_ceiling )
	{
	  cell_gas_energy(cell) += dU;
	  cell_gas_internal_energy(cell) += dU;
	  cell_gas_pressure(cell) += dU*(cell_gas_gamma(cell)-1);
	}
    }
#endif /* FEEDBACK */

#ifdef FEEDBACK_SNIa
  dteff = t_next - star_tbirth[ipart];
  if(dteff > 0.1*fbp_snIa_code.dt)
    {
      phi = f_SNIa(fbp_snIa_code.dt/dteff)*(delta_t/fbp_snIa_code.dt);

#ifdef ENRICH_SNIa
      cell_gas_metal_density_Ia(cell) += phi*fbp_snIa_code.metals*star_initial_mass[ipart];
#endif /* ENRICH_SNIa */

      dU = min(phi*fbp_snIa_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));

      /* limit energy release and don't allow to explode in hot bubble */
      if ( units->temperature*cell_gas_temperature(cell) < feedback_temperature_ceiling )
	{
	  cell_gas_energy(cell) += dU;
	  cell_gas_internal_energy(cell) += dU;
	  cell_gas_pressure(cell) += dU*(cell_gas_gamma(cell)-1);
	}
    }
#endif /* FEEDBACK_SNIa */
	
#ifdef STELLARMASSLOSS
  /* limit mass loss to 10% of star's current mass */
  dmloss = min( 0.1*particle_mass[ipart],
		star_initial_mass[ipart]*delta_t*ml.lost_fraction / 
		(particle_t[ipart] - star_tbirth[ipart] + dt_ml_code) );
					
  particle_mass[ipart] -= dmloss;

  /* convert to density for cell values */
  dmloss *= cell_volume_inverse[level];

  /* account for momentum change */
  rhor = 1.0 / cell_gas_density(cell);
  e_old = cell_gas_energy(cell) -
    0.5 * ( cell_momentum(cell,0)*cell_momentum(cell,0) +
	    cell_momentum(cell,1)*cell_momentum(cell,1) +
	    cell_momentum(cell,2)*cell_momentum(cell,2) ) * rhor; 
  cell_gas_density(cell) += dmloss;
  rhofact = rhor * cell_gas_density(cell);
  
  cell_momentum(cell,0) += dmloss * vx;
  cell_momentum(cell,1) += dmloss * vy;
  cell_momentum(cell,2) += dmloss * vz;
			
  cell_gas_energy(cell) = e_old + 
    0.5 * ( cell_momentum(cell,0)*cell_momentum(cell,0) +
	    cell_momentum(cell,1)*cell_momentum(cell,1) +
	    cell_momentum(cell,2)*cell_momentum(cell,2) ) /
    cell_gas_density(cell);
  cell_gas_internal_energy(cell) *= rhofact;
  cell_gas_pressure(cell) *= rhofact;

#ifdef ENRICH
  cell_gas_metal_density_II(cell) += dmloss*star_metallicity_II[ipart];
#ifdef ENRICH_SNIa
  cell_gas_metal_density_Ia(cell) += dmloss*star_metallicity_Ia[ipart];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
#endif /* STELLARMASSLOSS */
}


void setup_star_formation_feedback(int level)
{
  dUfact = feedback_temperature_ceiling/(units->temperature*constants->wmu*(constants->gamma-1));

#ifdef FEEDBACK

  fbp_snII_code.dt = fbp_snII_phys.dt*constants->yr/units->time;
  fbp_snII_code.energy = fbp_snII_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  fbp_snII_code.metals = fbp_snII_phys.metals*cell_volume_inverse[level]; 

#ifdef FEEDBACK_SNIa

  fbp_snIa_code.dt = fbp_snIa_phys.dt*constants->yr/units->time;
  fbp_snIa_code.energy = fbp_snIa_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  fbp_snIa_code.metals = fbp_snIa_phys.metals*cell_volume_inverse[level]; 

#endif /* FEEDBACK_SNIa */
#endif /* FEEDBACK */

#ifdef STELLARMASSLOSS
  dt_ml_code = ml.time_interval*constants->yr/units->time;
#endif /* STELLARMASSLOSS */
}

#endif /* HYDRO && PARTICLES */
#endif /* STARFORM && FEEDBACK */
