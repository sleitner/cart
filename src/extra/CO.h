#ifndef __CO_H__
#define __CO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#if defined(HYDRO) && defined(CO)

#include "control_parameter.h"

/*******************************************************************/
/* the CO & X-factor model is described in detail in               */
/*                                                                 */
/* Feldmann, Gnedin & Kravtsov APJ, 747, 124 (2012)                */
/*                                                                 */
/* please cite this reference, when using this module              */
/*                                                                 */
/* R. Feldmann 2011/12                                             */
/*******************************************************************/

typedef double (*Xfactor_fp)(int cell, float UMW, int level);
extern Xfactor_fp xfactor_cell;
extern int xfactor_method;
extern double T_excitation;
extern double z_CMB;
extern double const_DV;
extern double min_DV;
extern double f_grav;

typedef struct {
  float *data;
  unsigned *dims;
  unsigned num_dims;
} COTable;

struct {
   COTable AV_tab;
   COTable AVL_tab;
   COTable L_tab;
   COTable Z_tab;
   COTable UMW_tab;
      
   COTable xCO_tab;
   COTable xH2_tab; 
   
   COTable tau_tab;
   COTable escapefraction_tab;  
    
   char tablePath[CONTROL_PARAMETER_STRING_LENGTH];
} xfactor_method0;

double get_xCO(double AV, double Lcm,  double Z, double UMW);
double get_xH2(double AV, double Lcm,  double Z, double UMW);
double get_escape_fraction(double tau10);

void NH2_level(int level, int num_level_cells, int *level_cells, float *var_level);
void xfactor_level( int level, int num_level_cells, int *level_cells, float *var_level);
void xfactor_upto_level(int top_level, float *xf);
void xfactor_upto_level_hierarchy(int top_level, float *var);

double get_NH2(int cell, float UMW, int level);
double xfactor_cell_0(int cell, float UMW, int level);
double xfactor_cell_1(int cell, float UMW, int level);
double xfactor_cell_2(int cell, float UMW, int level);
double xfactor_cell_3(int cell, float UMW, int level);
double xfactor_cell_4(int cell, float UMW, int level);

void config_init_CO();
void config_verify_CO();
void setup_xfactor_methods(int level);
void unset_xfactor_methods(int level);

void CO_load_table(const char *filename, COTable * tab, int readln);
void CO_init_table(COTable * tab);
void CO_destroy_table(COTable * tab);

#endif /* HYDRO && CO */

#endif
