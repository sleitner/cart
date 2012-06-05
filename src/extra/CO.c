#include "config.h"
#if defined(HYDRO) && defined(CO)

/*******************************************************************/
/* the CO & X-factor model is described in detail in               */
/*                                                                 */
/* Feldmann, Gnedin & Kravtsov APJ, 747, 124 (2012)                */
/*                                                                 */
/* please cite this reference, when using this module              */
/*                                                                 */
/* R. Feldmann 2011/12                                             */
/*******************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <dirent.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "iterators.h"
#include "timing.h"
#include "parallel.h"
#include "hydro.h"
#include "particle.h"
#include "starformation.h"

#include "rt.h"

#include "CO.h"
#include "interpolation.h"

#include "tree.h"
#include "units.h"

#include "times.h"

#ifdef RADIATIVE_TRANSFER
#include "frt/frt_c.h"
#endif

#define READLN 1
#define READNORMAL 0

/* #define H2GLOVER */

#define XLSOB
#define FIXEDLVALUE 20
#define CLUMPING

Xfactor_fp xfactor_cell = NULL;
int xfactor_method = 2;
double T_excitation = 10;
double z_CMB = 0;
double f_grav = 1;
double min_DV = 1;   /* km/s */
double const_DV = 3; /* km/s */

/* compute X-factor at given level; taking into account this level and all higher res. levels */
void xfactor_upto_level(int top_level, float *xf) {
  MESH_RUN_DECLARE(level,cell);
  float *var;
  
  var = cart_alloc(float,num_cells);
  memset(var,0,sizeof(float)*num_cells);
  
  xfactor_upto_level_hierarchy(top_level, var);
    
  MESH_RUN_OVER_LEVELS_BEGIN(level,top_level,top_level);
  
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,var,cell_child_oct,cell_vars,level,xf)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

  xf[_Index]=var[cell];
  
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  MESH_RUN_OVER_LEVELS_END;
  
  cart_free(var);
}

/* var: allocated of size num_cells & filled with zeros */
void xfactor_upto_level_hierarchy(int top_level, float *var) {
  
  MESH_RUN_DECLARE(level,cell);
  int j,c;
  float *var_level;
  double aux1, aux2;
  
  float soblen;
  
  MESH_RUN_OVER_LEVELS_BEGIN(level,_MaxLevel,top_level);

  var_level = cart_alloc(float,_Num_level_cells);
  
  xfactor_level(level,_Num_level_cells,_Level_cells,var_level);
  
#pragma omp parallel for default(none), private(_Index,cell,j,aux1,aux2,soblen,c), shared(_Num_level_cells,_Level_cells,var_level,var,cell_child_oct,cell_vars,level, units,constants)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

  if(cell_is_leaf(cell))
    {
      var[cell] = MAX(0.0,var_level[_Index]);
    }
  else
    {
      aux1=0.; 
      aux2=0.;
      
#ifdef XLSOB
      for(j=0; j<num_children; j++) {
          c = cell_child(cell,j);
	  
	  soblen = cell_sobolev_length2(c,cell_level(c),NULL);
	  
	  aux1 += cell_H2_density(c)*soblen;
	  aux2 += cell_H2_density(c)*soblen/var[c];
      }
#else
      for(j=0; j<num_children; j++) {
	  aux1 += cell_H2_density(cell_child(cell,j));
	  aux2 += cell_H2_density(cell_child(cell,j))/var[cell_child(cell,j)];
      }
#endif 
     
      if (aux2==0)
        var[cell]= FLT_MAX;
      else	  
        var[cell] = aux1/aux2;
    }
  cart_assert( var[cell]>=0 && !(var[cell]!=var[cell]) );
  
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  
  cart_free(var_level);

  MESH_RUN_OVER_LEVELS_END;
}


void xfactor_level(int level, int num_level_cells, int *level_cells, float *var_level)
{
  int i;
  int cell;
  
  setup_xfactor_methods(level);

  fprintf(stderr,"Level %d - Tables loaded: num_cells=%d\n",level,num_level_cells);

  for(i=0; i<num_level_cells; i++)
    {
      cell = level_cells[i];
      var_level[i] = -1.0;
      if(cell_is_leaf(cell))
	  var_level[i] = xfactor_cell(cell,rtUmw(cell),level);     
    }
    
  unset_xfactor_methods(level);  
    
}

void NH2_level(int level, int num_level_cells, int *level_cells, float *var_level)
{
  int i;
  int cell;

  for(i=0; i<num_level_cells; i++)
    {
      cell = level_cells[i];
      var_level[i] = -1.0;
      if(cell_is_leaf(cell))
	  var_level[i] = get_NH2(cell,rtUmw(cell),level);
    }
}

double get_xCO(double AV, double Lcm,  double Z, double UMW) {
    
    double L = Lcm / 3.08568025e18; /* L in pc */
    
    /* interpolate in log space */
    double result = inter_linear_4d( log(UMW), log(Z), log(L), log(AV), 
             xfactor_method0.UMW_tab.data, xfactor_method0.UMW_tab.dims[0],
             xfactor_method0.Z_tab.data, xfactor_method0.Z_tab.dims[0],
	     xfactor_method0.L_tab.data, xfactor_method0.L_tab.dims[0],
             xfactor_method0.AV_tab.data, xfactor_method0.AV_tab.dims[0],
             xfactor_method0.xCO_tab.data);

    result = exp(result);

    if (result > 1.41e-4*Z)  
      result = 1.41e-4*Z;

    return result; /* this is now the mass-weighted CO abundance, based on simulations by Glover & Mac Low 2010 */
}

double get_xH2(double AV, double Lcm,  double Z, double UMW) {
   
    double L = Lcm / 3.08568025e18; /* L in pc */
   
    double result = inter_linear_2d( log(Z), log(AV/L),              
             xfactor_method0.Z_tab.data, xfactor_method0.Z_tab.dims[0],	     
             xfactor_method0.AVL_tab.data, xfactor_method0.AVL_tab.dims[0],
             xfactor_method0.xH2_tab.data);
    return exp(result); /* this is the mass-weighted H2 abundance, based on simulations by Glover & Mac Low 2010 */
}

double get_escape_fraction(double tau10) {
  
  double result = inter_linear_1d( log(tau10), 
                            xfactor_method0.tau_tab.data, xfactor_method0.tau_tab.dims[0],	 
			    xfactor_method0.escapefraction_tab.data);
  return exp(result);
}

double get_NH2(int cell, float UMW, int level) {
  double Z = rtDmw2(cell); /* this is gas-to-dust ratio in solar units */
  double n0 = units->number_density*cell_gas_density(cell)*constants->XH; /* hydrogen atom density in cm^(-3) */
  double nH2 = units->number_density*cell_H2_density(cell); /* molecular hydrogen density in cm^(-3) */
  double NH2, Lcm;
  float soblen, length;
  
#ifdef XLSOB  
  soblen = cell_sobolev_length(cell);
#else
  soblen = cell_size[cell_level(cell)];
#endif
 
#ifdef CLUMPING
  length = FIXEDLVALUE*constants->pc/units->length;
  n0 *= soblen/length;   
  nH2 *= soblen/length;
#else
  length = soblen; 
#endif

  Lcm = length*units->length;

#ifdef H2GLOVER
  NH2 = get_xH2(AV,Lcm,Z,UMW)/2 * n0 * Lcm; /* in cm^-2 */
#else
  NH2 = nH2*Lcm;
#endif
  
  return NH2;
}

/* constant CO line width */
double xfactor_cell_0(int cell, float UMW, int level)
{
  double Z = rtDmw2(cell); /* this is gas-to-dust ratio in solar units */  
  double n0 = units->number_density*cell_gas_density(cell)*constants->XH; /* hydrogen atom density in cm^(-3) */ 
  double nH2 = units->number_density*cell_H2_density(cell); /* molecular hydrogen density in cm^(-3) */
  double AV, Lcm;
  double TB, Tbg;
  double xCO, NCO, NH2;
  double prefac10, T10, f0, tau10, WCO;
  float soblen, length;
  double DV;

/* the variable soblen measures the full width of the cell
 * the variable length measures the width of the squashed cell */

#ifdef XLSOB  
  soblen=cell_sobolev_length(cell);
#else
  soblen = cell_size[cell_level(cell)];
#endif
 
#ifdef CLUMPING
  length = FIXEDLVALUE*constants->pc/units->length;
  n0 *= soblen/length;   
  nH2 *= soblen/length;
#else
  length = soblen; 
#endif

  Lcm = length*units->length;
  
  AV = 5.348e-22*Z*n0*Lcm;   
  
  xCO = get_xCO(AV,Lcm,Z,UMW);
  NCO = xCO * n0 * Lcm;   /* in cm^-2 */
 
#ifdef H2GLOVER
  NH2 = get_xH2(AV,Lcm,Z,UMW)/2 * n0 * Lcm; /* in cm^-2 */
#else
  NH2 = nH2*Lcm;
#endif
 
  /* these constants are taken from formula (2) in Glover and Mac Low 2010
     and with values A10, nu_10 from the Schoeier et al. 2005 database */
  prefac10 = 1.5125e-15; /* cm^2 km s^-1 */
  T10 = 5.532; /* K */
  
  if (z_CMB>=-1e-6) {
    if (z_CMB<0)
      z_CMB=0;
    Tbg = 2.725*(1+z_CMB);
  } else {
    Tbg = 2.725/auni[min_level];
  }
 
  TB = T10*(1./(exp(T10/T_excitation)-1)-1./(exp(T10/Tbg)-1));			
 
  DV = const_DV;
 
  f0 = T10 / (2*T_excitation); /* this only holds in LTE, and for Tex>>T0; i.e. for Tex>=10 K */
  tau10 = prefac10 * (1 - exp( - T10/T_excitation ) ) * f0 * NCO / DV;

  WCO = TB * DV * get_escape_fraction(tau10);    		    
  
  return NH2/WCO; 
}

/* not implemented */  
double xfactor_cell_1(int cell, float UMW, int level)
{ 
  return 0;
}

/* like method 0, but rescale CO line width according to the virial theorem */
/* default method */
double xfactor_cell_2(int cell, float UMW, int level) { 
  double Z = rtDmw2(cell); /* this is baryon fraction in dust relative to the solar neighbourhood; usually just Z/Zsun */  
  double n0 = units->number_density*cell_gas_density(cell)*constants->XH; /* hydrogen atom density in cm^(-3) */ 
  double nH2 = units->number_density*cell_H2_density(cell); /* molecular hydrogen density in cm^(-3) */
  double AV, Lcm, L;
  double TB, Tbg;
  double xCO, NCO, NH2, Nbar;
  double prefac10, T10, f0, tau10, WCO;
  float soblen, length;
  double DV;
  double convFact;

#ifdef XLSOB  
  soblen=cell_sobolev_length(cell);
#else
  soblen = cell_size[cell_level(cell)];
#endif
 
#ifdef CLUMPING
  length = FIXEDLVALUE*constants->pc/units->length;
  n0 *= soblen/length;   
  nH2 *= soblen/length;
#else
  length = soblen; 
#endif 

  Lcm = length*units->length;
  L  = Lcm / 3.08568025e18;
  
  AV = 5.348e-22*Z*n0*Lcm;   
  
  xCO = get_xCO(AV,Lcm,Z,UMW);
  NCO = xCO * n0 * Lcm;   /* in cm^-2 */

#ifdef H2GLOVER
  NH2 = get_xH2(AV,Lcm,Z,UMW)/2 * n0 * Lcm; /* in cm^-2 */
#else
  NH2 = nH2*Lcm;
#endif
   
  /* these constants are taken from formula (2) in Glover and Mac Low 2010
     and with values A10, nu_10 from the Schoeier et al. 2005 database */
  prefac10 = 1.5125e-15; /* cm^2 km s^-1 */
  T10 = 5.532; /* K */
  
  if (z_CMB>=-1e-6) {
    if (z_CMB<0)
      z_CMB=0;
    Tbg = 2.725*(1+z_CMB);
  } else {
    Tbg = 2.725/auni[min_level];
  }

  TB = T10*(1./(exp(T10/T_excitation)-1)-1./(exp(T10/Tbg)-1));	
  
  /* in cm^-2; this includes Helium, i.e. Sigma_bar = m_H * Nbar */
  Nbar =  (n0/constants->XH)*Lcm; 
  
  convFact = 5.847e-12; /* convert sqrt(amu*G*Nbar[cm^-2]*L[pc]) into km/s */ 
  DV = sqrt(0.5*f_grav*Nbar*L)*convFact;
  
  if (DV<min_DV) /* turbulence support */
     DV=min_DV;
    
  f0 = T10 / (2*T_excitation); /* this only holds in LTE, and for Tex>>T0; i.e. for Tex>=10 K */
  tau10 = prefac10 * (1 - exp( - T10/T_excitation ) ) * f0 * NCO / DV;

  WCO = TB * DV * get_escape_fraction(tau10);    
  
  return NH2/WCO; 
  
}

/* like method 0, but rescale velocity with Nbar */
double xfactor_cell_3(int cell, float UMW, int level) { 
  double Z = rtDmw2(cell); /* this is baryon fraction in dust relative to the solar neighbourhood; usually just Z/Zsun */  
  double n0 = units->number_density*cell_gas_density(cell)*constants->XH; /* hydrogen atom density in cm^(-3) */ 
  double nH2 = units->number_density*cell_H2_density(cell); /* molecular hydrogen density in cm^(-3) */
  double AV, Lcm, L;
  double TB, Tbg;
  double xCO, NCO, NH2, Nbar;
  double prefac10, T10, f0, tau10, WCO;
  float soblen, length;
  double DV;
  double convFact;

#ifdef XLSOB  
  soblen=cell_sobolev_length(cell);
#else
  soblen = cell_size[cell_level(cell)];
#endif
 
#ifdef CLUMPING
  length = FIXEDLVALUE*constants->pc/units->length;
  n0 *= soblen/length;   
  nH2 *= soblen/length;
#else
  length = soblen; 
#endif

  Lcm = length*units->length;
  L  = Lcm / 3.08568025e18;
  
  AV = 5.348e-22*Z*n0*Lcm;   
  
  xCO = get_xCO(AV,Lcm,Z,UMW);
  NCO = xCO * n0 * Lcm;   /* in cm^-2 */

#ifdef H2GLOVER
  NH2 = get_xH2(AV,Lcm,Z,UMW)/2 * n0 * Lcm; /* in cm^-2 */
#else
  NH2 = nH2*Lcm;
#endif
   
  /* these constants are taken from formula (2) in Glover and Mac Low 2010
     and with values A10, nu_10 from the Schoeier et al. 2005 database */
  prefac10 = 1.5125e-15; /* cm^2 km s^-1 */
  T10 = 5.532; /* K */
  
  /*
  z=1./auni[min_level]-1.;
  if (z<0) 
    z=0;
  */
  if (z_CMB>=-1e-6) {
    if (z_CMB<0)
      z_CMB=0;
    Tbg = 2.725*(1+z_CMB);
  } else {
    Tbg = 2.725/auni[min_level];
  }

  TB = T10*(1./(exp(T10/T_excitation)-1)-1./(exp(T10/Tbg)-1));	
  
  /* in cm^-2; this includes Helium, i.e. Sigma_bar = m_H * Nbar */
  Nbar =  (n0/constants->XH)*Lcm; 
  
  DV = 0.8 * (Nbar/1e22); /* in km/s; this results in an X-factor ~ constant for large AV & independent of Z */;
  
  if (DV<min_DV) /* turbulence support */
     DV=min_DV;
    
  f0 = T10 / (2*T_excitation); /* this only holds in LTE, and for Tex>>T0; i.e. for Tex>=10 K */
  tau10 = prefac10 * (1 - exp( - T10/T_excitation ) ) * f0 * NCO / DV;

  WCO = TB * DV * get_escape_fraction(tau10);    
 
  return NH2/WCO; 
  
}

/* Krumholz/Wolfire model */
double xfactor_cell_4(int cell, float UMW, int level) {

     double Z = rtDmw2(cell); /* this is baryon fraction in dust relative to the solar neighbourhood; usually just Z/Zsun */  
     double n0 = units->number_density*cell_gas_density(cell)*constants->XH; /* hydrogen atom density in cm^(-3) */ 
 
     double g0n=0.044*(1+3.1*pow(Z,0.365))/4.1;
     double deltaav=0.53-0.045*log(g0n)-0.097*log(Z);
     
     float soblen;
     double Lcm, AV;
     
#ifdef XLSOB  
  soblen=cell_sobolev_length(cell);
#else
  soblen = cell_size[cell_level(cell)];
#endif
     
     Lcm = soblen*units->length;
      
     AV = 5.348e-22*Z*n0*Lcm;  
     
     return 1e20*exp(4.0*deltaav/AV);
}

void config_init_CO()
{
  control_parameter_add2(control_parameter_int,&xfactor_method,"CO:xfactor-method",
  "xfactor_method",
  "method for xfactor. Available methods: \n  0,2,3,4 \nUse <CO:xfactor_method=2>");
  
  strcpy(xfactor_method0.tablePath,"CO_data");
  control_parameter_add2(control_parameter_string,xfactor_method0.tablePath,"CO:table-path",
  "xfactor_method0.tablePath",
  "directory containing the tables for the CO model \nUse <CO:table-path=CO_data>");

  control_parameter_add2(control_parameter_double,&T_excitation,"CO:T-excitation",
  "T_excitation",
  "excitation temperature of the gas in K. \nUse <CO:T-excitation=10>");

  control_parameter_add2(control_parameter_double,&z_CMB,"CO:z-CMB",
  "z_CMB",
  "redshift to compute CMB temperature; if z<0, code uses 1/auni-1. \nUse <CO:z-CMB=0>");

  control_parameter_add2(control_parameter_double,&const_DV,"CO:xfactor-method=0:const-DV",
  "const_DV",
  "constant CO line width in km/s. \nUse <CO:xfactor-method=0:const-DV=3>");

  control_parameter_add2(control_parameter_double,&min_DV,"CO:xfactor-method=2:min-DV",
  "min_DV",
  "minimum CO line width in km/s. \nUse <CO:xfactor-method=2:min_DV=1>");

  control_parameter_add2(control_parameter_double,&f_grav,"CO:xfactor-method=2:f-grav",
  "f_grav",
  "virialization parameter used to convert surface densities into line width. \nUse <CO:xfactor-method-2:f-grav=1>");
}

void config_verify_CO() {
  DIR *d;
  if((d = opendir(xfactor_method0.tablePath)) == NULL) {
                cart_error("Directory %s does not exist.",xfactor_method0.tablePath);
  } else closedir(d);

  VERIFY( CO:xfactor-method, xfactor_method >= 0 && xfactor_method <= 4 );
  VERIFY( CO:T-excitation, T_excitation >= 5. );
  VERIFY( CO:z-CMB, 1 );
  VERIFY( CO:xfactor-method=0:const-DV, const_DV >= 0. );
  VERIFY( CO:xfactor-method=2:min-DV, min_DV >= 0. );
  VERIFY( CO:xfactor-method=2:f-grav, f_grav >= 0. );
}

void setup_xfactor_methods(int level)
{
  char auxString[512];  

  CO_init_table(&xfactor_method0.AV_tab);
  CO_init_table(&xfactor_method0.Z_tab);
  CO_init_table(&xfactor_method0.UMW_tab);
  CO_init_table(&xfactor_method0.L_tab);
 
  if (xfactor_method>=0 && xfactor_method<=4) {

      if (xfactor_method==0)
	xfactor_cell = xfactor_cell_0;
      else if (xfactor_method==1)
       	xfactor_cell = xfactor_cell_1;
      else if (xfactor_method==2)
        xfactor_cell = xfactor_cell_2;	
      else if (xfactor_method==3)	
	xfactor_cell = xfactor_cell_3;
      else if (xfactor_method==4)
        xfactor_cell = xfactor_cell_4;
		
      /* load tables */
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_AV.mat");
      CO_load_table(auxString,&xfactor_method0.AV_tab,READLN);
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_AVL.mat");
      CO_load_table(auxString,&xfactor_method0.AVL_tab,READLN);      
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_L.mat");
      CO_load_table(auxString,&xfactor_method0.L_tab,READLN);
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_Z.mat");
      CO_load_table(auxString,&xfactor_method0.Z_tab,READLN);
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_UMW.mat");
      CO_load_table(auxString,&xfactor_method0.UMW_tab,READLN);
      
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_xCO_full.mat");
      CO_load_table(auxString,&xfactor_method0.xCO_tab,READLN);
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_xH2.mat");
      CO_load_table(auxString,&xfactor_method0.xH2_tab,READLN);            
      
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_tau.mat");
      CO_load_table(auxString,&xfactor_method0.tau_tab,READLN);
      sprintf(auxString,"%s/%s",xfactor_method0.tablePath,"rf_escapefraction.mat");
      CO_load_table(auxString,&xfactor_method0.escapefraction_tab,READLN);
      
   } else {
	xfactor_cell = NULL;
	cart_error("Invalid xfactor_method value: %d",xfactor_method);
   }
}

void unset_xfactor_methods(int level) {
    
   {
    CO_destroy_table(&xfactor_method0.AV_tab);
    CO_destroy_table(&xfactor_method0.AVL_tab);
    CO_destroy_table(&xfactor_method0.L_tab);
    CO_destroy_table(&xfactor_method0.Z_tab);
    CO_destroy_table(&xfactor_method0.UMW_tab);
    
    CO_destroy_table(&xfactor_method0.xCO_tab);
    CO_destroy_table(&xfactor_method0.xH2_tab);

    CO_destroy_table(&xfactor_method0.tau_tab);
    CO_destroy_table(&xfactor_method0.escapefraction_tab);
  }  

}


/* read CO table: note that it is written in column major format -> need to permute (invert) dimensions */
/* read log of the variable if readLn=1 */
void CO_load_table(const char *filename, COTable * tab, int readLn) {
  
 FILE *file;
 unsigned long int numel=1,ii;
 int idim;

 
 file = fopen( filename, "r" );
 fprintf(stderr,"Reading %s ---\n",filename);
 cart_assert(file != NULL);
 
 fread(&tab->num_dims, 4, 1, file);
 
 if (tab->dims)
   free(tab->dims);
 tab->dims = (unsigned*)malloc(sizeof(unsigned)*(tab->num_dims));  

 for (idim=0;idim<tab->num_dims;idim++) {
   fread( &tab->dims[tab->num_dims-1-idim] , 4, 1, file);
 }  
 
 for (idim=0;idim<tab->num_dims;idim++) {
   numel *= tab->dims[idim];
 }
 
 if (tab->data)
   free(tab->data);
 tab->data = (float*)malloc(sizeof(float)*numel);  

 fread(tab->data, 4, numel, file);
 
 if (readLn){
   for (ii=0;ii<numel;ii++) {
     tab->data[ii]=log(tab->data[ii]);
   }
 }
 
 fclose(file);
}

void CO_init_table(COTable * tab) {
  tab->dims=NULL;
  tab->data=NULL;
  tab->num_dims=0;
}

void CO_destroy_table(COTable * tab) {

 if (tab->dims) {
   free(tab->dims);
   tab->dims=NULL;
 }
 if (tab->data) {
   free(tab->data);
   tab->data=NULL;
 }
 tab->num_dims=0;
 
}

#endif /* HYDRO && CO */
