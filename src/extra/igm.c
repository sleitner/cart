#include "defs.h"

#include <stdio.h>
#include <mpi.h>

#include "auxiliary.h"
#include "parallel.h"
#include "timestep.h"


#ifdef RADIATIVE_TRANSFER
#include "rt_c2f.h"

void f2c_wrapper(frtquerybackground)(f2c_intg *n, f2c_real *wlen, f2c_real *jnu, f2c_real *aexp);

#endif



#ifdef RADIATIVE_TRANSFER
void extDumpRadiationBackground()
{
  int i;
  f2c_intg n;
  f2c_real *wlen, *jnu, aa;
  FILE *f;

  aa = aexp[0];

  if(local_proc_id == MASTER_NODE)
    {
      n = 0;
      f2c_wrapper(frtquerybackground)(&n,wlen,jnu,&aa);
      
      wlen = cart_alloc(n*sizeof(f2c_real));
      jnu = cart_alloc(n*sizeof(f2c_real));
      f2c_wrapper(frtquerybackground)(&n,wlen,jnu,&aa);

      f = fopen("jnu.res","w");
      cart_assert(f != NULL);

      fprintf(f,"# wlen[A]  jnu[cgs] at a=%lf\n",aexp[0]);
      for(i=0; i<n; i++)
	{
	  fprintf(f,"%9.3e %9.3e\n",wlen[i],jnu[i]);
	}

      fclose(f);

      cart_free(wlen);
      cart_free(jnu);
    }
}
#endif

