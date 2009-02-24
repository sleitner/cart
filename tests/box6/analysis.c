#include "defs.h"

#include <string.h>

#include "tree.h"

#include "extra/ism.h"
#include "extra/igm.h"

extern int num_options;
extern char **options;


void analysis()
{
  int i;
  float s, t;

  if(num_options < 1)
    {
      if(local_proc_id == MASTER_NODE)
	{
	  printf("An option is required. Available options:\n");
	  printf("  -chem:  dump checmical state of the gas at max_level\n");
	  printf("  -kenn:<scale>,<time>:  dump Kennicutt law at max_level\n");
	  printf("  -jnu:  dump average cosmic radiation background\n");
	}
    }

  for(i=0; i<num_options; i++)
    {

      if(strcmp(options[i],"-chem") == 0)
	{
#ifdef RADIATIVE_TRANSFER
	  extDumpChemicalState(max_level,max_level);
#else
	  cart_error("-chem requires RADIATIVE_TRANSFER on");
#endif
	}

      if(strncmp(options[i],"-kenn",5) == 0)
	{
#if defined(PARTICLES) && defined(STARFORM)
	  if(options[i][5]!=':' || sscanf(options[i]+6,"%e,%e",&s,&t) != 2)
	    {
	      cart_error("-kenn option format is: -kenn:<scale>,<time>");
	    }
	  extDumpKennicuttLaw(s,t);
#else
	  cart_error("-kenn requires PARTICLES && STARFORM on");
#endif
	}

      if(strcmp(options[i],"-jnu") == 0)
	{
#ifdef RADIATIVE_TRANSFER
	  extDumpRadiationBackground();
#else
	  cart_error("-jnu requires RADIATIVE_TRANSFER on");
#endif
	}
    }
}
