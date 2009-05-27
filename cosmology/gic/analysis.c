#include "defs.h"

#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "tree.h"

#include "extra/ism.h"
#include "extra/igm.h"

extern int num_options;
extern char **options;


void analysis()
{
  int i, level;
  float s, t;
  const char *tmp;

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
      tmp = check_option0(options[i],"chem");
      if(tmp != NULL)
	{
#ifdef RADIATIVE_TRANSFER
	  level = max_level_now_global(MPI_COMM_WORLD);
	  extDumpChemicalState(level,level);
#else
	  cart_error("-chem requires RADIATIVE_TRANSFER on");
#endif
	  continue;
	}

      tmp = check_option1(options[i],"kenn",NULL);
      if(tmp != NULL)
	{
#if defined(PARTICLES) && defined(STARFORM)
	  if(sscanf(tmp,"%e,%e",&s,&t) != 2)
	    {
	      cart_error("-kenn option format is: -kenn=<scale>,<time>");
	    }
	  extDumpKennicuttLaw(s,t);
#else
	  cart_error("-kenn requires PARTICLES && STARFORM on");
#endif
	  continue;
	}

      tmp = check_option0(options[i],"jnu");
      if(tmp != NULL)
	{
#ifdef RADIATIVE_TRANSFER
	  extDumpRadiationBackground();
#else
	  cart_error("-jnu requires RADIATIVE_TRANSFER on");
#endif
	  continue;
	}

      cart_error("Invalid option %s",options[i]);
    }
}
