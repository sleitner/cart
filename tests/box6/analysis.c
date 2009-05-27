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
  int i, id;
  float s, t;
  char *str;
  const char *tmp;

  if(num_options < 1)
    {
      if(local_proc_id == MASTER_NODE)
	{
	  printf("An option is required. Available options:\n");
	  printf("  -chem:  dump checmical state of the gas at max_level\n");
	  printf("  -kenn=<scale>,<time>:  dump Kennicutt law at max_level\n");
	  printf("  -jnu:  dump average cosmic radiation background\n");
	  printf("  -pz=<hlist-file>[#<halo-id>]:  compute proximity zones (optionally for a specific halo)\n");
	}
    }

  for(i=0; i<num_options; i++)
    {
      tmp = check_option0(options[i],"chem");
      if(tmp != NULL)
	{
	  extDumpChemicalState(max_level,max_level);
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

      tmp = check_option1(options[i],"pz",NULL);
      if(tmp != NULL)
	{
#ifdef RADIATIVE_TRANSFER
	  str = strchr(tmp,'#');
	  if(str == NULL)
	    {
	      id = 0;
	    }
	  else
	    {
	      *str = 0;
	      str++;
	      id = 0;
	      if(sscanf(str,"%d",&id)==0 || id<1)
		{
		  cart_error("-pz=<hlist-file>#<halo-id> requires positive integer <halo-id> as an argument.");
		}
	    }
	  extFindProximityZones(tmp,1,id);
#else
	  cart_error("-pz requires RADIATIVE_TRANSFER on");
#endif
	  continue;
	}

      cart_error("Invalid option %s",options[i]);
    }
}
