#include "defs.h"

#include "tree.h"

#include "extra/ism.h"


void analysis()
{
#ifdef RADIATIVE_TRANSFER
  extDumpChemicalState(max_level,max_level);
#endif

#if defined(PARTICLES) && defined(STARFORM)
  extDumpKennicuttLaw(0.5,30.0);
#endif
}
