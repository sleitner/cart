#include "config.h"

#include <stdio.h>

#include "auxiliary.h"

/*
// In the analysis mode the code is run as
//
//   analysis file.cfg [--fast] <snapshot1> [snapshot2] ... [snapshorN] [analysis options] 
//
// where <snapshot> labels are either scale factors or step numbers depending
// on whether COSMOLOGY is declared or not.
// For each snapshot main_analysis is called with its (argc,argv) arguments 
// including only [analysis options].
// Global option --fast will result in fast initialization, when densities
// are not assigned and RT is not initialized.
*/
int main_analysis(int argc, char **argv)
{

  cart_debug("In analysis mode");

  return 0;
}
