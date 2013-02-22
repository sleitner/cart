#ifndef __CELL_FEEDBACK_STEP_H__
#define __CELL_FEEDBACK_STEP_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif

void cell_feedback(int level);
#ifdef RT_OTVET_SAVE_FLUX
void RaP_longrange_setup(int level);   
#endif /* RT_OTVET_SAVE_FLUX */


#endif
