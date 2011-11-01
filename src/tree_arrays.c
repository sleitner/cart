#include "config.h"

#include "tree.h"

#ifdef GCC_COMPILER
      int external_direction[num_children][nDim]
#else
const int external_direction[num_children][nDim]
#endif 
= {
	#if nDim == 1
		{ 0 }, { 1 }
	#elif nDim == 2
		{ 0, 1, 0, 1 }, {  1, 0, 0, 1 },
		{ 0, 1, 1, 0 }, {  1, 0, 1, 0 }
	#elif nDim == 3
		{ 0, 2, 4 }, { 1, 2, 4 }, { 0, 3, 4 }, { 1, 3, 4 },
		{ 0, 2, 5 }, { 1, 2, 5 }, { 0, 3, 5 }, { 1, 3, 5 }
	#else
		#error "Unsupported number of dimensions for external_neighbors"
	#endif
};
#ifdef GCC_COMPILER
      int reverse_direction[num_neighbors] 
#else
const int reverse_direction[num_neighbors] 
#endif
= {
	#if nDim == 1
		1, 0
	#elif nDim == 2
		1, 0, 3, 2
	#elif nDim == 3
		1, 0, 3, 2, 5, 4
	#else
		#error "No valid reverse_direction for that number of dimensions!"
	#endif
};

#ifdef GCC_COMPILER
      int uniform_stencil[num_stencil][nDim] 
#else
const int uniform_stencil[num_stencil][nDim] 
#endif
= {
	#if nDim == 3 
		{ -1, -1, -1 }, {  0, -1, -1 }, {  1, -1, -1 }, 
		{ -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 }, 
		{ -1,  1, -1 }, {  0,  1, -1 }, {  1,  1, -1 }, 	
		{ -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 }, 
		{ -1,  0,  0 },                 {  1,  0,  0 }, 
		{ -1,  1,  0 }, {  0,  1,  0 }, {  1,  1,  0 },  
		{ -1, -1,  1 }, {  0, -1,  1 }, {  1, -1,  1 }, 
		{ -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 }, 
		{ -1,  1,  1 }, {  0,  1,  1 }, {  1,  1,  1 },  
		{  2,  0,  0 }, { -2,  0,  0 }, {  0,  2,  0 },
		{  0, -2,  0 }, {  0,  0,  2 }, {  0,  0, -2 }
	#else
		#error "No valid stencil for that number of dimesions!"
	#endif
};

#ifdef GCC_COMPILER
      int secondary_neighbors[num_secondary_neighbors][2] 
#else
const int secondary_neighbors[num_secondary_neighbors][2] 
#endif
= {
    #if nDim == 3
        { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 },
        { 4, 0 }, { 4, 1 }, { 4, 2 }, { 4, 3 }, 
        { 5, 0 }, { 5, 1 }, { 5, 2 }, { 5, 3 }, 
    #else
        #error "No valid secondary_neighbors for that number of dimensions"
    #endif
};

#ifdef GCC_COMPILER
      int tertiary_neighbors[num_tertiary_neighbors][2] 
#else
const int tertiary_neighbors[num_tertiary_neighbors][2] 
#endif
= {
    #if nDim == 3
        { 0, 4 }, { 0, 5 },
        { 1, 4 }, { 1, 5 },
        { 2, 4 }, { 2, 5 },
        { 3, 4 }, { 3, 5 }
    #else
        #error "No valid tertiary_neighbors for that number of dimensions"
    #endif
};


#ifdef GCC_COMPILER
      int secondary_external_neighbors[num_children][nDim] 
#else
const int tertiary_neighbors[num_tertiary_neighbors][2] 
#endif
= {
	#if nDim == 3
		{ 2, 4, 0 }, { 2, 4, 1 },
		{ 3, 4, 0 }, { 3, 4, 1 },
		{ 2, 5, 0 }, { 2, 5, 1 },
		{ 3, 5, 0 }, { 3, 5, 1 }
	#else
		#error "No valid secondary_external_neighbors for that number of dimensions"
	#endif
};

#ifdef GCC_COMPILER
      int pyramid_vertices[num_children][nDim] 
#else
const int pyramid_vertices[num_children][nDim] 
#endif
= {
	#if nDim == 3
		{ 0, 2, 4 }, { 1, 2, 4 }, { 0, 3, 4 }, { 1, 3, 4 },
		{ 0, 2, 5 }, { 1, 2, 5 }, { 0, 3, 5 }, { 1, 3, 5 }
	#else
		#error Unsupported number of dimensions for pyramid_vertices 
	#endif
};

#ifdef GCC_COMPILER
      int local[num_children][num_neighbors] 
#else
const int local[num_children][num_neighbors] 
#endif
= {
        #if nDim == 1
                { 1, 1 }, { 0, 0 }
        #elif nDim == 2
                { 1, 1, 2, 2 }, {  0, 0, 3, 3 },
                { 3, 3, 0, 0 }, {  2, 2, 1, 1 }
        #elif nDim == 3
                { 1, 1, 2, 2, 4, 4 }, { 0, 0, 3, 3, 5, 5 },
                { 3, 3, 0, 0, 6, 6 }, { 2, 2, 1, 1, 7, 7 },
                { 5, 5, 6, 6, 0, 0 }, { 4, 4, 7, 7, 1, 1 },
                { 7, 7, 4, 4, 2, 2 }, { 6, 6, 5, 5, 3, 3 }
        #else
                #error "Unsupported number of dimensions for local"
        #endif
};
                                                                                                                                                             
                                                                                                                                                             
#ifdef GCC_COMPILER
      int in_local_oct[num_children][num_neighbors] 
#else
const int in_local_oct[num_children][num_neighbors] 
#endif
= {
        #if nDim == 1
                { 0, 1 }, { 1, 0 }
        #elif nDim == 2
                { 0, 1, 0, 1 }, {  1, 0, 0, 1 },
                { 0, 1, 1, 0 }, {  1, 0, 1, 0 }
        #elif nDim == 3
                { 0, 1, 0, 1, 0, 1 }, { 1, 0, 0, 1, 0, 1 },
                { 0, 1, 1, 0, 0, 1 }, { 1, 0, 1, 0, 0, 1 },
                { 0, 1, 0, 1, 1, 0 }, { 1, 0, 0, 1, 1, 0 },
                { 0, 1, 1, 0, 1, 0 }, { 1, 0, 1, 0, 1, 0 }
        #else
                #error "Unsupported number of dimensions for in_local_oct"
        #endif
};
                                                                                                                                                             
                                                                                                                                                             
#ifdef GCC_COMPILER
      int ishift[num_neighbors][nDim] 
#else
const int ishift[num_neighbors][nDim] 
#endif
= {
        #if nDim == 1
                { -1 }, { 1 }
        #elif nDim == 2
                { -1,  0 }, {  1,  0 },
                {  0, -1 }, {  0, -1 }
        #elif nDim == 3
                { -1,  0,  0 }, {  1,  0,  0 }, {  0, -1,  0 },
                {  0,  1,  0 }, {  0,  0, -1 }, {  0,  0,  1 }
        #else
                #error "Unsupported number of dimensions for ishift"
        #endif
};

/* array which describes how child cells are offset from
* the center of their parent oct */
#ifdef GCC_COMPILER
      double cell_delta[num_children][nDim] 
#else
const double cell_delta[num_children][nDim] 
#endif
= {
	#if nDim == 1
		{ -0.5 }, { 0.5 }
	#elif nDim == 2
		{ -0.5, -0.5 }, { 0.5, -0.5 }, { -0.5, 0.5 }, { 0.5, 0.5 }
	#elif nDim == 3
		{ -0.5, -0.5, -0.5 }, {  0.5, -0.5, -0.5 }, { -0.5,  0.5, -0.5 },
		{  0.5,  0.5, -0.5 }, { -0.5, -0.5,  0.5 }, {  0.5, -0.5,  0.5 },
		{ -0.5,  0.5,  0.5 }, {  0.5,  0.5,  0.5 }
	#else
		#error "No valid cell_delta for that number of dimensions!"
	#endif
};

/* this array is used to find the lower left corner cell */
#ifdef GCC_COMPILER
      int neighbor_moves[num_children][nDim] 
#else
const int neighbor_moves[num_children][nDim] 
#endif
= {
	{ 0, 2, 4 }, { 2, 4, -1 }, { 0, 4, -1 }, { 4, -1, -1 },
	{ 0, 2, -1 }, { 2, -1, -1 }, { 0, -1, -1 }, { -1, -1, -1 }
};

