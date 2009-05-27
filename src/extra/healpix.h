/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2008 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */

/*
//  05202009 Modified by NG: names changed to avoid name clash
//                           eliminated extra calls
*/

#ifndef __HEALPIX_H__
#define __HEALPIX_H__

/* -------------------- */
/* Constant Definitions */
/* -------------------- */

#ifndef HEALPIX_NULLVAL
#define HEALPIX_NULLVAL (-1.6375e30)
#endif /* HEALPIX_NULLVAL */

/* --------------------- */
/* Function Declarations */
/* --------------------- */

/* pixel operations */
/* ---------------- */
void hp_ang2pix_nest(const long nside, double theta, double phi, long *ipix);
void hp_pix2ang_nest(long nside, long ipix, double *theta, double *phi);

#endif /* __HEALPIX_H__ */

