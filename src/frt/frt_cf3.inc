C
C  Configuration is in this include
C
	include 'frt_cf_LHO13.inc'
C
	logical MetalPartOnly
	parameter(MetalPartOnly = .true.)
C
C  Number of components per T-D bin
C
	parameter(NC = 6)
C
C  Data
C
	common/RT_CF_A/ indx(NP1,NP2,NP3)
	common/RT_CF_B/ data(NN,NC,NT,ND,NX)
	common/RT_CF_C/ qmin(3), qmax(3), qstp(3)
	common/RT_CF_D/ alt(NT), altmin, altstp
	common/RT_CF_E/ np(3)
	common/RT_CF_F/ xmin, xmax, xstp
