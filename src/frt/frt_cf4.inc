C
C  Configuration is in this include
C
	include 'frt_cf_LHOSI13.inc'
C
	logical MetalPartOnly
	parameter(MetalPartOnly = .false.)
C
C  Number of components per T-D bin
C
	parameter(NC = 6)
C
C  Data
C
	common/RT_CF_A/ indx(NP1,NP2,NP3,NP4)
	common/RT_CF_B/ data(NN,NC,NT,ND,NX)
	common/RT_CF_C/ qmin(4), qmax(4), qstp(4)
	common/RT_CF_D/ alt(NT), altmin, altstp
	common/RT_CF_E/ np(4)
	common/RT_CF_F/ xmin, xmax, xstp
