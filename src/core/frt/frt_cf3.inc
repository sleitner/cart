C
C  Configuration is in this include
C
	include 'frt_cf3.I2.inc'
C
C  Internally used unit number
C
	parameter(IOCF = 97)
C
C  Number of components per T-D bin
C
	parameter(NC = 6)
	parameter(NICH = 12)
	parameter(NRCH = 13)
C
C  Data
C
	common/RT_CF_A/ indx(NP1,NP2,NP3)
	common/RT_CF_C/ data(NC,NT,NX,ND)
	common/RT_CF_D/ np(3), mode
	common/RT_CF_E/ altval(NT), altmin, altstp,
     .			xval(NX), xmin, xmax, xstp,
     .			qmin(3), qmax(3), qstp(3)
