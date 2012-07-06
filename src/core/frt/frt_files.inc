C
C  Unfortunately, I have to hardcode these file locations here,
C  as g77 has a bug that prevents proper string manipulation
C
	character*(*) fileCF, fileLW, fileHM
	parameter(fileCF = 'drt/cf_table.I2.dat')
	parameter(fileLW = 'drt/lwbands.dat')
	parameter(fileHM = 'drt/bkg_HM.dat')

