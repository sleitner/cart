C
C  Sampling frequencies
C
	parameter(NSETS = 3)
	common /RT_TRF_A/ LRSETS(NSETS)
C
C  Radiation Field parameters (for internal packing)
C
        parameter(irfNorm =  1)
        parameter(irfH1lf =  2)
        parameter(irfG1lf =  3)
        parameter(irfG2lf =  4)
        parameter(irfH1gf =  5)
        parameter(irfG1gf =  6)
        parameter(irfG2gf =  7)
        parameter(irfH1hw =  8)
        parameter(irfG1hw =  9)
        parameter(irfG2hw = 10)
        parameter(irfDim  = 10)
