C
C  Functions needed for secondary electrons block
C
	se1fun(x,a,b,c) = a*(1.00001-x**b)**c
	se2fun(x,a,b,c) = a*(x*1.0e4)**b*(1.0-x**c)**2
	se1der(x,b,c) = -c*b*x**(b-1.0)/(1.00001-x**b)
	se2der(x,b,c) = (b-2*c*x**c/(1.00001-x**c))/x

	
	
