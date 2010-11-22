
#include <math.h>
#include <stdio.h>


void accel(double x[2], double a[2])
{
  int j;
  double r, p;

  r = sqrt(1.0e-4+x[0]*x[0]+x[1]*x[1]);
  p = 4*pow(M_PI,3.0)/(3*r);

  for(j=0; j<2; j++)
    {
      a[j] = -x[j]*p/(r*r);
    }
}


int main()
{
  const int level = 3;
  const float vel0 = 1.35;
  const float fac0 = 2.0;
  double x[2], v[2], a[2],  t, dt;
  FILE *f;
  int j;

  x[0] = 1.0/fac0;
  x[1] = 0.0;

  v[0] = 0.0;
  v[1] = 6.2830*vel0*sqrt(fac0);

  f = fopen("orbit.res","w");

  t = 0.0;
  dt = 0.1/(1 << level);

  accel(x,a);

  while(t < 100.0)
    {
      printf("%9.3le %6.4lf %6.4lf\n",t,x[0],x[1]);
      fprintf(f,"%9.3le %6.4lf %6.4lf %6.4lf %le %le %le %d\n",t,x[0],x[1],0.0,v[0],v[1],0.0,level);

      for(j=0; j<2; j++)
	{
	  v[j] += 0.5*dt*a[j];
	  x[j] += dt*v[j];
	}

      accel(x,a);

      for(j=0; j<2; j++)
	{
	  v[j] += 0.5*dt*a[j];
	}

      t += dt;
      
    }

  fclose(f);

}

