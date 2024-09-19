#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define KBT 1
#define NUMRUNS 500
#define TIME 1000000
#define PI 2*asin(1)
#define dt 0.001
#define MBIG 10070030850
#define MSEED 123785692
#define MZ 0
#define FAC (1.0/MBIG)
#define dim 1
#define bins 500000
#define deltar 0.1
double alpha(double x,double y)
{
  if(x/y==1)
    return pow(fabs(x),-0.25);
  return fabs((4.0/3.0*(pow(fabs(x),0.75)-pow(fabs(y),0.75)))/(x-y));
}
double coefb(double malpha)
{
   return 1/(1+malpha*dt/2);
}
double coefa(double malpha)
{
  return 1-malpha*dt/2;
}
double ran3(int* idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i=0,ii=0,k=0;
  if (*idum<0 || iff==0){
    iff=1;
    mj=MSEED;
    mj%=MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i)%55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk<MZ) mk+=MBIG;
	mj=ma[ii];
      }
    for (k=1;k<=4;k++)
      for(i=1;i<=55;i++){
      ma[i]-=ma[1+(i+30)%55];
      if (ma[i]<MZ) ma[i]+=MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if(++inext==56) inext=1;
  if (++inextp==56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if(mj<MZ) mj+=MBIG;
  ma[inext]=mj;
  return mj*FAC;  
}
double  myrand(int *idum)
{
double  x1, x2, w, y1, y2;
         do {
                 x1 = 2.0 * ran3(idum) - 1.0;
                 x2 = 2.0 * ran3(idum) - 1.0;
                 w = x1 * x1 + x2 * x2;
	         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;
	 return y1;

}
int main()
{
  double xnew=0, xlast=0, vnew=0,vlast=0;
  double tempx=0,malpha=0,b=0,a=0,r=0,sigma=0;
  int i=0, j=0,k=0,m=0,seedn=-1;
  int *idum=NULL;
  FILE *ifp,*ifp2;
  idum=&seedn;
   ifp=fopen("rs0.25","w");
   ifp2=fopen("v0.25","w");
  for(m=0;m<NUMRUNS;m++){
    for (k=0;k<NUMRUNS;k++){
      vnew=myrand(idum);
          xnew=0;
      for (i=0;i<TIME;i++){
	xlast=xnew;
	  vlast=vnew;
     tempx=xlast+vlast*dt;
     malpha=alpha(xlast,tempx);
      b=coefb(malpha);
     a=b*coefa(malpha);
     sigma=myrand(idum);
     xnew=xlast+b*dt*vlast+b*(dt/2)*sqrt(2*malpha*KBT*dt)*sigma;
     vnew=a*vlast+b*sqrt(2*malpha*KBT*dt)*sigma;
 if ((i==100) || (i==200) || (i==500) || (i==1000) || (i==2000) || (i==5000) || (i==10000) || (i==20000) || (i==50000) || (i==100000) || (i==200000) || (i==500000)){
       fprintf(ifp,"%f\n",xnew);
 fprintf(ifp2,"%f\n",vnew);
         }
    }
    fprintf(ifp,"%f\n",xnew);
 fprintf(ifp2,"%f\n",vnew);
 }
 }
   fclose(ifp);
}
