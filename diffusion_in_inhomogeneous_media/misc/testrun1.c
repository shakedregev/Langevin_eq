#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define KBT 1
#define m 1.0
#define malpha 2
#define nalpha 16
#define leftdt (float)2*m/malpha
#define rightdt (float)2*m/nalpha
#define V0 sqrt(KBT/m)
#define NUMRUNS 1
#define TIME 100000
#define min 256
#define PI 2*asin(1)
#define MBIG 10070030850
#define MSEED 127385694
#define MZ 0
#define FAC (1.0/MBIG)
#define dim 1
#define bins 640
#define deltar 0.1
double alpha(double x,double y,int flag)
{ 
  if (flag==0){
    if (x<0){
      if (y<0)
	return malpha;
      return ((-x*malpha+y*nalpha)/(y-x));
  }
    if (y>0)
      return nalpha;
    return  ((x*nalpha-y*malpha)/(x-y));
  }
  if (x<y){
    return (((32.0-x)*nalpha+(y-32.0)*malpha)/(y-x));
  }
   return (((x+32.0)*malpha+(-32.0-y)*nalpha)/(x-y));
}
double coefb(double dt,double alp)
{
   return 1/(1+alp*dt/2);
}
double coefa(double dt,double alp)
{
  return 1-alp*dt/2;
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
void advance(double x[2],double v[2],double dt, int* idum){
  double b=0,a=0,sigma=0,tempx=0,alp=0;
  x[0]=x[1];
  v[0]=v[1];
  tempx=x[0]+v[0]*dt;
  if(fabs(tempx)>32)
    alp=alpha(x[0],tempx,1);
  else
  alp=alpha(x[0],tempx,0);
  b=coefb(dt,alp);
  a=b*coefa(dt,alp);
  sigma=myrand(idum);
  x[1]=x[0]+b*dt*v[0]+b*(dt/2)*sqrt(2*alp*KBT*dt)*sigma;
  v[1]=a*v[0]+b*sqrt(2*alp*KBT*dt)*sigma;
  if(x[1]>32)
    x[1]-=64;
  if(x[1]<-32)
    x[1]+=64;
  return;
}
void timestep(double x[2],double v[2],double dt,int* idum){
  if(dt<=leftdt/min){
    advance(x,v,dt,idum);
    return;
  }

  if ((fabs(x[1])>8*V0*dt) && (fabs(x[1])<32-8*V0*dt)){
      advance(x,v,dt,idum);
      return;
    }

  timestep(x,v,dt/2,idum);
  timestep(x,v,dt/2,idum);
  return;
}
int main()
{
  double x[2]={0,0}, v[2]={0,0};
  int i=0, j=0,k=0,n=0,seedn=-1,counter=0,temp,hist[bins]={0};
  int *idum=NULL;
  FILE *ifp1;
  idum=&seedn;
  for(j=0;j<NUMRUNS;j++){
    for (k=0;k<NUMRUNS;k++){
      v[1]=myrand(idum);
      x[1]=0;
      for (i=0;i<TIME;i++)
	for(n=0;n<TIME;n++){
	 if (x[1]<=0)	
	   timestep(x,v,leftdt,idum);
	 else{
	   for (counter=0;counter<8;counter++)
	    timestep(x,v,rightdt,idum);
	 }
	  temp=(int)(bins/2+round(x[1]/deltar-deltar/0.2));
	  hist[temp]++;
	}
    }
 }
   ifp1=fopen("xscon","w");
   for (i=0;i<bins;i++)
     fprintf(ifp1,"%d\n",hist[i]);
   fclose(ifp1);
}
