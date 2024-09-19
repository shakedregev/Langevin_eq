#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define KBT 1
#define m 1.0
#define malpha 700.0
#define nalpha 1.0
#define leftdt (float)2*m/malpha
#define rightdt (float)2*m/nalpha
#define V0 sqrt(KBT/m)
#define len1 1.0
#define len2 20.0
#define pos 21
#define runs 10000
#define TIME 10000
#define maxtime 1200.0
#define mu 10
#define start 0.4999
#define min 100
#define mydt 0.0001
#define PI 2*asin(1)
#define MBIG 10070030850
#define MSEED 127384609
#define MZ 0
#define FAC (1.0/MBIG)
#define dim 1
#define bins 18
#define deltar 0.1
double time=0;
double alpha(double x,double y)
{ 
    if (x<len1/2){
      if (y<len1/2)
	return malpha;
      return ((-(x-len1/2)*malpha+(y-len1/2)*nalpha)/(y-x));
  }
    if (y>len1/2)
      return nalpha;
    return  (((x-len1/2)*nalpha-(y-len1/2)*malpha)/(x-y));
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
double expdev(int *idum){
double dum;
do
dum=ran3(idum);
while (dum == 0.0);
return -log(dum);
}
void advance(double x[2],double v[2],double dt, int* idum){
  double b=0,a=0,sigma=0,tempx=0,alp=0;
  x[0]=x[1];
  v[0]=v[1];
  tempx=x[0]+v[0]*dt;
  alp=alpha(x[0],tempx);
  b=coefb(dt,alp);
  a=b*coefa(dt,alp);
  sigma=myrand(idum);
  x[1]=x[0]+b*dt*v[0]+b*(dt/2/m)*sqrt(2*alp*KBT*dt)*sigma;
  v[1]=a*v[0]+b/m*sqrt(2*alp*KBT*dt)*sigma;
  if (x[1]<-len1/2){
    x[1]=-len1-x[1];
    v[1]*=-1;
  }
  return;
}
int main()
{
  double x[2]={0,0}, v[2]={0,0};
  int i=0, j=0,k=0,n=0,seedn=-1,counter=0,temp,hist[2][bins]={0},sample[bins]={2,3,4,6,8,12,20,30,40,60,80,120,200,300,400,600,800,1200};
  int *idum=NULL;
  FILE *out,*out2;
  idum=&seedn;
   out=fopen("coat3","w");
   out2=fopen("artery3","w");
  for(j=0;j<pos;j++){
    for (k=0;k<runs;k++){
        time=0;
        counter=0;
      v[1]=myrand(idum);
      x[1]=-len1/2+0.05*j;
      if(x[1]==-len1/2)
          x[1]+=0.0001;
    if(x[1]==len1/2)
        x[1]-=0.0001;
	while(1){
        advance(x,v,mydt,idum);
        if (x[1]>len1/2 && x[0]<=len1/2){
        time+=round(expdev(idum)*mu/mydt)*mydt;
        while ((int)time>=sample[counter])
        {
            hist[0][counter]++;
            counter++;
        }
        }
        time+=mydt;
        if(x[1]>=len2+len1/2)
        {
            break;
        }
        if((int)time==sample[counter])
        {
         if (x[1]<=len1/2)
             hist[0][counter]++;
         else
             hist[1][counter]++;
        counter++;
        }
         if(time>=maxtime)
        {
            break;
        }
	}
    }
 }
 for (j=0;j<bins;j++)
 {
     fprintf(out,"%d\n",hist[0][j]); 
     fprintf(out2,"%d\n",hist[1][j]); 
 }
  fclose(out);
  fclose(out2);
}
