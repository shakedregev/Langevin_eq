#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define num 8
#define KBT 1
#define m 1.0
#define d1 1.0/700
#define d2 1.0
#define len1 1.0
#define len2 20.0
#define tmax 180001
#define mydt 0.01
#define PI 2*asin(1)
double flux(double x,double t,double D)
{ 
    if (t==0)
        return 0;
    return x/2/t/2/sqrt(PI*t*D)*exp(-x*x/4/t/D);
}

int main()
{
  double j1[tmax]={0}, j2[tmax]={0},jt[tmax]={0},integ1[tmax]={0},sum=0,time=0;
  int i=0, k=0,n=0,counter=0;
FILE *out;
   out=fopen("data4","w");
for (n=1;n<tmax;n++)
{
    sum=0;
    for(i=1;i<=n;i++)
    {
    for (k=-num;k<=num;k++)
    {
        j1[i]+=flux(len1*k+len1/2,i*mydt,d1)*pow(-1,((k+num)%4)/2);
        j2[i]+=flux(2*len2*k+len2,(n-i)*mydt,d2)*pow(-1,k);
    }
    jt[i]=j1[i]*j2[i];
    j1[i]=j2[i]=0;
    sum+=(jt[i]+jt[i-1])*mydt/2;
    }
    integ1[n]=sum;
}
sum=0;
for(i=1;i<tmax;i++)
{
   sum+=(integ1[i]+integ1[i-1])*mydt/2; 
   if(i%5000==0)
       fprintf(out,"%f\n",sum); 
}
  fclose(out);
}
