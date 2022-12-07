/* construction regresion line Y = b0 + b1*X               */
/* calculation parameters b0,b1 by samples Xi,Yi i=1,...N  */

 int ParRegres(int N, double *X,double *Y,double *pb0,double *pb1)
 //int N;
 //float X[],Y[];
 //float *pb0,*pb1;
{
int ierr;  /* my universal error code */
double b0,b1;                /* regression parameters */
double sy,sx,sxy,sx2,znamen; /* Sum of Y-s, Sum of X-s, Sum of X*Y, Sum squares of X, 'znamenatel' */
double xi,yi;
int i;

if(N<2) ierr=177;
if(ierr!=1) return(ierr);

/* PRACTICE FORMULAES: -----------------------------  */
/*                                                    */ 
/*       (Sum Yi)*(Sum Xi**2)-(Sum Xi)*(Sum Xi*Yi)    */
/* b0 =  -----------------------------------------    */
/*             N*(Sum Xi**2) - (Sum Xi)**2            */
/*                                                    */ 
/*            N* (Sum Xi*Yi) - (Sum Xi)*(Sum Yi)      */
/* b1 =  -----------------------------------------    */
/*             N*(Sum Xi**2) - (Sum Xi)**2            */
/*____________________________________________________*/

xi=0;
yi=0;
sx=0;
sy=0;
sxy=0;
sx2=0;
for(i=0;i<N;i++)
  {
  xi=X[i];             /*\ local for calculation speed */
  yi=Y[i];             /*/                             */
  sx+=xi;              /*\                             */
  sy+=yi;              /*-\ calculation sums           */
  sxy+=xi*yi;          /*-/ need for formulaes         */
  sx2+=xi*xi;          /*/                             */
  }

 znamen=N*sx2-sx*sx;
  if(znamen==0) { ierr=178; return(ierr); }  /* wrong sample X : all Xi=const */
 b0=(sy*sx2-sx*sxy)/znamen;
 b1=(N*sxy-sx*sy)/znamen;
/*-- regression parameters (Y=b0+b1*X) are calculated --*/

*pb0=b0;               /*\ pass results */
*pb1=b1;               /*/              */ 
return(0);
}

double X[10001];//10K
double X_i_n[10001];//10K
double Y[10001];//10K
double RS[10001];//10K
int maskRY[4]={0,1,0,1};//RY -check ATGC!


double mean0(int N, double *X)
{
int i;double tmp=0;
for (i=0;i<N;i++)
	tmp+=X[i];
tmp/=N;
return(tmp);
}

double expHurst(int N, double *X)
{
int i,j,n_;
double tmp,b0=-2,b1=-3;
double R,S,maxX_i_n,minX_i_n,xmean;

for(n_=2;n_<N;n_++) //srart from at least 2 points X
{
xmean=mean0(n_,X);//mean

tmp=0;
for(j=0;j<n_;j++) tmp+=(X[j]-xmean);
X_i_n[n_]=tmp;//deviation

maxX_i_n=X_i_n[0];minX_i_n=X_i_n[0];//just to start
for(j=0;j<n_;j++) 
	{
	if(maxX_i_n<X_i_n[j])maxX_i_n=X_i_n[j];
	if(minX_i_n>X_i_n[j])minX_i_n=X_i_n[j];
	}
R=maxX_i_n-minX_i_n;
S=0;
for(i=0;i<n_;i++)
S+=(X[i]-xmean)*(X[i]-xmean);
S/=n_;
S=sqrt((double)(S));

RS[n_]=R/S;
}

for(i=0;i<N;i++) 
	{
	RS[i]=log(RS[i]);
	Y[i]=log((double)(i));
	}

ParRegres(N, RS, Y, &b0, &b1);
printf("deb: N=%i\n",N);
for(i=0;i<N;i++) printf("%i ",i);printf("===\n");
for(i=0;i<N;i++) printf("%f ",RS[i]);printf("===\n");
for(i=0;i<N;i++) printf("%f ",Y[i]); printf("===\n");
printf("deb: b0=%f, b1=%f\n",b0,b1);
return(b0);
}

int prepHurst(int N, byte *buffATGC)
{
int i;
for(i=0;i<N;i++)
X[i]=maskRY[(int)(buffATGC[i])];
return(0);
}


double propHurst(int N, byte *buffATGC)
{
double tmp=-1;
if(N<=1) return(-1);
prepHurst(N,buffATGC);
tmp=expHurst(N,X);

return(tmp);
}
