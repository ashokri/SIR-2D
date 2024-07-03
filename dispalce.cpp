// quasi_1D_final.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"

#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20
#define BIGNO 1.0e16
#define BIGNI 1.0e-9
#define XMIN 2.0
#define EPS 1.0e-16
#define EPS2 3.0e-11
#define FPMIN 1.0e-30
#define PI 3.141592653589793
#define MAXIT 10000
#define NUSE1 7
#define NUSE2 8
#define RTPIO2 1.2533141
#define freq 0.002
#define alpha 7.00
#define jm 14

typedef struct { double r, i; } doublecomplex;

double factorial(int n);
double sign(int m);
double z2_norm(doublecomplex a);
double reflection(double E,doublecomplex **f);
double transmition(double E,doublecomplex **f);
double incident(double E);
double Refle_rate(double E,doublecomplex **f);
double Trans_rate(double E,doublecomplex **f);
double plgndr(int l,int mm,double x);
double potential(double r, double x,double V0,double displace,double ph);
double PhaseShift(double kin,int am,int n,double rm,double V0);
double chebev(double a, double b, double c[], int m, double x);
double jl(int n,double x); 
double nl(int n,double x);
double phi(int n, double r,double x,double ph,double displace);
double snrm(int n,double sx[],int itol);
void ludcmp(double **H,double **aa,int n,int *indx,double *vv,double **a);
void LUslme(double **C,int m,int n,double **A,int *indx,double **B);
double **D_Matrix_alloc(int m,int n);
double Ta(int i,int n,double rm,double k,double x,double ph,double displace);
double Ts(int i,int n,double rm,double k,double x,double ph,double displace);
double delta(int i,int j);
double SIGN(double x,double y);
double pl(int l,int m,double x);

doublecomplex complex_add(doublecomplex a, doublecomplex b);
doublecomplex complex_mult(doublecomplex a, doublecomplex b);
doublecomplex scaler_mult(double a, doublecomplex b);
doublecomplex complex_sub(doublecomplex a, doublecomplex b);
doublecomplex complex_div(doublecomplex a, doublecomplex b);
doublecomplex *doublecomplexMalloc(int n);
doublecomplex **DCom_Matrix_alloc(int m,int n);
doublecomplex **LUc_slme(int m,int n,doublecomplex **A,doublecomplex **B);
doublecomplex a1(int i,int n,double rm,double x,double k,double ph,double displace);
doublecomplex T1(int i,int n,double rm,double x,double k0,double k2,double k4,double k6);
doublecomplex gama_sym(int i,int n,double rm,double x,double w,double k,double ph,double displace);
doublecomplex gama_asy(int i,int n,double rm,double x,double w,double k,double ph,double displace);
doublecomplex **F(double E,int m,int n,int n_t,int n_phi,double x,double w,doublecomplex *ss,doublecomplex *sa,double rm,double ph,double displace);


void initialize(double **A,int m,int n);
void construction(double thresh,double *a,int *ija,int n, int n_t, int nmax, double E,double rm,double x[],double w[],double V0,int *nnz);
void gauleg(double x1, double x2, double x[], double w[], int n);
void av(doublecomplex *v,doublecomplex *w,doublecomplex *a,int *asub,int *xa,int n);
void L2_construct(double L2[],int n_t,int n_phi, double x[],double w[],double ww[],double ph[]);
void SphBes(int n,double x,double *sj,double *sy,double *sjp,double *syp);
void bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp); 
void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
void DMUT(double V0,double E,double rm,double x[],int n,int n_t,doublecomplex **f,FILE *output,double ScatterinLenght3D,double Vp);
void solve(int n,int n_t,int n_phi,double ***C1,double ***C2,double ***C3,doublecomplex *ss,doublecomplex *sa,double E,double rm,double L2[],double x[],double w[],double ph[],double displace);
void asolve(int n,double b[],double x[],int itrnsp,double a[], int ija[]);
int *intMalloc(int n);
void *LU_MALLOC(size_t size);
void LU_FREE(void *addr);
void C_Construction(double ***C1,double ***C2,double ***C3,int n,int n_t,int n_phi,double En,double rm,double L2[],double x[],double w[],double V0,double displace,double ph[]);

double CurentConservation(double E,doublecomplex **f);
int IMAX(int i,int j);

int main()
{

	FILE /**input,*/*output;
		
	double k;
	double r;
	double E;
	double rm;
	double V0;
	double max;
	double Dif;
	double ScatterinLenght3D,Vp;
	double tol=0.0001;
	double errp;
	double invn_phi;
	double *phy;
	double *x;
	double *w;
	double *L2;
        double *err;
	double ***C1;
	double ***C2;
	double ***C3;
	
    	
	int itol=2;
	int itmax=100000;
	int infinity;
	int N;
	int n;
	int n_t;
	int n_phi;
	int l;
	int i;
	int ii;
	int conti = 0;
        int *indx;
	int am=0;
	char which;
	double *ww;
	double *ph;
	double d_on_a_0;
	double displace;
		
	doublecomplex zero={0.0,0.0};
	doublecomplex **f;
	doublecomplex *ss;
	doublecomplex *sa;

	err = &errp;
		
	/*if((input = fopen("input.txt","r")) == NULL)
	{
		printf("  can't open input file\n");
		return 0;
	}*/

	if((output = fopen("output.txt","w")) == NULL)
	{
		printf("  can't open output file\n");
		return 0;
	}
	
	/*fscanf(input, "%lf\n", &E);  // E=k*k/2 
	fscanf(input, "%c\n", &which);
	fscanf(input, "%lf\n", &max);
	fscanf(input, "%lf\n", &Dif);
	fscanf(input, "%d\n", &n);
	fscanf(input, "%d\n", &n_t);
        fscanf(input, "%d\n", &n_phi);
	fscanf(input, "%lf\n", &rm);*/
    
	E=2.0e-5;
	which='P';
  	V0 = 0.923;
	max=0.925;
	Dif=0.0005;
	n=50000;
	n_t=30;
	n_phi=9;
	rm=150.0;
	d_on_a_0=0.1;
	
   	displace=d_on_a_0/sqrt(freq);
	double P_devide=(2.0*PI)/(n_phi+0.0);
	
	N =  n*n_t*n_phi;
	
	E += freq;

	if( !(L2 = (double *) LU_MALLOC(n_t*n_phi*n_t*n_phi*sizeof(double))) ) printf("  LU_MALLOC fails for L2[].");
        if( !(x = (double *) LU_MALLOC(n_t*sizeof(double))) ) printf("  LU_MALLOC fails for x[].");
        if( !(w = (double *) LU_MALLOC(n_t*sizeof(double))) ) printf("  LU_MALLOC fails for w[].");
	if( !(ww = (double *) LU_MALLOC(n_t*sizeof(double))) ) printf("  LU_MALLOC fails for ww[].");
        if( !(indx = intMalloc(2*n_t*n_phi)) ) printf("  Malloc fails for indx.");
        if(!(ph = (double *) LU_MALLOC(n_phi*sizeof(double))) ) printf("  LU_MALLOC fails for ph[].");
		
  	if( !(ss = doublecomplexMalloc(7*n_t*n_phi)) ) printf("  doublecomplexMalloc fails for s[].");
	if( !(sa = doublecomplexMalloc(7*n_t*n_phi)) ) printf("  doublecomplexMalloc fails for s[].");
    	
	/***************************** CONSTRUCTION OF MATRICES A *********************************/
	
       double fi=(2.0*PI)/(n_phi+0.0);
       for(i=0;i<n_phi;i++){ph[i]=i*fi;}


	gauleg(-1,+1,x,w,n_t);
	for (i=0;i<n_t;i++){ww[i]=w[i]*P_devide;}
	printf("  size of matrix A = %d x %d\n", N,N);
			
	printf("  E = %e   \n",E);
	printf("  z_max = %e\n",rm*x[n_t-1]);
	printf("  \n\n");
		
	fprintf(output,"  E = %e\n  freq= %e\n  rm = %e\n  n = %d\n  n_t = %d\n  n_phi = %d\n\n  alpha = %e\n d/a_0=%e \n",E,freq,rm,n,n_t,n_phi,alpha,d_on_a_0);
	fprintf(output,"  -------------------------------------------------------------------------------------\n");
        
	L2_construct(L2,n_t,n_phi,x,w,ww,ph);

	k = 1.0e-16;//sqrt(2.0*(E-freq));
	
    
	
        C1 = (double ***)malloc( n * sizeof(double *) );
	C2 = (double ***)malloc( n * sizeof(double *) );
	C3 = (double ***)malloc( n * sizeof(double *) );
	
	for (i=0;i<10;i++)
	{ 
		if( !(C1[i] = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for C1[i][][].");
		if( !(C2[i] = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for C2[i][][].");
		if( !(C3[i] = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for C3[i][][].");
	}

	do
	{		
	
		ScatterinLenght3D = -PhaseShift(k,am,n,rm,2.0*V0)/k;
	    Vp = -PhaseShift(k,1,n,rm,2.0*V0)/(k*k*k);
	
		/********************   SOLVE THE LINEAR (SYMMETRIC) SYSTEMS PROBLEM   **********************/ 	
	
		C_Construction(C1,C2,C3,n,n_t,n_phi,E,rm,L2,x,w,2.0*V0,displace,ph);
		solve(n,n_t,n_phi,C1,C2,C3,ss,sa,E,rm,L2,x,w,ph,displace);	
		f = F(E,n,n,n_t,n_phi,x[n_t-1],w[n_t-1],ss,sa,rm,ph[n_phi-1],displace);
			
		DMUT(V0,E,rm,x,n,n_t,f,output,ScatterinLenght3D,Vp);
		
		/*for (i=1;i<8;i++) 
		{	
			r= rm*(exp(alpha*(n-7+i)/n)-1.0)/(exp(alpha)-1.0);;
			fprintf(output,"  %1.5e  %1.14e \n",r,ss[i*n_t-1]);
		}*/
	 
		if (which == 'P')
		    {
	              V0 = V0 + Dif;
                      if (V0 > max) conti = 1;
		    }
		else 
		   {
		     E = E + Dif;
		     if (E >= max) conti = 1;
		   }

	} while (conti == 0);
	
	free (C1);
	free (C2);
	free (C3);
	free (x);
	free (w);
	free (ww);
	free (ss);
	free (sa);
	free (indx);
	free (L2);
	free (ph);

	fclose(output);
	return 0;
}

void DMUT(double V0,double E,double rm,double x[],int n,int n_t,doublecomplex **f,FILE *output,double ScatterinLenght3D,double Vp)
{
	double T,g_1D;
	int i,j;
	doublecomplex **fe,**fo;
	
	if( !(fe = DCom_Matrix_alloc(4,2))) printf("  DMUT: Malloc fails for fe[][].");
	if( !(fo = DCom_Matrix_alloc(4,2))) printf("  DMUT: Malloc fails for fo[][].");

	printf("  \n");
   
	printf("  E = %e\n  s-scatteringlength = %e\n  p-scatteringvolume = %e\n",E,ScatterinLenght3D,Vp);		
	fprintf(output,"  %1.5e  %1.14e  %1.14e  %1.14e",V0,E,ScatterinLenght3D,Vp);		

	printf("  \n");

	printf("  f_e = %e   ",f[0][0].r);
	fprintf(output,"  %1.12e",f[0][0].r);
	
	printf("%e\n",f[0][0].i);
	fprintf(output," %1.12e",f[0][0].i);
		
	printf("  f_o = %e   ",f[0][1].r);
	fprintf(output,"  %1.12e",f[0][1].r);
	
	printf("%e\n",f[0][1].i);
	fprintf(output," %1.12e",f[0][1].i);
		
	fe[0][0]=f[0][0];
	fo[0][1]=f[0][1];
	
    g_1D = sqrt(2.0*(E-freq))*f[0][0].r/f[0][0].i;
	fprintf(output,"    %1.14e  %1.14e  ",g_1D,Trans_rate(E,fe));	

	printf("  Bosons:\n");
	printf("  Current Conservation (percentage) = %e \n",CurentConservation(E,fe));
	printf("  Transmition Rate = %f \n",Trans_rate(E,fe));
	printf("  Reflection Rate = %f \n",Refle_rate(E,fe));

	printf("  Fermions:\n");
	printf("  Current Conservation (percentage) = %e \n",CurentConservation(E,fo));
	printf("  Transmition Rate = %f \n",Trans_rate(E,fo));
	printf("  Reflection Rate = %f \n",Refle_rate(E,fo));

	printf("  Mix:\n");
	printf("  Current Conservation (percentage) = %e \n",CurentConservation(E,f));
	printf("  Transmition Rate = %f \n",Trans_rate(E,f));
	printf("  Reflection Rate = %f \n",Refle_rate(E,f));

	T=(1.0+f[0][0].r/2.0)*(1.0+f[0][0].r/2.0)+f[0][0].i*f[0][0].i/4.0;
	fprintf(output,"  b %1.4e",T);
		
	T=(1.0+f[0][1].r/2.0)*(1.0+f[0][1].r/2.0)+f[0][1].i*f[0][1].i/4.0;
	fprintf(output,"  f %1.4e",T);
		
	T=(1.0+f[0][0].r/2.0+f[0][1].r/2.0)*(1.0+f[0][0].r/2.0+f[0][1].r/2.0)+(f[0][0].i+f[0][1].i)*(f[0][0].i+f[0][1].i)/4.0;
	fprintf(output,"  m %1.4e\n",T);
		
}

double potential(double r,double x,double V0,double displace,double ph)
{
	return (-V0*exp(-r)/r + freq*freq*(r*r*(1.0-x*x)-2*r*displace*cos(ph)*sqrt(1.0-x*x)+displace*displace));
}

double potential0(double r,double V0)
{
	return(-V0*exp(-r)/r);
}


void C_Construction(double ***C1,double ***C2,double ***C3,int n,int n_t,int n_phi,double En,double rm,double L2[],double x[],double w[],double V0,double displace,double ph[])

{
	int j,jj,j1,j2,j3,I,II,J,JJ,K,KK;
	int *indx;
	double X,Y,O,P,Q,h=1.0/n,r0,r1,r2,r,F0,F1,F2;
	long double V;
	double **a00,**a11,**a22,**ajj,**A,**B,**C,**D,**E,**F,**G,**H,*vv,**a;
        
	
	double Co0, A0=(4.0/(360.0*h*h)+alpha*12.0/(720.0*h)),
		   Co1, A1=-(54.0/(360.0*h*h)+alpha*108.0/(720.0*h)),
		   Co2, A2=(536.0/(360.0*h*h)+alpha*528.0/(720*h)),
		   Co3, A3=(540.0/(360.0*h*h)+alpha*540.0/(720.0*h)),
		   Co4, A4=(-926.0/(360.0*h*h)+alpha*108.0/(720.0*h)),
		   Co5, 
		   Co6, A6=(536.0/(360.0*h*h)-alpha*552.0/(720.0*h)),
		   Co7, A7=(540.0/(360.0*h*h)-alpha*540.0/(720.0*h)),
		   Co8, A8=-(54.0/(360.0*h*h)-alpha*108.0/(720.0*h)),
		   Co9, A9=(4.0/(360.0*h*h)-alpha*12.0/(720.0*h));


        r0=rm*(exp(alpha*(0+1.0)/n)-1.0)/(exp(alpha)-1.0);
	r1=rm*(exp(alpha*(1+1.0)/n)-1.0)/(exp(alpha)-1.0);
	r2=rm*(exp(alpha*(2+1.0)/n)-1.0)/(exp(alpha)-1.0);
	F0=(exp(alpha)-1.0)/(alpha*(rm+(exp(alpha)-1.0)*r0)); 
	F1=(exp(alpha)-1.0)/(alpha*(rm+(exp(alpha)-1.0)*r1));
	F2=(exp(alpha)-1.0)/(alpha*(rm+(exp(alpha)-1.0)*r2));

	if( !(a00 = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for a00[][].");
	if( !(a11 = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for a11[][].");
	if( !(a22 = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for a22[][].");	
	if( !(indx = intMalloc(n_t*n_phi))) printf("  Malloc fails for indx.");	

	if( !(A = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for C[][].");
	if( !(B = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for B[][].");
	if( !(C = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for C[][].");
	
	if( !(H = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for H[][].");

	if( !(vv = (double *) LU_MALLOC(n_t*n_phi*sizeof(double))) ) printf("  ludcmp: LU_MALLOC fails for vv[].");
	if( !(a = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  ludcmp: Malloc fails for a[][].");

	//----------------------------------  C_0  ----------------------------------------//


	Co4=F0*F0*A4;
	Co6=F0*F0*A6;
	Co8=F0*F0*A8;
	Co9=F0*F0*A9;

	for (I=0;I<n_t;I++)
	     {
	      for (II=0;II<n_phi;II++)
	           {
	  	    A[I*n_phi+II][I*n_phi+II]=Co6;
		    B[I*n_phi+II][I*n_phi+II]=Co8; 
		    C[I*n_phi+II][I*n_phi+II]=Co9;
	           }
	     }

	 for(I=0;I<n_t;I++)
		{
	      for(II=0;II<n_phi;II++)
		  {
	     	   for(J=0;J<n_t;J++)
		        {
		          for(JJ=0;JJ<n_phi;JJ++)
		               {

				a00[I*n_phi+II][J*n_phi+JJ]= -L2[I*n_phi*n_t*n_phi+II*n_t*n_phi+J*n_phi+JJ]/(r0*r0);
				a11[I*n_phi+II][J*n_phi+JJ]= -L2[I*n_phi*n_t*n_phi+II*n_t*n_phi+J*n_phi+JJ]/(r1*r1);
				a22[I*n_phi+II][J*n_phi+JJ]= -L2[I*n_phi*n_t*n_phi+II*n_t*n_phi+J*n_phi+JJ]/(r2*r2);
			
				if (I == J && II==JJ)  
				    {
				     V= potential(r0,x[I],V0,displace,ph[II]);
				     a00[I*n_phi+II][J*n_phi+JJ] += Co4 + (2*En-V);
												
				     V= potential(r1,x[I],V0,displace,ph[II]);
				     a11[I*n_phi+II][J*n_phi+JJ] += -F1*F1*980.0/(360.0*h*h) + (2*En-V); 

				     V= potential(r2,x[I],V0,displace,ph[II]);
				     a22[I*n_phi+II][J*n_phi+JJ] += -F2*F2*980.0/(360.0*h*h)+ (2*En-V); 
			            }	
		               }
	                }
                      }
		   }
	ludcmp(H,a00,n_t*n_phi,indx,vv,a);

	LUslme(C1[0],n_t*n_phi,n_t*n_phi,H,indx,A);
	LUslme(C2[0],n_t*n_phi,n_t*n_phi,H,indx,B);
	LUslme(C3[0],n_t*n_phi,n_t*n_phi,H,indx,C);

	
	//-----------------------------------  C_1  ---------------------------------------//		

	
	Co2=F1*F1*A2;
	Co7=F1*F1*A7;
	Co8=F1*F1*A8;
	Co9=F1*F1*A9;

	initialize(A,n_t*n_phi,n_t*n_phi);
	initialize(B,n_t*n_phi,n_t*n_phi);
	initialize(C,n_t*n_phi,n_t*n_phi);
		
	if( !(D = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for D[][].");

	for (I=0;I<n_t;I++)
	     {
	      for(II=0;II<n_phi;II++)
 	          {
		   for (J=0;J<n_t;J++)
		        {
		         for (JJ=0;JJ<n_phi;JJ++)
		              {		  
				D[I*n_phi+II][J*n_phi+JJ] = -(a11[I*n_phi+II][J*n_phi+JJ]+Co2*C1[0][I*n_phi+II][J*n_phi+JJ]);
				A[I*n_phi+II][J*n_phi+JJ] = Co2*C2[0][I*n_phi+II][J*n_phi+JJ];
				B[I*n_phi+II][J*n_phi+JJ] = Co2*C3[0][I*n_phi+II][J*n_phi+JJ];
		              }
		        }
                   A[I*n_phi+II][I*n_phi+II] += Co7;
	           B[I*n_phi+II][I*n_phi+II] += Co8;
	           C[I*n_phi+II][I*n_phi+II] = Co9;
	          }
	     }

	ludcmp(H,D,n_t*n_phi,indx,vv,a);

	LUslme(C1[1],n_t*n_phi,n_t*n_phi,H,indx,A);
	LUslme(C2[1],n_t*n_phi,n_t*n_phi,H,indx,B);
	LUslme(C3[1],n_t*n_phi,n_t*n_phi,H,indx,C);

	
	//-------------------------------------  C_2  --------------------------------------------//

	Co1=F2*F2*A1;
	Co3=F2*F2*A3;
	Co7=F2*F2*A7;
	Co8=F2*F2*A8;
	Co9=F2*F2*A9;

	initialize(A,n_t*n_phi,n_t*n_phi);
	initialize(B,n_t*n_phi,n_t*n_phi);
	initialize(C,n_t*n_phi,n_t*n_phi);
	initialize(D,n_t*n_phi,n_t*n_phi);
	
	for (I=0;I<n_t;I++)
	      {
		for (II=0;II<n_phi;II++)
		     {
		      for (J=0;J<n_t;J++)
			   {
			    for (JJ=0;JJ<n_phi;JJ++)
				 {
				  for (K=0;K<n_t;K++)
				       {
			                for (KK=0;KK<n_phi;KK++)
			                     {
					      D[I*n_phi+II][J*n_phi+JJ] -= Co1*C1[0][I*n_phi+II][K*n_phi+KK]*C1[1][K*n_phi+KK][J*n_phi+JJ];
					      A[I*n_phi+II][J*n_phi+JJ] += Co1*C1[0][I*n_phi+II][K*n_phi+KK]*C2[1][K*n_phi+KK][J*n_phi+JJ];
					      B[I*n_phi+II][J*n_phi+JJ] += Co1*C1[0][I*n_phi+II][K*n_phi+KK]*C3[1][K*n_phi+KK][J*n_phi+JJ];
					     }
				       }
				  D[I*n_phi+II][J*n_phi+JJ] -= Co1*C2[0][I*n_phi+II][J*n_phi+JJ]+Co3*C1[1][I*n_phi+II][J*n_phi+JJ]+a22[I*n_phi+II][J*n_phi+JJ];
				  A[I*n_phi+II][J*n_phi+JJ] += Co1*C3[0][I*n_phi+II][J*n_phi+JJ]+Co3*C2[1][I*n_phi+II][J*n_phi+JJ];
				  B[I*n_phi+II][J*n_phi+JJ] += Co3*C3[1][I*n_phi+II][J*n_phi+JJ];
				 }
			   }
		      A[I*n_phi+II][I*n_phi+II] += Co7;
		      B[I*n_phi+II][I*n_phi+II] += Co8;
		      C[I*n_phi+II][I*n_phi+II] += Co9;
		     }
	      }

	ludcmp(H,D,n_t*n_phi,indx,vv,a);

	LUslme(C1[2],n_t*n_phi,n_t*n_phi,H,indx,A);
	LUslme(C2[2],n_t*n_phi,n_t*n_phi,H,indx,B);
	LUslme(C3[2],n_t*n_phi,n_t*n_phi,H,indx,C);
	
	free(a00);
	free(a11);
	free(a22);
	

	//-------------------------------------  C_I  --------------------------------------------//

	if( !(ajj = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for ajj[][].");

	if( !(E = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for E[][].");
        if( !(F = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for F[][].");
	if( !(G = D_Matrix_alloc(n_t*n_phi,n_t*n_phi))) printf("  Malloc fails for G[][].");


	j3=0;
	j2=1;
	j1=2;	
	jj=3;

	for (j=3;j<n-3;j++)
	     {
	      r0=rm*(exp(alpha*(j+1.0)/n)-1.0)/(exp(alpha)-1.0);
	      F0=(exp(alpha)-1.0)/(alpha*(rm+(exp(alpha)-1.0)*r0));
	
	      Co0=F0*F0*A0;
	      Co1=F0*F0*A1;
	      Co3=F0*F0*A3;
	      Co7=F0*F0*A7;
	      Co8=F0*F0*A8;
	      Co9=F0*F0*A9;

		for (I=0;I<n_t;I++)
	             {
		     for (II=0;II<n_phi;II++)
			  {
		          for (J=0;J<n_t;J++)
			       {
			       for (JJ=0;JJ<n_phi;JJ++)
			            {	
			             ajj[I*n_phi+II][J*n_phi+JJ] = -L2[I*n_phi*n_t*n_phi+II*n_t*n_phi+J*n_phi+JJ]/(r0*r0);
			             if (I == J && II==JJ)  
			                 {					
			                  V= potential(r0,x[I],V0,displace,ph[II]);
			                  ajj[I*n_phi+II][J*n_phi+JJ] += -F0*F0*980.0/(360.0*h*h)+(2*En-V);  
			                 }
            		            }
			        }
			     }
                           }	      		     		
		initialize(A,n_t*n_phi,n_t*n_phi);
		initialize(B,n_t*n_phi,n_t*n_phi);
		initialize(C,n_t*n_phi,n_t*n_phi);
		initialize(D,n_t*n_phi,n_t*n_phi);
		initialize(E,n_t*n_phi,n_t*n_phi);
		initialize(F,n_t*n_phi,n_t*n_phi);
		initialize(G,n_t*n_phi,n_t*n_phi);
	             
		for (I=0;I<n_t;I++)
		     {
		      for (II=0;II<n_phi;II++)
		           {
			    for (J=0;J<n_t;J++)
			         {
			          for (JJ=0;JJ<n_phi;JJ++)
			               {
				        X = 0;
				        Y = 0;
				        O = 0;
				        P = 0;	
				        Q = 0;

				        for (K=0;K<n_t;K++)
				             {
				              for (KK=0;KK<n_phi;KK++)
				                   {
					            X += C1[j3][I*n_phi+II][K*n_phi+KK]*C2[j2][K*n_phi+KK][J*n_phi+JJ];
					            Y += C2[j3][I*n_phi+II][K*n_phi+KK]*C1[j1][K*n_phi+KK][J*n_phi+JJ];
					            O += C1[j3][I*n_phi+II][K*n_phi+KK]*C3[j2][K*n_phi+KK][J*n_phi+JJ];
						    P += C2[j3][I*n_phi+II][K*n_phi+KK]*C2[j1][K*n_phi+KK][J*n_phi+JJ];	
						    Q += C2[j3][I*n_phi+II][K*n_phi+KK]*C3[j1][K*n_phi+KK][J*n_phi+JJ];

						    E[I*n_phi+II][J*n_phi+JJ] += C1[j2][I*n_phi+II][K*n_phi+KK]*C1[j1][K*n_phi+KK][J*n_phi+JJ];
						    F[I*n_phi+II][J*n_phi+JJ] += C1[j2][I*n_phi+II][K*n_phi+KK]*C2[j1][K*n_phi+KK][J*n_phi+JJ];
						    G[I*n_phi+II][J*n_phi+JJ] += C1[j2][I*n_phi+II][K*n_phi+KK]*C3[j1][K*n_phi+KK][J*n_phi+JJ];
						    }
					      }

					D[I*n_phi+II][J*n_phi+JJ] = Co0*(X+Y)+Co1*E[I*n_phi+II][J*n_phi+JJ]+Co3*C1[j1][I*n_phi+II][J*n_phi+JJ]+Co1*C2[j2][I*n_phi+II][J*n_phi+JJ]+Co0*C3[j3][I*n_phi+II][J*n_phi+JJ]+ajj[I*n_phi+II][J*n_phi+JJ];
					A[I*n_phi+II][J*n_phi+JJ] = -(Co0*(O+P)+Co1*F[I*n_phi+II][J*n_phi+JJ]+Co3*C2[j1][I*n_phi+II][J*n_phi+JJ]+Co1*C3[j2][I*n_phi+II][J*n_phi+JJ]);
					B[I*n_phi+II][J*n_phi+JJ] = -(Co0*Q+Co1*G[I*n_phi+II][J*n_phi+JJ]+Co3*C3[j1][I*n_phi+II][J*n_phi+JJ]);
				       }
		                 }

			    A[I*n_phi+II][I*n_phi+II] -= Co7;
			    B[I*n_phi+II][I*n_phi+II] -= Co8;
		           }
	             }
	             
		for (I=0;I<n_t;I++)	
		     {
		      for (II=0;II<n_phi;II++)
			   {
			    for (J=0;J<n_t;J++)
			         {
			          for (JJ=0;JJ<n_phi;JJ++)
				       {
				        for (K=0;K<n_t;K++)
				             {
				              for (KK=0;KK<n_phi;KK++)
				                   {
					             D[I*n_phi+II][J*n_phi+JJ] += Co0*C1[j3][I*n_phi+II][K*n_phi+KK]*E[K*n_phi+KK][J*n_phi+JJ];
					             A[I*n_phi+II][J*n_phi+JJ] -= Co0*C1[j3][I*n_phi+II][K*n_phi+KK]*F[K*n_phi+KK][J*n_phi+JJ];
					 	     B[I*n_phi+II][J*n_phi+JJ] -= Co0*C1[j3][I*n_phi+II][K*n_phi+KK]*G[K*n_phi+KK][J*n_phi+JJ];
				                   }
				             }
				       }
			         }

			    C[I*n_phi+II][I*n_phi+II] = -Co9;
		           }
		     }

		ludcmp(H,D,n_t*n_phi,indx,vv,a);

		LUslme(C1[jj],n_t*n_phi,n_t*n_phi,H,indx,A);
		LUslme(C2[jj],n_t*n_phi,n_t*n_phi,H,indx,B);
		LUslme(C3[jj],n_t*n_phi,n_t*n_phi,H,indx,C);

		j3++;if(j3>9) j3=0;
		j2++;if(j2>9) j2=0;
		j1++;if(j1>9) j1=0;	
		jj++;if(jj>9) jj=0;	
		
	      }

	free(ajj);
	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(G);
	free(indx);
	free(H);
	free(vv);
	free(a);
}

void solve(int n,int n_t,int n_phi,double ***C1,double ***C2,double ***C3,doublecomplex *ss,doublecomplex *sa,double E,double rm,double L2[],double x[],double w[],double ph[],double displace)
{	
	doublecomplex **A,**B,**S,one={1.0,0.0};
	int i,I,II,J,JJ;
	double k;

	k=sqrt(2*(E-freq));  
	
	if( !(A = DCom_Matrix_alloc(7*n_t*n_phi,7*n_t*n_phi))) printf("  solve routine: Malloc fails for A[][].");
	if( !(B = DCom_Matrix_alloc(7*n_t*n_phi,2))) printf("  solve routine: Malloc fails for B[][].");
	if( !(S = DCom_Matrix_alloc(7*n_t*n_phi,2))) printf("  solve routine: Malloc fails for S[][].");

	for (i=0;i<4;i++)
		{
	     for (I=0;I<n_t;I++)
		  {
		    for (II=0;II<n_phi;II++)
		         {
			  A[i*n_t*n_phi+I*n_phi+II][i*n_t*n_phi+I*n_phi+II].r=-1.0;

			  for (J=0;J<n_t;J++)
			       {
				for (JJ=0;JJ<n_phi;JJ++)
				     {
				  
				       A[i*n_t*n_phi+I*n_phi+II][i*n_t*n_phi+n_t*n_phi+J*n_phi+JJ].r= C1[10-7+i][I*n_phi+II][J*n_phi+JJ];
				       A[i*n_t*n_phi+I*n_phi+II][i*n_t*n_phi+2*n_t*n_phi+J*n_phi+JJ].r= C2[10-7+i][I*n_phi+II][J*n_phi+JJ];
				       A[i*n_t*n_phi+I*n_phi+II][i*n_t*n_phi+3*n_t*n_phi+J*n_phi+JJ].r=C3[10-7+i][I*n_phi+II][J*n_phi+JJ];
				     }
			       }
		         }
		  }
		}
	for (i=4;i<7;i++)
		{
	      for (I=0;I<n_t;I++)
		    {
		      for (II=0;II<n_phi;II++)
		            {
			      A[i*n_t*n_phi+I*n_phi+II][i*n_t*n_phi+I*n_phi+II].r=1.0;
			
			      /*A[i*n_t+I][(i-4)*n_t+I]=a4(n-7+i,n,rm,x[I],k);
			      A[i*n_t+I][(i-3)*n_t+I]=a3(n-7+i,n,rm,x[I],k);
			      A[i*n_t+I][(i-2)*n_t+I]=a2(n-7+i,n,rm,x[I],k);*/
			      
			      A[i*n_t*n_phi+I*n_phi+II][(i-1)*n_t*n_phi+I*n_phi+II]=a1(n-7+i,n,rm,x[I],k,ph[II],displace);
			      B[i*n_t*n_phi+I*n_phi+II][0]=gama_sym(n-7+i,n,rm,x[I],w[I],k,ph[II],displace);
			      B[i*n_t*n_phi+I*n_phi+II][1]=gama_asy(n-7+i,n,rm,x[I],w[I],k,ph[II],displace);
			      
			    }
		    }
		}
	S = LUc_slme(7*n_t*n_phi,2,A,B);

	for (i=0;i<7;i++)
		{
		for (I=0;I<n_t;I++)
		     {
		      for (II=0;II<n_phi;II++)
		           {
			    ss[i*n_t*n_phi+I*n_phi+II] = S[i*n_t*n_phi+I*n_phi+II][0];
			    sa[i*n_t*n_phi+I*n_phi+II] = S[i*n_t*n_phi+I*n_phi+II][1];
		           }
		     }
 		}
	free(A);
	free(B);
	free(S);
}


double Trans_rate(double E,doublecomplex **f)
{
	return (transmition(E,f)/incident(E));
}

double Refle_rate(double E,doublecomplex **f)
{
	return (reflection(E,f)/incident(E));
}

double CurentConservation(double E,doublecomplex **f)
{
	double In,Tr,Re;
	double out;

	In = incident(E);
	Tr = transmition(E,f);
	Re = reflection(E,f);

	out = (double) 100*((double)(Tr+Re)/In);
	return out;
}

double incident(double E)
{
	double a=0.0;
	
	a+=sqrt(2*(E-1.0*freq)); 
	
	return a;
}

double transmition(double E,doublecomplex **f)
{
	doublecomplex b={0.0,0.0};
	double a=0.0;
	int i;

	b.r=2.0; 
	a+=z2_norm(complex_add(b,complex_add(f[0][0],f[0][1])))*sqrt(2*(E-freq))/4.0;
			
	return (a);
}

double reflection(double E,doublecomplex **f)
{
	double a=0.0;
	int i;

	a+=z2_norm(complex_sub(f[0][0],f[0][1]))*sqrt(2*(E-freq))/4.0;
		
	return (a);
}

doublecomplex a1(int i,int n,double rm,double x,double k,double ph,double displace)
{
	doublecomplex a,b,zero={0.0,0.0},one={1.0,0.0};
	double y,r,rr,p;

	y = sqrt(1.0-x*x);
	
	a = one;//complex_add(one,T1(i,n,rm,x,k));

	r = rm*(exp(alpha*(double)(i+1.0)/n)-1.0)/(exp(alpha)-1.0);
	rr = rm*(exp(alpha*(double)i/n)-1.0)/(exp(alpha)-1.0);
	
	p = r*phi(0,r,x,ph,displace)/(rr*phi(0,rr,x,ph,displace));
	b.r = -p*cos(k*(r-rr)*fabs(x));
	b.i = -p*sin(k*(r-rr)*fabs(x));
	
	return(complex_mult(a,b));
}

doublecomplex gama_sym(int i,int n,double rm,double x,double w,double k,double ph,double displace)
{
	double y,r,rr,p;
	doublecomplex T,a,b,c,d,e,out={0.0,0.0};//,one={1.0,0.0};

	y = sqrt(1.0-x*x);

	r = rm*(exp(alpha*(double)(i+1.0)/n)-1.0)/(exp(alpha)-1.0);
	p = 2*cos(k*r*x);
	a.r = p;
	a.i = 0.0;

	rr = rm*(exp(alpha*(double)i/n)-1.0)/(exp(alpha)-1.0);
	p = 2*cos(k*rr*x);	
	b.r = -cos(k*(r-rr)*fabs(x))*p;
	b.i = -sin(k*(r-rr)*fabs(x))*p;
			
	out.r += sqrt(w)*r*phi(0,r,x,ph,displace)*(a.r+b.r);
	out.i += sqrt(w)*r*phi(0,r,x,ph,displace)*(a.i+b.i);
	
	return out;
}

doublecomplex gama_asy(int i,int n,double rm,double x,double w,double k,double ph,double displace)
{	
	double y,r,rr,p;
	doublecomplex T,a,b,c,d,e,out={0.0,0.0};//,one={1.0,0.0};

	y = sqrt(1.0-x*x);

	r = rm*(exp(alpha*(double)(i+1.0)/n)-1.0)/(exp(alpha)-1.0);
	p = 2*sin(k*r*x);
	a.r = p;
	a.i = 0.0;

	rr = rm*(exp(alpha*(double)i/n)-1.0)/(exp(alpha)-1.0);
	p = 2*sin(k*rr*x);	
	b.r = -cos(k*(r-rr)*fabs(x))*p;
	b.i = -sin(k*(r-rr)*fabs(x))*p;
			
	out.r -= sqrt(w)*r*phi(0,r,x,ph,displace)*(a.i+b.i);
	out.i += sqrt(w)*r*phi(0,r,x,ph,displace)*(a.r+b.r);
	 
	return out;
}


doublecomplex **F(double E,int m,int n,int n_t,int n_phi,double x,double w,doublecomplex *ss,doublecomplex *sa,double rm,double ph,double displace)
{
	doublecomplex **G,**z,**fe,**fo,**f;
	double k,g,r,rr,y=sqrt(1-x*x);
	int i,j;

	if( !(f = DCom_Matrix_alloc(4,2))) printf("  F: Malloc fails for f[][].");
	
	j=1;
	
	if( !(G = DCom_Matrix_alloc(j,j))) printf("  F: Malloc fails for G[][].");
	if( !(z = DCom_Matrix_alloc(j,1))) printf("  F: Malloc fails for z[][].");

	k=sqrt(2*(E-freq));

	for (i=0;i<j;i++)
	{
		g = Ts(m-1-i,n,rm,k,x,ph,displace);
		r = rm*(exp(alpha*(double)(m-i)/n)-1.0)/(exp(alpha)-1.0);
		rr = rm*(exp(alpha*(double)(m-i-1.0)/n)-1.0)/(exp(alpha)-1.0);
		
		G[i][0].r=(cos(k*r*x)*phi(0,r,x,ph,displace)-g*cos(k*rr*x)*phi(0,rr,x,ph,displace));
		G[i][0].i=(sin(k*r*x)*phi(0,r,x,ph,displace)-g*sin(k*rr*x)*phi(0,rr,x,ph,displace));
					    
		z[i][0]=complex_sub(scaler_mult(1.0/(sqrt(w)*r),ss[(7-i)*n_t*n_phi-1]),scaler_mult(g/(sqrt(w)*rr),ss[(7-i-1)*n_t*n_phi-1]));
	}

	f[0][0]=complex_div(z[0][0],G[0][0]);//fe=LUc_slme(j,1,G,z);

	for (i=0;i<j;i++)
	{
		g = Ta(m-1-i,n,rm,k,x,ph,displace);
		r = rm*(exp(alpha*(double)(m-i)/n)-1.0)/(exp(alpha)-1.0);
		rr = rm*(exp(alpha*(double)(m-i-1)/n)-1.0)/(exp(alpha)-1.0);

		G[i][0].r=cos(k*r*x)*phi(0,r,x,ph,displace)-g*cos(k*rr*x)*phi(0,rr,x,ph,displace);
		G[i][0].i=sin(k*r*x)*phi(0,r,x,ph,displace)-g*sin(k*rr*x)*phi(0,rr,x,ph,displace);
		
		z[i][0]=complex_sub(scaler_mult(1.0/(sqrt(w)*r),sa[(7-i)*n_t*n_phi-1]),scaler_mult(g/(sqrt(w)*rr),sa[(7-i-1)*n_t*n_phi-1]));
	}

	f[0][1]=complex_div(z[0][0],G[0][0]);//fo=LUc_slme(j,1,G,z);
	
	/*for (i=0;i<j;i++)
	{	
		f[i][0]=fe[i][0];
		f[i][1]=fo[i][0];
	}/**/
	free(G);
	free(z);
	//free(fo);
	//free(fe);

	return(f);

}


double Ts(int i,int n,double rm,double k,double x,double ph,double displace)
{
	double r,rr,y=sqrt(1.0-x*x),a=0.0,b=0.0;

	r = rm*(exp(alpha*(double)(i+1.0)/n)-1.0)/(exp(alpha)-1.0);
	rr = rm*(exp(alpha*(double)i/n)-1.0)/(exp(alpha)-1.0);

	a = cos(k*r*x)*phi(0,r,x,ph,displace); 
	b = cos(k*rr*x)*phi(0,rr,x,ph,displace); 

	return(a/b);
}


double Ta(int i,int n,double rm,double k,double x,double ph,double displace)
{
	double r,rr,y=sqrt(1.0-x*x),a=0.0,b=0.0;

	r = rm*(exp(alpha*(double)(i+1.0)/n)-1.0)/(exp(alpha)-1.0);
	rr = rm*(exp(alpha*(double)i/n)-1.0)/(exp(alpha)-1.0);

	a = sin(k*r*x)*phi(0,r,x,ph,displace); 
	b = sin(k*rr*x)*phi(0,rr,x,ph,displace); 

	return(a/b);
}



void L2_construct(double L2[],int n_t,int n_phi,double x[], double w[],double ww[],double ph[])
{
int Nt=n_t;// baze taghirate l dar chand jomleihaie legendre
int Nph=n_phi;

int M=(Nph-1)/2;
int L=Nt+M-1;
int m;

double P_inv=1.0/(2.0*PI);



//double ph[Nph];
int j,i,l,ll,k,I,J;                   //ll=l_prime
double pl_mad[2*M+1][L+1][Nt];
double pl_bar[2*M+1][L+1][Nt];
double sigma,norm,sqnorm,s;

/*double fi=(2.0*PI)/(n_phi+0.0);
       for(i=0;i<Nph;i++){ph[i]=i*fi;
		}
//gauleg(-1, 1, x, w,Nt);

         /*--------------------------tolide jomalate legndre be ezaie l=n-teta--------------------------------*/

  for(m=-M;m<=M;m++)
    {
      for(l=abs(m);l<Nt;l++)
          {
	    for(j=0;j<Nt;j++)
	        {
	         pl_bar[m+M][l][j] = pl(l,m,x[j]);
		}
	  }
    }
  

for(m=-M;m<=M;m++)
    {
     for(i=0;i<Nt;i++)
         {
	  sigma = 0.0;
	  if(m!=0)
             {
              for(ll=abs(m);ll<Nt;ll++)
                  {
	           s =0.0;
                   for(j=0;j<Nt;j++)
                       {
	                 s += pl(ll,m,x[j])*pl(Nt,m,x[j])*w[j];
                       }
                   sigma += s*pl(ll,m,x[i]);
                  } 
              pl_mad[m+M][Nt][i] = pl(Nt,m,x[i]) - sigma;
             }
	 }
    }

for(m=-M;m<=M;m++)
    {
     if(m!=0)
        {
         norm = 0.0;
         sqnorm = 0.0;
         for(i=0;i<Nt;i++)
             {
               norm += pl_mad[m+M][Nt][i]*pl_mad[m+M][Nt][i]*w[i];
             }
         sqnorm = sqrt(norm);
         for(j=0;j<Nt;j++)
             {
               pl_bar[m+M][Nt][j] =  pl_mad[m+M][Nt][j]/sqnorm;

	     }
        }
        else
	{
         for(j=0;j<Nt;j++)
             {
              pl_bar[M][Nt][j] = pl(Nt,0,x[j]);
             }
        }
    }

  

   //-----------------------------------------end of determination pbar for l=Nt-----------------------------------------

   //------------------------------------determination of pbar for upper than n_tetha l's--------------------------------------     

  
       for(m=-M;m<=M;m++)
           {
            for(l=(Nt+1);l<Nt+abs(m);l++)
                {
		  for(i=0;i<Nt;i++)
                      {
 	               sigma = 0.0;
                       for(ll=abs(m);ll<l;ll++)
                           {
		            if(ll<Nt)
                               {
			        s = 0.0;
		                for(j=0;j<Nt;j++)
                                    {
                                     s +=pl(l,m,x[j])*pl(ll,m,x[j])*w[j];
                                    }
                                sigma +=s*pl(ll,m,x[i]);
			       }
		            if(ll>=Nt)
                               {
			        s = 0.0;
			        for(j=0;j<Nt;j++)
                                    {
                                     s +=pl(l,m,x[j])*pl_bar[m+M][ll][j]*w[j];
                                    }
                                sigma +=s*pl_bar[m+M][ll][i];
                               }
                           }
                       pl_mad[m+M][l][i] = pl(l,m,x[i])-sigma;
//       	       printf("%e        %e\n",/*pmad[m+M][l][i]*/sigma,plgndr(l,m,x[i]));
                      }
                 norm = 0.0;
	         for(i=0;i<Nt;i++)
                     {
                      norm += pl_mad[m+M][l][i]*pl_mad[m+M][l][i]*w[i];
                     }
	         sqnorm =(1.0)/sqrt(norm);
	         for(j=0;j<Nt;j++)
                     {
                      pl_bar[m+M][l][j] =  pl_mad[m+M][l][j]*sqnorm;
//    		      printf("l=%d    m=%d    %f    %f\n",l,m,pbar[m+M][l][j],plgndr(l,m,x[j]));
                     }
                }
           }





		/*--------------------------------tolide vizhe maghadire L^2--------------------------*/

for(i=0;i<Nt;i++){
	for(j=0;j<Nph;j++){
		for(I=0;I<Nt;I++){
			for(J=0;J<Nph;J++){
				L2[i*Nt*Nph*Nph+j*Nph*Nt+I*Nph+J]=0.0;
			}}}}

double h;
for(i=0;i<Nt;i++){
	for(j=0;j<Nph;j++){
		for(I=0;I<Nt;I++){
			for(J=0;J<Nph;J++){
                                h=0.0;
				for(m=-M;m<=M;m++){
					
						for(l=abs(m);l<Nt+abs(m);l++){
							h+=P_inv*l*(l+1)*cos(m*(ph[j]-ph[J]))*pl_bar[m+M][l][i]*pl_bar[m+M][l][I]*sqrt(ww[i]*ww[I]);
						}
				}
			
			  L2[i*Nt*Nph*Nph+j*Nph*Nt+I*Nph+J]=h;
			} 
		}
	}	
}




 

}




double sign(int m)
{
    int i,s;
    s=1;
    for(i=1;i<=abs(m);i++)
         {
	   s *=-1;
         }
    return s;
}


double phi(int n, double r,double x,double ph,double displace)
{
	double X,Y,R,r_prime;
	R=r*sqrt(1.0-x*x);
	X=R*cos(ph);
	Y=R*sin(ph);
	r_prime=sqrt((X-displace)*(X-displace)+Y*Y);
	double out=0.0;
	
	if (n == 0)
		out = sqrt(2.0*freq)*exp(-freq*r_prime*r_prime/2.0);
	else if (n == 2)
		out = sqrt(2.0*freq)*exp(-freq*r_prime*r_prime/2.0)*(-freq*r_prime*r_prime + 1.0);
	else if (n == 4) 
		out = sqrt(2.0*freq)*exp(-freq*r_prime*r_prime/2.0)*(freq*freq*r_prime*r_prime*r_prime*r_prime/2.0 -2.0*freq*r_prime*r_prime + 1.0);
	else if (n == 6)
		out = sqrt(2.0*freq)*exp(-freq*r_prime*r_prime/2.0)*(-freq*freq*freq*r_prime*r_prime*r_prime*r_prime*r_prime*r_prime/6.0 + 3.0*freq*freq*r_prime*r_prime*r_prime*r_prime/2.0 -3.0*freq*r_prime*r_prime + 1.0);
	else printf("  error phi: n should be 0, 2, 4, or 6.\n");/**/
	
	return out;
}


void ludcmp(double **H,double **aa,int n,int *indx/*,double *d*/,double *vv,double **a)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			a[i][j]=aa[i][j];

	for (i=0;i<n;i++)
	{
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;

		if (big == 0.0) printf("  Singular matrix in routine ludcmp");

		vv[i]=1.0/big;
	}
	
	for (j=0;j<n;j++)
	{
		for (i=0;i<j;i++)
		{
			sum=a[i][j];
			for (k=0;k<i;k++)	sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}

		big=0.0;
		for (i=j;i<n;i++)
		{
			sum=a[i][j];
			for (k=0;k<j;k++)	sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if( (dum=vv[i]*fabs(sum)) >= big)
			{
				big=dum;
				imax=i;
			}
		}
		if (j != imax)
		{
			for (k=0;k<n;k++)
			{
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
//			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;

		if (j != n-1)
		{
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++)	a[i][j] *= dum;
		}
	}

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			H[i][j]=a[i][j];

 }

void lubksb(double **a, int n,int *indx,double b[])
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++)
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];

		if(ii>= 0)
			for (j=ii;j<=i-1;j++)	sum -= a[i][j]*b[j];
		else if (sum) ii=i;

		b[i]=sum;
	}

	for (i=n-1;i>=0;i--)
	{
		sum=b[i];

		for (j=i+1;j<n;j++)	sum -= a[i][j]*b[j];

		b[i]=sum/a[i][i];
	}
}

void LUslme(double **C,int m,int n,double **A,int *indx,double **B)	//  Solves a linear double matrix equation AX=B
																//  Here A is m x m.  X is m x n. B is m x n.
{
	int i,j;
	double *col;//,**X;
	
	//if( !(X = D_Matrix_alloc(m,n))) printf("  LUslme: LU_Malloc fails for X[][].");
	
	if( !(col = (double *) LU_MALLOC(m*sizeof(double))) ) printf("  LUslme: LU_MALLOC fails for col[].");
	
	for (j=0;j<n;j++)
	{
		for (i=0;i<m;i++)	col[i]=B[i][j];

		lubksb(A,m,indx,col);

		for (i=0;i<m;i++)	C[i][j]=col[i];
	}

	free(col);

	//return (X);
}

doublecomplex **LUc_slme(int m,int n,doublecomplex **A,doublecomplex **B)	
																//  Solves a linear doublecomplex matrix equation CX=B
																//  Here C is m x m.  X is m x n. B is m x n.
{
	int i,j,*indx;
	doublecomplex **X;
	double **x,**c,**b,**H,*vv,**a;
	
	if( !(X = DCom_Matrix_alloc(m,n))) printf("  LUc_slme: LU_Malloc fails for X[][].");
	if( !(c = D_Matrix_alloc(2*m,2*m))) printf("  LUc_slme: LU_Malloc fails for a[][].");
	if( !(H = D_Matrix_alloc(2*m,2*m))) printf("  LUc_slme: LU_Malloc fails for H[][].");
	if( !(b = D_Matrix_alloc(2*m,n))) printf("  LUc_slme: LU_Malloc fails for b[][].");
	if( !(x = D_Matrix_alloc(2*m,n))) printf("  LUc_slme: LU_Malloc fails for x[][].");
	if ( !(indx = intMalloc(2*m)) ) printf("  LUc_slme: Malloc fails for indx.");	
    if( !(vv = (double *) LU_MALLOC(2*m*sizeof(double))) ) printf("  ludcmp: LU_MALLOC fails for vv[].");
	if( !(a = D_Matrix_alloc(2*m,2*m))) printf("  ludcmp: Malloc fails for a[][].");

	for (i=0;i<m;i++)
		for (j=0;j<m;j++)
		{
			c[i][j]=A[i][j].r;
			c[i+m][j]=A[i][j].i;
			c[i][j+m]=-A[i][j].i;
			c[i+m][j+m]=A[i][j].r;
		}
	
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
		{
			b[i][j]=B[i][j].r;
			b[i+m][j]=B[i][j].i;
		}

	ludcmp(H,c,2*m,indx,vv,a);
    LUslme(x,2*m,n,H,indx,b);
			
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
		{
			X[i][j].r = x[i][j];
			X[i][j].i = x[i+m][j];
		}
			
	free(c);
	free(b);
	free(x);
	free(H);
	free(vv);
	free(a);

	return (X);
}

double z2_norm(doublecomplex a)
{
	return (a.r*a.r + a.i*a.i);
	
}

doublecomplex complex_add(doublecomplex a, doublecomplex b)
{
	doublecomplex out;

	out.r = a.r + b.r;
        out.i = a.i + b.i;
	
	return out;
	
}

doublecomplex scaler_mult(double a, doublecomplex b)
{
	doublecomplex out;

	out.r = a*b.r;
        out.i = a*b.i;
	
	return out;
	
}

doublecomplex complex_mult(doublecomplex a, doublecomplex b)
{
	doublecomplex out;
	
	out.r = a.r*b.r - a.i*b.i;
        out.i = a.i*b.r + a.r*b.i;

	return out;
	
}

doublecomplex complex_sub(doublecomplex a, doublecomplex b)
{
	doublecomplex out;
	
	out.r = a.r - b.r;
        out.i = a.i - b.i;

	return out;
	
}

doublecomplex complex_div(doublecomplex a, doublecomplex b)
{
    double ratio, den;
    double abr, abi, cr, ci;
	doublecomplex out;

    if( (abr = b.r) < 0.)
	abr = - abr;
    if( (abi = b.i) < 0.)
	abi = - abi;
    if( abr <= abi ) {
		if (abi == 0) {
		  fprintf(stderr, "complex_div.c: division by zero");
		  exit (-1);
		}	  
		
		ratio = b.r / b.i ;
		den = b.i * (1 + ratio*ratio);
		cr = (a.r*ratio + a.i) / den;
		ci = (a.i*ratio - a.r) / den;
    } 
	else {
		ratio = b.i / b.r ;
		den = b.r * (1 + ratio*ratio);
		cr = (a.r + a.i*ratio) / den;
		ci = (a.i - a.r*ratio) / den;
    }
    
	out.r = cr;
        out.i = ci;

	return out;
}

double PhaseShift(double kin,int am,int n,double rm,double V0)
{
	FILE *wave;
	double x,y,dum,h=1.0/n,r,rr,E=kin*kin,F1,V,out;
	double Co0, A0=(4.0/(360.0*h*h)+alpha*12.0/(720.0*h)),
		   Co1, A1=-(54.0/(360.0*h*h)+alpha*108.0/(720.0*h)),
		   Co2, A2=(536.0/(360.0*h*h)+alpha*528.0/(720*h)),
		   Co3, A3=(540.0/(360.0*h*h)+alpha*540.0/(720.0*h)),
		   Co4, A4=(-926.0/(360.0*h*h)+alpha*108.0/(720.0*h)),
		   Co5, 
		   Co6, A6=(536.0/(360.0*h*h)-alpha*552.0/(720.0*h)),
		   Co7, A7=(540.0/(360.0*h*h)-alpha*540.0/(720.0*h)),
		   Co8, A8=-(54.0/(360.0*h*h)-alpha*108.0/(720.0*h)),
		   Co9, A9=(4.0/(360.0*h*h)-alpha*12.0/(720.0*h));
	double *b;
	double *L;
	double *a;
		
	int *indx;
	int i,j,II,k,l,m,mm,m1,m2;

        if((wave = fopen("wave.txt","w")) == NULL)
	{
		printf("  can't open wavet file\n");
		return 0;
	}

	if( !(b = (double *) LU_MALLOC(n*sizeof(double))) ) printf("  LU_MALLOC fails for b[].");
        if( !(L = (double *) LU_MALLOC(/*4*12*/3*n*sizeof(double))) ) printf("  LU_MALLOC fails for L[].");
        if( !(a = (double *) LU_MALLOC(/*4*25*/7*n*sizeof(double))) ) printf("  LU_MALLOC fails for a[].");
   
	if ( !(indx = intMalloc(n)) ) printf("  Malloc fails for indx.");	
	
	m1=3; 
	m2=3;
	m=m1+m2+1;  // The number of columns of matrix a

									 //		Make the 2-dimentional array a as follows:
									 //		The diagonal elements are in a[0..n-1][3].  Subdiagonal elements are in a[i..n-1][0..2]
									 //     (with i>0 appropriate to the number of elements on each subdiagonal).  Superdiagonal 
									 //		elements are in a[0..i][4..6] with i<n-1 appropriate to the number of elements on each 
									 //     superdiagonal. 

	II=0;
	for (i=0;i<n;i++)
	{
        r=rm*(exp(alpha*(i+1.0)/n)-1.0)/(exp(alpha)-1.0);
		F1=(exp(alpha)-1.0)/(alpha*(rm+(exp(alpha)-1.0)*r));
		F1=F1*F1;
		Co0=F1*A0;
		Co1=F1*A1;
		Co2=F1*A2;
		Co3=F1*A3;
		Co4=F1*A4;
		Co5=-F1*980.0/(360.0*h*h);
		Co6=F1*A6;
		Co7=F1*A7;
		Co8=F1*A8;
		Co9=F1*A9;
			
		
		for (j=i-3;j<i+4;j++)  //  7 is because of 7-point approximation
		{
			if ((j<0) || (j>n-1))	a[II]= 0.0;
			else
			{	
				if (j==i-3) 
				{
					if (i<n-3)	a[II]= Co0;
					else	a[II]= 0.0;
				}
				else if (j==i-2) 
				{
					if (i<n-3)	a[II]= Co1;
					else	a[II]= 0.0;
				}
				else if (j==i-1)
				{ 
        			if(i==1)	a[II]= Co2;
					else if(i > n-4)	a[II]= nl(am,kin*r);
					else	a[II]= Co3;
				}
				else if (j==i)
				{
					V=potential0(r,V0)+ am*(am+1.0)/(/*2*M*/r*r);
					if(i==0)	a[II]= Co4 +/*2*M*/(2.0*E-V);
					else if(i > n-4)
					{
						rr=rm*(exp(alpha*(double)i/n)-1.0)/(exp(alpha)-1.0);
						a[II]= -nl(am,kin*rr);
					}
					else a[II]= Co5 +/*2*M*/(2.0*E-V);
									
				}
				else if (j==i+1) 
				{
					if(i==0)	a[II]= Co6;
					else if (i>=n-3)	a[II]= 0.0;	
					else	a[II]= Co7;
				}
				else if (j==i+2) 		
				{
					if (i<n-3)	a[II]= Co8;
					else	a[II]= 0.0;
				}
				else if (j==i+3) 
				{
					if (i<n-3)	a[II]= Co9;
					else a[II]= 0.0;
				}
				else a[II]= 0.0;
			}

			II++;
		}

					
	}								 //	   Factorization of a 
	
	mm=m-1;
	l=m1;

	for (i=0;i<m1;i++)  // Rearrange the storage a bit.
	{
		for (j=m1-i;j<=mm;j++) a[i*m+j-l]=a[i*m+j];
		l--;
		for (j=mm-l;j<=mm;j++) a[i*m+j]=0;
	}
	
	l=m1-1;
	for (k=0;k<=n-1;k++)  // For each row ...    
	{
		dum=a[k*m+0];
		i=k;
		if (l<n-1) l++;
		for (j=k+1;j<=l;j++)  // Find the pivot element.
		{
			if (fabs(a[j*m+0]) > fabs(dum))
			{
				dum=a[j*m+0];
				i=j;
			}
		}
		indx[k]=i;
		if (dum==0.0) a[k*m+0]=TINY;
		if (i!=k)  // Interchange rows.
		{
			
			for (j=0;j<=mm;j++)  SWAP(a[k*m+j],a[i*m+j]);
		}
		for (i=k+1;i<=l;i++)  // Do the elimination.
		{
			dum=a[i*m+0]/a[k*m+0];
			L[k*m1+i-k-1]=dum;
			for (j=1;j<=mm;j++) a[i*m+j-1]=a[i*m+j]-dum*a[k*m+j];
			a[i*m+mm]=0.0;
		}
	}
		

									//     Make the 1-dimentional array b
	for(i=0;i<n;i++) b[i]=0.0;

	x=rm*(exp(alpha*(1.0-2.0/n))-1.0)/(exp(alpha)-1.0);
	y=rm*(exp(alpha*(1.0-3.0/n))-1.0)/(exp(alpha)-1.0);
	b[n-3]=(-nl(am,kin*y)*jl(am,kin*x)+nl(am,kin*x)*jl(am,kin*y));  
		
	x=rm*(exp(alpha*(1.0-1.0/n))-1.0)/(exp(alpha)-1.0);
	y=rm*(exp(alpha*(1.0-2.0/n))-1.0)/(exp(alpha)-1.0);
	b[n-2]=(-nl(am,kin*y)*jl(am,kin*x)+nl(am,kin*x)*jl(am,kin*y));  
	
	x=rm;
	y=rm*(exp(alpha*(1.0-1.0/n))-1.0)/(exp(alpha)-1.0);
	b[n-1]=(-nl(am,kin*y)*jl(am,kin*x)+nl(am,kin*x)*jl(am,kin*y));

								   //     Solve the linear system problem ax = b	 							
	
	l=m1-1;
	for (k=0;k<=n-1;k++)
	{
		i=indx[k];
		if (i!=k) SWAP(b[k],b[i])
		if(l<n-1) l++;
		for (i=k+1;i<=l;i++) b[i] -= L[k*m1+i-k-1]*b[k];
	}
	l=0;
	for (i=n-1;i>=0;i--) 
	{
		dum=b[i];
		for (k=1;k<=l;k++) dum -= a[i*m+k]*b[k+i];
		b[i]=dum/a[i*m+0];
		if (l < mm) l++;
	}
	
        for (i=n-1;i>=0;i--) 
	{
		r=rm*(exp(alpha*(1.0-(double)i/n))-1.0)/(exp(alpha)-1.0);
		fprintf(wave,"  %1.7e  %1.7e  \n ",r,b[n-i-1]);
	}

        x=rm;
	y=rm*(exp(alpha*(1.0-1.0/n))-1.0)/(exp(alpha)-1.0);

	out = -(b[n-1]*jl(am,kin*y)-b[n-2]*jl(am,kin*x))/(nl(am,kin*x)*jl(am,kin*y)-nl(am,kin*y)*jl(am,kin*x)); 

	free (b);
	free (L);
	free (indx);
	free (a);

	fclose(wave);

	return out;
}

void SphBes(int n,double x,double *sj,double *sy,double *sjp,double *syp)  // Returns spherical Bessel functions jn(x), yn(x), and their derivatives jn'(x), yn'(x)for integer n.
{
	double factor,order,rj,rjp,ry,ryp;

	if (n < 0 || x <= 0.0) printf("  bad arguments in SphBes\n");
	order=n+0.5;
	bessjy(x,order,&rj,&ry,&rjp,&ryp);
	factor=RTPIO2/*    /    */*sqrt(x);
	*sj=factor*rj;
	*sy=factor*ry;
	*sjp=factor*rjp-(*sj)/(2.0*x);
	*syp=factor*ryp-(*sy)/(2.0*x);

}

void bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp) 
// Returns the Bessel functions rj=Jv and ry=Yv and their derivatives rjp=Jv' and ryp=Yv', for positive x and for xnu = v >= 0.
{

	int i,isign,l,nl;
	double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,
		   rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2;

	if(x<= 0.0 || xnu<0.0) printf("  bad arguments in bessjy\n");

	nl=(x<XMIN ? (int) (xnu+0.5) : IMAX(0,(int)(xnu-x+1.5)));
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	w=xi2/PI;
	isign=1;
	h=xnu*xi;
	if(h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for(i=1;i<= MAXIT;i++)
	{
		b += xi2;
		d = b-d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b-1.0/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=c*d;
		h=del*h;
		if (d < 0.0) isign = -isign;
		if (fabs(del-1.0) < EPS) break;
	}

	if (i>MAXIT) printf("  x too large in bessjy; try asymptotic expansion\n");
	rjl=isign*FPMIN;
	rjpl=h*rjl;
	rjl1=rjl;
	rjp1=rjpl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--)
	{
		rjtemp=fact*rjl+rjpl;
		fact -= xi;
		rjpl=fact*rjtemp-rjl;
		rjl=rjtemp;
	}

	if (rjl == 0.0) rjl=EPS;
	f=rjpl/rjl;
	if (x < XMIN)
	{
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);  // Chebyshev evaluation of Gamma functions 1 and 2
		ff=2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);   //  f0.  
		e=exp(e);
		p=e/(gampl*PI);      //  p0.
		q=1.0/(e*PI*gammi);  //  q0.
		pimu2=0.5*pimu;
		fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
		r=PI*pimu2*fact3*fact3;
		c=1.0;
		d = -x2*x2;
		sum=ff+r*q;
		sum1=p;
		for (i=1;i<=MAXIT;i++)
		{
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*(ff+r*q);
			sum += del;
			del1=c*p-i*del;
			sum1 += del1;
			if (fabs(del) < (1.0+fabs(sum))*EPS) break;
		}

		if (i>MAXIT) printf("  bessy series failed to converge/n");
		rymu = -sum;
		ry1 = -sum1*xi2;
		rymup=xmu*xi*rymu-ry1;
		rjmu=w/(rymup-f*rymu);
	}
	else
	{
		a=0.25-xmu2;
		p = -0.5*xi;
		q=1.0;
		br=2.0*x;
		bi=2.0;
		fact=a*xi/(p*p+q*q);
		cr=br+q*fact,
		ci=bi+p*fact;
		den=br*br+bi*bi;
		dr=br/den;
		di = -bi/den;
		dlr=cr*dr-ci*di;
		dli=cr*di+ci*dr;
		temp=p*dlr-q*dli;
		q=p*dli+q*dlr;
		p=temp;
		for (i=2;i<=MAXIT;i++)
		{
			a += 2*(i-1);
			bi += 2.0;
			dr=a*dr+br;
			di=a*di+bi;
			if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
			fact=a/(cr*cr+ci*ci);
			cr=br+cr*fact;
			ci=bi-ci*fact;
			if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
			den=dr*dr+di*di;
			dr /= den;
			di /= -den;
			dlr=cr*dr-ci*di;
			dli=cr*di+ci*dr;
			temp=p*dlr-q*dli;
			q=p*dli+q*dlr;
			p=temp;
			if (fabs(dlr-1.0)+fabs(dli) < EPS) break;
		}
		if (i > MAXIT) printf("  cf2 failed in bessjy\n");
		gam=(p-f)/q;
		rjmu=sqrt(w/((p-f)*gam+q));
		rjmu=SIGN(rjmu,rjl);
		rymu=rjmu*gam;
		rymup=rymu*(p+q/gam);
		ry1=xmu*xi*rymu-rymup;
	}
	fact=rjmu/rjl;
	*rj=rjl1*fact;
	*rjp=rjp1*fact;
	for (i=1;i<=nl;i++)
	{
		rytemp=(xmu+i)*xi2*ry1-rymu;
		rymu=ry1;
		ry1=rytemp;
	}
	*ry=rymu;
	*ryp=xnu*xi*rymu-ry1;

}

double chebev(double a, double b, double c[], int m, double x)  //  Chebyshev evaluation: All arguments are input.
{
	double d=0.0, dd=0.0, sv, y, y2;
	int j;

	if ((x-a)*(x-b) > 0.0) printf("  x not in range in routine chebev\n");
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--)
	{
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)  //  Evaluates Gamma functions 1 and 2 by Chebyshev expansion
                                                                                 //  for |x|<= 1/2. 
{
	double xx;
	static double c1[]={  -1.142022680371168e0,6.5165112670737e-3,3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,3.67795e-11,-1.356e-13};
	static double c2[]={   1.843740587300905e0,-7.68528408447867e-2,1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,2.423096e-10,-1.702e-13,-1.49e-15};
	xx=8.0*x*x-1.0;
	*gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
	*gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);
	
}
	
double jl(int n,double x)  // Returns the spherical Bessel function Jn for integer n.
{
	double sj,sy,sjp,syp;

	SphBes(n,x,&sj,&sy,&sjp,&syp);
	return sj;
	
}

double nl(int n,double x)  // Returns the spherical Bessel function Yn for integer n.
{
	double sj,sy,sjp,syp;

	SphBes(n,x,&sj,&sy,&sjp,&syp);
	return sy;

}

double plgndr(int l,int mm,double x)		// Computes the associated Legendre polynomial P{l,m}(x).
{										// Here m and l are integers satisfying 0<= m <= l, while x 
	int i,ll,m;							// lies in the range -1<= x <= 1.
	double fact,pll,pmm,pmmp1,somx2,output;

        m = abs(mm);

	if((m<0) || (m>l) || (fabs(x)>1.0))
		printf(" Eror: Bad arrguments in routine plgndr!\n");

	pmm=1.0;		// Compute P{m,m}.
	if(m>0)
	{
		somx2 = sqrt((1.0-x)*(1.0+x));
		fact = 1.0;
		for(i=1;i<=m;i++)
		{
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}

	if(l== m)
		output = pmm;
	else
	{
		pmmp1 = x*(2*m+1)*pmm;
		if(l==(m+1))
			output = pmmp1;
		else		// Compute P{l,m}, l > m+1.
		{
			for(ll=m+2;ll<=l;ll++)
			{
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			output = pll;
		}
	}

	if (mm < 0)
	{
		for(i=1; i<m;i++)
		    output *= (-1)/((l+i)*(l-i));
	    output *= (-1)/(l*(l+m));
	}
	
	return output;
}

double pl(int l,int m,double x){
          double out;
		if(m<0){out=sign(m)*(sqrt((2.0*l+1.0)*factorial(l-abs(m))/(2.0*factorial(l+abs(m))))*plgndr(l,abs(m),x));}
		else out=sqrt((2.0*l+1.0)*factorial(l-m)/(2.0*factorial(l+m)))*plgndr(l,m,x);
	return out;
}
void gauleg(double x1, double x2, double x[], double w[], int n) 
// Given the lower and upper limits of integration x1 and x2, and given n, this routin 
// returns arrays x[0..n-1] and w[0..n-1] of length n, containing the abscissas and weights 
// of the Gauss-Legendre n-point quadrature formula.
{
	double xm,xl,p1,p2,p3,pp,z,z1;
	int i,j,m;
	
	for(i=0;i<n;i++)	x[i]=w[i]=0.0;

	m = (n+1)/2;
	xm = 0.5*(x2 + x1);
	xl = 0.5*(x2 - x1);
	for(i=1;i<=m;i++)
	{
		z=cos(PI*(i-0.25)/(n+0.5));
		do		// Starting with the above approximation to the ith root, we enter the 
		{		// main loop of the refinement by Newton's method.
			p1=1.0;
			p2=0.0;
			for(j=1;j<=n;j++)		// loop up the recurrence relation to get the 
			{						// Legendre polynomial evaluated at z.
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}		// p1 is now the desired Legendre polynomial.  pp is its derivative 
					// which is computed by a standard relation involving p2, the 
					// polynomial of one lower order.
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while(fabs(z-z1)> EPS2);
		x[i-1]=xm-xl*z;		// The abscisses
		x[n-i]=xm+xl*z;
		w[i-1]=2.0*xl/((1-z*z)*pp*pp);		// The weights
		w[n-i]=w[i-1];
	}
}

doublecomplex *doublecomplexMalloc(int n)
{
    doublecomplex *buf;
    doublecomplex zero = {0.0, 0.0};
	register int i;

	buf = (doublecomplex *) LU_MALLOC(n * sizeof(doublecomplex)); 
        if ( !buf ) printf("  LU_MALLOC failed for buf in doublecomplexMalloc()\n");
   
	for (i = 0; i < n; ++i) buf[i] = zero;

        return (buf);
}

int *intMalloc(int n)
{
    int *buf;
	register int i;

    buf = (int *) LU_MALLOC(n * sizeof(int));
    if ( !buf )		printf("LU_MALLOC fails for buf in intMalloc()");
     
	for (i = 0; i < n; ++i) buf[i] = 0;

        return (buf);
}


void *LU_MALLOC(size_t size)
{
    void *buf;
    buf = (void *) malloc(size);
    return (buf);
}

double **D_Matrix_alloc(int m,int n)
{
	double **matrix;
	int row,col;

	matrix = (double **)malloc( m * sizeof(double *) );

    for ( row = 0; row < m; row++ )
		matrix[row] = (double *)malloc( n * sizeof(double) );
        
    for ( row = 0; row < m; row++ )
        for ( col = 0; col < n; col++ )
               matrix[row][col] = 0.0;

	return (matrix);
}

doublecomplex **DCom_Matrix_alloc(int m,int n)
{
	doublecomplex **matrix;
	int row,col; 
	doublecomplex zero = {0.0, 0.0};

    matrix = (doublecomplex **)LU_MALLOC(m * sizeof(doublecomplex*) );
		
	for ( row = 0; row < m; row++ )
		matrix[row] = doublecomplexMalloc(n);
      
    for ( row = 0; row < m; row++ )
        for ( col = 0; col < n; col++ )
               matrix[row][col] = zero;

	return (matrix);
}

void initialize(double **A,int m,int n)
{
	int row,col;

	for (row=0;row<m;row++)
		for (col=0;col<n;col++)
			A[row][col]=0.0;
}

void c_initialize(doublecomplex **A,int m,int n)
{
	int row,col;
	doublecomplex zero={0.0,0.0};

	for (row=0;row<m;row++)
		for (col=0;col<n;col++)
			A[row][col]=zero;
}

double delta(int i,int j)
{
	if (i==j)	return 1.0;
	else return 0.0;
}
int IMAX(int i,int j)
{
	if (i>j) return(i);
	else return(j);
}

double SIGN(double x,double y)
{
	if (y>0 ) return x;
	else return (-1)*x;
}

double factorial(int n)
{
  double f;
  if(n==0)
     return 1;
  else
      {
       if(n>0)
       f=n*factorial(n-1);
      }
      return f;
}

int factoril_gen(int n, int m)
{
	int j,out=1;
	if(n>m) for (j=n;j>m;j--)	out = out*j;
	
	return out;
}
