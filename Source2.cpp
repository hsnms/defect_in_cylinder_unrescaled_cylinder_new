//Source.cpp helps to find the transition lines equilibrium nematic defect structures in an unrescale cylindrical nematic bridge given specific parameters including the Frank constant ratio and aspect ratio. Case (1): keep the thin layer near the lateral surface.

//We consider homeotropic boundary conditions and proper inner boundary conditions (after we apply a cut-off of the defect core), and solve the Euler-Lagrange equation. Because of the cylindrical symmetry, we only consider the rectangular diametrical plane.

//The finite difference method (particularly the over-relaxation method) is used. The rectangular region is not rescaled.

//Input: n: number of lattices horizontally and vertically; R0: aspect ratio; K:Frank constant ratio; b: the number of the lattice spaces the length of the defect core occupies; m:types of defect structures (radial:1.0, hyperbolic:-1.0).

//Output: 1,"Bridge.txt": different energies when the defect core is at different locations, error, all the input values; 2, "BridgeVectorField.txt": the vector field.

//How to compile: 1, make sure the in the "makefile", we have the relevant "Source2"; 2, M-X compile (option + X + compile); 3, Compile command: make -k output; 4, Save file, y; 5, M-X shell; 6, bash-3.2$ ./output





#include <stdio.h>
#include <errno.h>
#include <math.h>
#define Up 1.570796326794897
#define Down1 1.570796326794897
#define Down2 1.570796326794897
#define Right 1.570796326794897
#define Left 1.570796326794897
#define Pi 3.141592653589793
#define S 0.00000000000001
#include "Header1new.h"
#include "Header2new.h"
/*#include "Header3new.h"*/
/*#include "Header4new.h"*/
#include "Header5new.h"
#include "Header6new.h"
#include "Header7new.h"
int main()
{
    int i,j,k,n,ncycle=5,ii,ii2,ii3,ii4,sos,kk,radius,min,n1;
    double **u,*x,*y,*z,*z2,*q,**v,**f;
	double result,zero,rev2,rev3,v1,v2,v4,K,t,sum,result22,sum2,fff,dfff;
	double m,point,point2,H,h,R0,R1,R2,K1,K2,RR,KK; int a,b; 
	FILE *fp,*fp2,*fp3;

	printf ("Number of lattices vertically, n:");
	scanf("%d", &n);
	//n=17;

        printf ("Initial aspect ratio, R1:");
	scanf("%lf", &R1);
	printf ("Final aspect ratio, R2:");
	scanf("%lf", &R2);
        printf ("Aspect ratio increment, RR:");
	scanf("%lf", &RR);
	
	H=1.0;

	printf ("Initial Frank constant ratio, K1:");
	scanf("%lf", &K1);
	printf ("Final Frank constant ratio, K2:");
	scanf("%lf", &K2);
	printf ("Frank constant ratio increment, KK:");
	scanf("%lf", &KK);
	//K=0.65;

	printf ("How many of the lattice spaces does the length of the defect core occupy? b:");
	scanf("%d", &b);
	//b=1;/*radius of the point or small radius of the ring*/
	
	  kk=1+2*b;/*location of the smallest ring*/
	  // printf ("Types of defect structures (radial:1.0, hyperbolic:-1.0), m:");
	  //	scanf("%lf", &m);
	  m=1.0;
	
	fp = fopen ("Transition.txt","w");
	  if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	  //printf("Number of lattices vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\nDefect type m=%f\n",n,b,m);

	  printf("Number of lattices vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\n",n,b);//something new
	  
//fprintf(fp,"Number of lattices vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\nDefect type m=%f\n",n,b,m);

	  fprintf(fp,"Number of lattices vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\n",n,b);//something new


   for (K=K1;K<=K2;K=K+KK) //different Frank constant ratios
   for(R0=R1;R0<=R2;R0=R0+RR)//different aspect ratios
	     {	  

	       //	m=-1.0;
	  
    u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	v=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	q=dvectornew(1,(int)(n*H-H+1));
	
	x=dvectornew(1,(int)(n*R0-R0+1));
	y=dvectornew(1,(int)(n*H-H+1));
	z=dvectornew(1,(int)(n*R0-R0+1));
	z2=dvectornew(1,(int)(n*H-H+1));
	for (i=1;i<=(int)(n*R0-R0+1);i++) /*set the length*/
	{
        
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)/*set the width*/
	{
        
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}

	for (j=1;j<=(int)(n*H-H+1);j++)/*the side boundary*/
		{
			
			q[j]=R0;
		
		}
	
	for(i=1;i<=(int)(n*R0-R0+1);i++)/*initial condition and boundary condition*/
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=b;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+b;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			   
		  u[(int)(n*R0-R0+1)][j]=Right/*+1.0-1.0*(j-1)/((int)(n*H-H))*/;   
		
			u[1][j]=Left-0.5*m*Pi;
	    }


	if (m==-1.0) /*hyperbolic type*/
	  { /*start of the first major part*/
	    for (j=2;j<=1+b;j++) /*inner boundary condition*/
	{
	  //	u[1+b][j]=atan((y[j]-y[1])/((x[1+b]-x[1])))+0.5*Pi-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*atan((y[j]-y[1])/((x[1+b]-x[1])))));
		  u[1+b][j]=atan((y[j]-y[1])/((x[1+b]-x[1])))+0.5*Pi-((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*atan((y[j]-y[1])/((x[1+b]-x[1])))));	
	}
	for (i=2;i<=1+b;i++)
	{
	  //	u[i][1+b]=atan((y[1+b]-y[1])/((x[i]-x[1])))+0.5*Pi-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*atan((y[1+b]-y[1])/((x[i]-x[1])))));
		   u[i][1+b]=atan((y[1+b]-y[1])/((x[i]-x[1])))+0.5*Pi-((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*atan((y[1+b]-y[1])/((x[i]-x[1])))));	
	}
	
		for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
		u[i][(int)(n*H-H+1)-b]=Pi-atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])));
		
	}
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
    {
		u[(int)(n*R0-R0+1)-b][j]=Pi-atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])));
		
		}
	
	mgfas5new(u,n,ncycle,m,b,q,x,y, R0,K, H); /*solve the PDE*/
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc4new(f,u,n,q,b,R0,K,H); /*get the free energy density*/
	
n1=1;
 z[1]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, n1,n1+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0); /*get the free energy*/
	zero=Stat3new(u,n,b,R0,K, H);
	result=z[1];
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	radius=1;
	
	printf("z[1]=%f\t\nzero=%f\n",z[1],zero);
	//	fprintf(fp,"z[1]=%f\t\n{0, %f},\nzero=%f\n",z[1],z[1],zero);
	//	fprintf(fp,"{0, %f},\t",z[1]);


  //Ring defect with smallest radius
		a=1+2*b;kk=a;
	    //initial and outer boundary conditions
	for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }	
		for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
		{
			for(k=2;k<=(int)(n*H-H+1)-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			u[(int)(n*R0-R0+1)][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }
	//inner boundary
	for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	    u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
	    u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
	}
	//	u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
		u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
     
      //	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
           u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
           u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[j]-y[1])/((x[a]-x[a-b]))))));
	}
	for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
		u[i][(int)(n*H-H+1)-b]=Pi-atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])));
	}
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
    {
     
      u[(int)(n*R0-R0+1)-b][j]=Pi-atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])));
	}
	for (j=1;j<=(int)(n*H-H+1);j++)
		q[j]=R0;
       mgfas4new(u,n,ncycle,m,a,b,q,x,y, R0,K, H);
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
	nrfunc3new(f,u,n,q,kk,b,R0,K,H);
		
n1=1;
	z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);
	zero=Stat2new(u,n,a,b,R0,K,H);


	if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}

	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	//	fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);

		// ring defect with different radii

	for(kk=1+2*b;kk<=(int)(n*R0-R0+1)-b-1;kk++)
	{
	 a=kk;
	

	u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	q=dvectornew(1,(int)(n*H-H+1));
	
	
	for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
        
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
        
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}



	//initial and outer boundary conditions
	
	for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }
for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
		{
			for(k=2;k<=(int)(n*H-H+1)-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			u[(int)(n*R0-R0+1)][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }
	//inner boundary conditions
	for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	   u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
	     u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
	}
	//	u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
		u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      // u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
             u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
        u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[j]-y[1])/((x[a]-x[a-b]))))));
	}
	for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
	 
		u[i][(int)(n*H-H+1)-b]=Pi-atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])));
	}
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
    {
     
		u[(int)(n*R0-R0+1)-b][j]=Pi-atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])));
	}

	for (j=1;j<=(int)(n*H-H+1);j++)
		q[j]=R0;
       mgfas4new(u,n,ncycle,m,a,b,q,x,y, R0,K, H);
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
nrfunc3new(f,u,n,q,kk,b,R0,K,H);
		
	n1=1;
	z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);
	zero=Stat2new(u,n,a,b,R0,K,H);

	if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}
	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	
	//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	//	fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);
	}
	  } /*end of the first major*/

	//	printf("Defect type m=%f, Frank constant ratio K=%f, aspect ratio R0=%f, ring radius=%d, total energy=%f\n",m, K,R0,radius,result);	
	//	fprintf(fp,"Defect type m=%f, Frank constant ratio K=%f, aspect ratio R0=%f, ring radius=%d, total energy=%f\n",m, K,R0,radius,result);

<<<<<<< HEAD
			m=1.0;
=======
		m=1.0;
>>>>>>> b541eeffca17024a3f1e3fcda23d313d5931029c
	if (m==1.0)
	    {  /*start of the first major*/
               //point defect
	  for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
	       
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
        
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}

for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }

 	for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=b;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+b;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;
	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			   
		  u[(int)(n*R0-R0+1)][j]=Right/*+1.0-1.0*(j-1)/((int)(n*H-H))*/;   
		
			u[1][j]=Left-0.5*m*Pi;
	    }
	//inner boundary
		for (j=2;j<=1+b;j++)
	{
		u[1+b][j]=-atan((y[j]-y[1])/((x[1+b]-x[1])))+0.5*Pi;
		
	}
	for (i=1;i<=1+b;i++)
	{
		u[i][1+b]=-atan((y[1+b]-y[1])/((x[i]-x[1])))+0.5*Pi;
		
	}
	
	for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)//important
	{
	  //	u[i][(int)(n*H-H+1)-b]=atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))-Pi)));
		u[i][(int)(n*H-H+1)-b]=atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))-Pi)));	
		
	}
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
    {
      //	u[(int)(n*R0-R0+1)-b][j]=atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-Pi)));
      	u[(int)(n*R0-R0+1)-b][j]=atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-Pi)));
		
	}
	for (j=1;j<=(int)(n*H-H+1);j++)
		q[j]=R0;
        mgfas5new(u,n,ncycle,m,b,q,x,y, R0,K, H);
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
	nrfunc4new(f,u,n,q,b,R0,K,H);
		
n1=1;
	z[1]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, n1,n1+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);
	zero=Stat3new(u,n,b,R0,K, H);

	result=z[1];
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	radius=1;
	
	printf("z[1]=%f\t\nzero=%f\n",z[1],zero);
	//	fprintf(fp,"z[1]=%f\t\n{0, %f},\nzero=%f\n",z[1],z[1],zero);
	//	fprintf(fp,"{0, %f},\t",z[1]);
	


	a=1+2*b;//ring defect with smallest radius
		kk=a;
			for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }
		for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
		{
			for(k=2;k<=(int)(n*H-H+1)-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			u[(int)(n*R0-R0+1)][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }
	//inner boundary
for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	    u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));
	  	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));
	}
	//	u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(0.5*Pi));
		u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
       u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
      	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
	}
	for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
	  //	u[i][(int)(n*H-H+1)-b]=atan(H*(y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/(R0*(x[i]-x[(int)(n*R0-R0+1)])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/(R0*(x[i]-x[(int)(n*R0-R0+1)])))-Pi)));
	  	u[i][(int)(n*H-H+1)-b]=atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))-Pi)));
		
	}
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
    {
      //	u[(int)(n*R0-R0+1)-b][j]=atan(H*(y[j]-y[(int)(n*H-H+1)])/(R0*(x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[j]-y[(int)(n*H-H+1)])/(R0*(x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-Pi)));
      	u[(int)(n*R0-R0+1)-b][j]=atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-Pi)));
		
	}
	for (j=1;j<=(int)(n*H-H+1);j++)
		q[j]=R0;
        mgfas4new(u,n,ncycle,m,a,b,q,x,y, R0,K, H);
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
	nrfunc3new(f,u,n,q,kk,b,R0,K,H);
		
n1=1;
  z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);
	zero=Stat2new(u,n,a,b,R0,K,H);

	if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}

	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	//	fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);

		//ring defect with different radii

	for(kk=1+2*b;kk<=(int)(n*R0-R0+1)-b-1;kk++)
	{
	 a=kk;
	

	u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	q=dvectornew(1,(int)(n*H-H+1));
	
	
	for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}


	for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }

for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
		{
			for(k=2;k<=(int)(n*H-H+1)-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			u[(int)(n*R0-R0+1)][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }
	//inner boundary
	for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	   u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));
	  	  	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));

	}
	//	u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(0.5*Pi));
		u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
       u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
      	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
	}
	for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
	  //	u[i][(int)(n*H-H+1)-b]=atan(H*(y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/(R0*(x[i]-x[(int)(n*R0-R0+1)])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))-Pi)));
		u[i][(int)(n*H-H+1)-b]=atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])))-Pi)));	
	}
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
    {
      //	u[(int)(n*R0-R0+1)-b][j]=atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-Pi)));
	  	u[(int)(n*R0-R0+1)-b][j]=atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])))-Pi)));	
	}
	for (j=1;j<=(int)(n*H-H+1);j++)
		q[j]=R0;
        mgfas4new(u,n,ncycle,m,a,b,q,x,y, R0,K, H);
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
	nrfunc3new(f,u,n,q,kk,b,R0,K, H);
		
	n1=1;

	z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);
	zero=Stat2new(u,n,a,b,R0,K,H);

	if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}
	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	
	//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	//	fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);
	}
	} /*end of the first major part*/
	


			
	//	printf("result%f\t\nradius=%d\n",result,radius);	
	//	fprintf(fp,"m=%f\tR0=%f\nresult=%f\tradius=%d\n",m,R0,result,radius);
	//	printf("Total energy=%f\t\nRing radius=%d\n",result,radius);

	//	printf("Defect type m=%f, Frank constant ratio K=%f, aspect ratio R0=%f, ring radius=%d, total energy=%f\n",m,K,R0,radius,result);
	//	fprintf(fp,"\nNumber of lattices vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\nDefect type m=%f\nAspect ratio R0=%f\nFrank constant ratio K=%f\nTotal energy=%f\nRing radius=%d\n",n,b,m,R0,K,result,radius);

		printf("Defect type m=%f, Frank constant ratio K=%f, aspect ratio R0=%f, ring radius=%d, total energy=%f\n",m,K,R0,radius,result);	
	fprintf(fp,"Defect type m=%f, Frank constant ratio K=%f, aspect ratio R0=%f, ring radius=%d, total energy=%f\n",m,K,R0,radius,result);	
	     }

		fclose (fp);
                   
		
		

		/*	fprintf(fp2,"{");
		for(i=1;i<=(int)(n*R0-R0+1)-1;i++)
		{
			fprintf(fp2,"{");
			for(k=1;k<=(int)(n*H-H+1)-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(v[i][k]),cos(v[i][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(v[i][(int)(n*H-H+1)]),cos(v[i][(int)(n*H-H+1)]));
			fprintf(fp2,"},");
		}
			fprintf(fp2,"{");
			for(k=1;k<=(int)(n*H-H+1)-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(v[(int)(n*R0-R0+1)][k]),cos(v[(int)(n*R0-R0+1)][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(v[(int)(n*R0-R0+1)][(int)(n*H-H+1)]),cos(v[(int)(n*R0-R0+1)][(int)(n*H-H+1)]));
			fprintf(fp2,"}");
		fprintf(fp2,"}");
		fclose (fp2);*/
		
		printf ("File created okay\n");

	
	free_dmatrix2new(u,1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	free_dmatrix2new(v,1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	free_dvectornew(q,1,(int)(n*H-H+1));
	
	free_dvectornew(x,1,(int)(n*R0-R0+1));
	free_dvectornew(y,1,(int)(n*H-H+1));
	free_dvectornew(z,1,(int)(n*R0-R0+1));
	free_dvectornew(z2,1,(int)(n*H-H+1));
    return 0;
}
