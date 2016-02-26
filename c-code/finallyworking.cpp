#include <valarray>
#include <iostream>
#include <math.h>
const double PI = 3.1415926535897;

int N = 4;
double P=1;
double C0=0;
double scale=pow(P,-1.0/3.0);
double Sigma = -1.1*pow(P,2.0/3.0);

typedef std::valarray<double> Vector;

/*-------------------------------------------------------------------------*/
/*-----------------------start of main program-----------------------------*/
/*-------------------------------------------------------------------------*/
int main()
{
	Vector f( const double t, const Vector M );
	
	FILE *fp;
	char fiche[100];
	sprintf(fiche,"day2.dat");
	fp=fopen(fiche,"w");
	fclose(fp);
	fp=fopen(fiche,"a");
  for (double u0=-2/scale;u0<=2/scale;u0=u0+0.001)
	{  	
  	Vector X(N);
  	X[0]=0;X[1]=u0;X[2]=0;X[3]=0.0001;
 	
  	double Psi,U,Gamma,Xaxis,time,t = 0;
  	double dt = 0.001;
  	
	while (X[0]<PI)
		{
			Psi=X[0];U=X[1];Gamma=X[2];Xaxis=X[3];time=t;
			Vector k1 = dt*f( t, X );								/*--------RK4-----------*/
   			Vector k2 = dt*f( t + dt / 2.0,  X + dt / 2.0 * k1 );	/*--------RK4-----------*/
    		Vector k3 = dt*f( t + dt / 2.0, X + dt / 2.0 * k2 );	/*--------RK4-----------*/
    		Vector k4 = dt*f( t + dt, X + dt * k3 );				/*--------RK4-----------*/
 	 		X = X + ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6.0 ;		/*--------RK4-----------*/
    		t = t + dt;												/*--------RK4-----------*/
    		
//			X=X+ dt*f( t, X ); 										/*-----Euler Method-----*/
//			t = t + dt;		   										/*-----Euler Method-----*/
    	} 
  	printf("\nu0=%f Psi=%f U=%f Gamma=%f X=%f t=%f Sigma=%f",u0,Psi,U,Gamma,Xaxis,time,Sigma);
  	
  	fprintf(fp,"%f %f %f %f %f %f %f \n",u0,Psi,U,Gamma,Xaxis,time,Sigma);
  	
	}
	fclose(fp);
}
/*-------------------------------------------------------------------------*/
/*-----------------------end of main program-------------------------------*/
/*-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/
/*------------------------------ODES---------------------------------------*/
/*-------------------------------------------------------------------------*/
Vector f( const double t, const Vector M )
{
  Vector output(N);
  output[0]=M[1]; 
  output[1]=-M[1]*cos(M[0])/M[3]+cos(M[0])*sin(M[0])/(M[3]*M[3])+M[2]*sin(M[0])/M[3]+P*M[3]*cos(M[0])/2.;
  output[2]=(M[1]-C0)*(M[1]-C0)/2.-sin(M[0])*sin(M[0])/(2.*M[3]*M[3])+P*M[3]*sin(M[0])+Sigma;
  output[3]=cos(M[0]);
  return output;
}

/*-------------------------------------------------------------------------*/
/*------------------------------End of ODES--------------------------------*/
/*-------------------------------------------------------------------------*/
