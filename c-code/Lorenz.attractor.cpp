#include <math.h>
#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

double C0=0.0,P=1.0,sigma=-1.1,u0;

typedef std::vector< double > state_type;
           
void odemodel (const state_type &x , state_type &dxdt ,double t )
{
    dxdt[0] =  x[1];
    dxdt[1] = -1*x[1]*cos(x[0])/x[3]+cos(x[0])*sin(x[0])/(x[3]*x[3])+x[2]*sin(x[0])/x[3]+P*x[3]*cos(x[0])/2.;
    dxdt[2] = (x[1]-C0)*(x[1]-C0)/2.-sin(x[0])*sin(x[0])/(2*x[3]*x[3])+P*x[3]*sin(x[0])+sigma;
    dxdt[3] = cos(x[0]);
}

void write_output(const state_type &x , double t )
{
    cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << x[3] << endl;
}

int main(int argc, char **argv)

{	
	FILE *fp;
	char fiche[100];
	sprintf(fiche,"results.dat");
	fp=fopen(fiche,"w");
	fclose(fp);	
	
for (u0=0.;u0<0.1;u0=u0+0.001){	

    state_type yyy = {0.0,u0,0.0,0.0};
    integrate(odemodel, yyy ,0.0 ,25.0 ,0.001 , write_output);
    
  	fp=fopen(fiche,"a");
  	fprintf(fp,"%f %f %f %f %f %f \n",u0,x[0],x[1],x[2],x[3],sigma);
  	fclose(fp);
}}
