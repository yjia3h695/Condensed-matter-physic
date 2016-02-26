#include <iostream>
#include <boost/array.hpp>
#include <fstream>
#include <math.h>

#include <boost/numeric/odeint.hpp>

const double PI = 3.1415926535897;
using namespace std;
using namespace boost::numeric::odeint;

double error = 1e-2;
double C0 = 0;

typedef std::vector< double > state_type;

struct lorenz_class
{ 
	double m_P;
	lorenz_class(double P) : m_P(P) {}
 	
	void operator()( state_type &x , state_type &dxdt , double t)
	{
	double Sigma = -1.1*pow(m_P,2.0/3.0);

  	dxdt[0]=x[1]; 
 	dxdt[1]=-x[1]*cos(x[0])/x[3]+cos(x[0])*sin(x[0])/(x[3]*x[3])+x[2]*sin(x[0])/x[3]+m_P*x[3]*cos(x[0])/2.;
 	dxdt[2]=(x[1]-C0)*(x[1]-C0)/2.-sin(x[0])*sin(x[0])/(2.*x[3]*x[3])+m_P*x[3]*sin(x[0])+Sigma;
 	dxdt[3]=cos(x[0]);
	}
};

struct streaming_observer
{
     std::ostream &m_out;
     double m_u0; double m_P;
     streaming_observer( std::ostream &out, double u0, double P) : m_out(out),m_u0(u0),m_P(P) {}
	 double Sigma = -1.1*pow(m_P,2.0/3.0);	
     void operator()( const state_type &x , double t ) const
     {
     	if (x[0]<PI+error && x[0]>PI-error)
     	{
        	m_out << m_u0 <<"\t" << m_P <<"\t" << Sigma <<"\t" << t;        
        	for( size_t i=0 ; i < x.size() ; ++i )
            m_out << "\t" << x[i];
       	 	m_out << "\n";  
		} 
     }
};

string intToString(double t)
{
  std::string ch;
  ostringstream outs; 
  outs << t;   // Convert value into a string.
  ch = outs.str();   
     
  return ch;
} 

int main(int argc, char **argv)
{
    state_type x(4);
    std::remove("output50-100.txt"); // delete file
    for (double P=51;P<=100;P+=1){
//		std::string name = "output" + intToString(P) + ".txt";  
//		ofstream file (name.c_str(),std::ofstream::out | std::ofstream::app);
		double scale = pow(P,-1.0/3.0);
    	for (double u0=-2.0/scale;u0<=2.0/scale;u0+=0.01)	{
    		x[0] = x[2]=0.0,x[3]=1e-4;x[1]=u0;	
			ofstream file ("output.txt",std::ofstream::out | std::ofstream::app);	
	
//			typedef controlled_runge_kutta< runge_kutta_cash_karp54< state_type > > stepper_type;
//			integrate_adaptive( stepper_type() , lorenz_class() , x , 0.0 , 30.0 , 0.00001 , streaming_observer( file, u0) );

 			typedef runge_kutta_dopri5< state_type > dopri5_type;
    		typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
    		typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
    		dense_output_dopri5_type dopri5 = make_dense_output( 1E-12 , 1E-12 , dopri5_type() );
    		integrate_adaptive( dopri5 , lorenz_class(P) , x , 0.0 , 30.0 , 0.01 , streaming_observer( file, u0, P) );
		}
	}
}
