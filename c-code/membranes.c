#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*#include <sys/times.h>
*/
#include "ran1.c"
#include "nrutil.c"
#define PI 3.141592654
double xinew,unew,ganmanew,xnew,deltas,sigmap=-1.1;
main()
{
    double xi0,xi,xiold,x0,x,ganma0,ganma,u0,u,s,error=0.1;
    int contador;
    void RK4();
    FILE *fp;
    char fiche[100];
    sprintf(fiche,"results.dat");
    fp=fopen(fiche,"w");
    fclose(fp);
    xi0=0;x0=0.0001;ganma0=0;deltas=0.001;
    for (sigmap=-1.1;
    sigmap>=-1.2+0.001;
    sigmap=sigmap-0.1)
    {
        for (u0=2;u0>=-2;u0=u0-0.001)
        {
            xiold=xi0;xi=xi0;u=u0;ganma=ganma0;x=x0;s=0;contador=0;
            //  printf("nu0=%f xi=%f x=%f u=%f s=%f",u0,xi,x,u,s);
            do
            {
                s=s+deltas;
                RK4(xi,u,ganma,x);
                xiold=xi;
                xi=xinew;
                u=unew;
                ganma=ganmanew;
                x=xnew;
            }
            while (xinew<PI);
            printf("nu0=%f xi=%f u=%f ganma=%f x=%f s=%f sigmap=%f \n",u0,xinew,unew,ganmanew,xnew,s,sigmap);
            fp=fopen(fiche,"a");
            fprintf(fp,"%f %f %f %f %f %f %f \n",u0,xinew,unew,ganmanew,xnew,s,sigmap);
            fclose(fp);
        }
    }

            }
            /* .............................................................. */
            /* ............... Fin Programa y comienzo de Funciones ......... */
            /* .............................................................. */
            /* .......................... RK4 ................................*/
            void RK4(xi,u,ganma,x)
            double xi,u,ganma,x;
            {
                double xik,uk,ganmak,xk;
                double Kxi1,Ku1,Kganma1,Kx1,Kxi2,Ku2,Kganma2,Kx2,Kxi3,Ku3,Kganma3,Kx3,Kxi4,Ku4,Kganma4,Kx4;
                
                xik=xi;uk=u;ganmak=ganma;xk=x;
                
                Kxi1=deltas*(uk);
                Ku1=deltas*(-uk*cos(xik)/xk+cos(xik)*sin(xik)/(xk*xk)+ganmak*sin(xik)/xk+xk*cos(xik)/2);
                Kganma1=deltas*(uk*uk/2.-sin(xik)*sin(xik)/(2.*xk*xk)+xk*sin(xik)+sigmap);
                Kx1=deltas*(cos(xik));
                
                xik=xi+Kxi1/2.;uk=u+Ku1/2.;ganmak=ganma+Kganma1/2.;xk=x+Kx1/2.;
                
                Kxi2=deltas*(uk);
                Ku2=deltas*(-uk*cos(xik)/xk+cos(xik)*sin(xik)/(xk*xk)+ganmak*sin(xik)/xk+xk*cos(xik)/2);
                Kganma2=deltas*(uk*uk/2.-sin(xik)*sin(xik)/(2.*xk*xk)+xk*sin(xik)+sigmap);
                Kx2=deltas*(cos(xik));
                
                xik=xi+Kxi2/2.;uk=u+Ku2/2.;ganmak=ganma+Kganma2/2.;xk=x+Kx2/2.;
				
				Kxi3=deltas*(uk);                
                Ku3=deltas*(-uk*cos(xik)/xk+cos(xik)*sin(xik)/(xk*xk)+ganmak*sin(xik)/xk+xk*cos(xik)/2);
                Kganma3=deltas*(uk*uk/2.-sin(xik)*sin(xik)/(2.*xk*xk)+xk*sin(xik)+sigmap);
                Kx3=deltas*(cos(xik));
                
                xik=xi+Kxi3;uk=u+Ku3;ganmak=ganma+Kganma3;xk=x+Kx3;
                
                Kxi4=deltas*(uk);
                Ku4=deltas*(-uk*cos(xik)/xk+cos(xik)*sin(xik)/(xk*xk)+ganmak*sin(xik)/xk+xk*cos(xik)/2);
                Kganma4=deltas*(uk*uk/2.-sin(xik)*sin(xik)/(2.*xk*xk)+xk*sin(xik)+sigmap);
                Kx4=deltas*(cos(xik));
                
                xinew=xi+Kxi1/6.+Kxi2/3.+Kxi3/3.+Kxi4/6.;
                unew=u+Ku1/6.+Ku2/3.+Ku3/3.+Ku4/6.;
                ganmanew=ganma+Kganma1/6.+Kganma2/3.+Kganma3/3.+Kganma4/6.;
                xnew=x+Kx1/6.+Kx2/3.+Kx3/3.+Kx4/6.;
            }
            /*................................................................*/
            #undef PI
