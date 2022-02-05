#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double beta = 0.001;
const double N = 500;
const double y = 0.1;
const double tmax = 100;
const double dt = 0.1;
const double u0 = 1;
const double TOL = pow(10,-6);
double alpha = beta*N-y;

double a11 = 1/4;
double a12 = 1/4 - sqrt(3)/6;
double a21 = 1/4 + sqrt(3)/6;
double a22 = 1/4;

void Picard(FILE *f);
void IterNewton(FILE *f);

double F1(double U1, double U2, double u_n);
double F2(double U1, double U2, double u_n);
double m11(double U1);
double m21(double U1);
double m12(double U1);
double m22(double U1);
double f(double U);
void RK2(FILE *f);

int main(void)
{
    FILE *f1 = fopen("Picard.dat","w");
    Picard(f1);
    fclose(f1);
    FILE *f2 = fopen("IterNewton.dat","w");
    IterNewton(f2);
    fclose(f2);
    FILE *f3 = fopen("RK2.dat","w");
    RK2(f3);
    fclose(f3);
}

void Picard(FILE *f)
{
    double u1 = u0, mi = 0, mi_max=20;;
    double u_0 = u0, u_mi0=0;
    for(double i=0; i<=tmax; i+=dt)
    {   
        fprintf(f,"%lf\t%lf\t%lf\n", i, u_0, N-u_0);

        while(1)
        {
            u_mi0=u1;
            mi++;
            u1 = u_0 + (dt/2)*(alpha*u_0-beta*u_0*u_0+alpha*u_mi0-beta*(u_mi0*u_mi0));
            if(abs(u1-u_mi0)<TOL ||  mi>mi_max)
            {
                break;
            }
             
        }
        u_0=u1;
        u_mi0=0;
        mi=0;
    }
}

void IterNewton(FILE *f)
{
    double u1 = u0, mi = 0, mi_max=20;;
    double u_0 = u0, u_mi0=0;
    for(double i=0; i<=tmax; i+=dt)
    {   
        fprintf(f,"%lf\t%lf\t%lf\n", i, u_0, N-u_0);
        while(1)
        {
            mi++;
            u1 =u_mi0 -(u_mi0 -u_0 - (dt/2)*(alpha*u_0-beta*u_0*u_0+alpha*u_mi0-beta*(u_mi0*u_mi0)))/(1-(dt/2)*(alpha - 2*beta*u_mi0));
            if(abs(u1-u_mi0)<TOL ||  mi>mi_max)
            {
                break;
            }
             u_mi0=u1;
        }
        u_0=u1;
        u_mi0=0;
        mi=0;
    }
}

double F1(double U1, double U2, double u_n)
{   
    return U1-u_n-dt*(a11*(alpha*U1-beta*pow(U1,2))+a12*(alpha*U2-beta*pow(U2,2)));
}

double F2(double U1, double U2, double u_n)
{   
    return U2-u_n-dt*(a21*(alpha*U1-beta*pow(U1,2))+a22*(alpha*U2-beta*pow(U2,2)));
}

double m11(double U1)
{
    return 1- dt*a11*(alpha-2*beta*U1);
}

double m12(double U1)
{
    return -dt*a12*(alpha-2*beta*U1);
}

double m21(double U1)
{
    return -dt*a21*(alpha-2*beta*U1);
}

double m22(double U1)
{
    return 1- dt*a22*(alpha-2*beta*U1);
}

double f(double U)
{
    return (beta*N - y)*U - beta*pow(U,2);
}

void RK2(FILE *f1)
{
    double u1 = u0, mi = 0, mi_max=20;
    double u_0 = u0, u_mi0=0;

    double U1,U2,dU1,dU2;

    for(double i=0; i<=tmax; i+=dt)
    {   
        U1=U2=u_0;
        fprintf(f1,"%lf\t%lf\t%lf\n", i, u_0, N-u_0);
        while(1)
        {
            mi++;
            dU1 = (F2(U1,U2,u_0)*m12(U2)-F1(U1,U2,u_0)*m22(U2))/(m11(U1)*m22(U2) - m12(U2)*m21(U1));
            dU1 = (F1(U1,U2,u_0)*m21(U1)-F2(U1,U2,u_0)*m11(U1))/(m11(U1)*m22(U2) - m12(U2)*m21(U1));
            U1 = U1 + dU1;
            U2 = U2 + dU2;

            u1 = u_0 + dt*(f(U1)/2+f(U2)/2);
            if(abs(u1-u_mi0)<TOL ||  mi>mi_max)
            {
                break;
            }
             u_mi0=u1;
        }
        u_0=u1;
        u_mi0=0;
        mi=0;
    }
}