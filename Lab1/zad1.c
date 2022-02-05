#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void euler(FILE *f,double krok);
void RK2(FILE *f, double krok);
void RK4(FILE *f, double krok);
void RRZ(FILE *f, double w_podane);
double V(double w, double t);


int main(void)
{
    double* krok = malloc(3*sizeof(double));
    krok[0]=0.01;
    krok[1]=0.1;
    krok[2]=1.0;


    // ----------- Euler-----------
    
    FILE *f1 = fopen("euler1.dat","w");
    euler(f1,krok[0]);
    FILE *f2 = fopen("euler2.dat","w");
    euler(f2,krok[1]);
    FILE *f3 = fopen("euler3.dat","w");
    euler(f3,krok[2]);
    fclose(f1);
    fclose(f2);
    fclose(f3);

    // ----------- RK2 --------------
    
    f1 = fopen("RK2_1.dat","w");
    RK2(f1,krok[0]);
    f2 = fopen("RK2_2.dat","w");
    RK2(f2,krok[1]);
    f3 = fopen("RK2_3.dat","w");
    RK2(f3,krok[2]);
    fclose(f1);
    fclose(f2);
    fclose(f3);

    // ----------- RK4 --------------

    f1 = fopen("RK4_1.dat","w");
    RK4(f1,krok[0]);
    f2 = fopen("RK4_2.dat","w");
    RK4(f2,krok[1]);
    f3 = fopen("RK4_3.dat","w");
    RK4(f3,krok[2]);
    fclose(f1);
    fclose(f2);
    fclose(f3);

    // ----------- RRZ 2rzedu -----------

    f1 = fopen("RRZ_1.dat","w");
    RRZ(f1,0.5);
    f2 = fopen("RRZ_2.dat","w");
    RRZ(f2,0.8);
    f3 = fopen("RRZ_3.dat","w");
    RRZ(f3,1.0);
    FILE *f4 = fopen("RRZ_4.dat","w");
    RRZ(f4,1.2);
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);


}

void euler(FILE *f,double krok)
{   

    double lambda = -1, tmin = 0, tmax = 5;
    double* y = malloc(2*sizeof(double)); 
    y[0]=1.0;
    y[1]=1.0;
        for(double j=tmin; j<tmax; j+=krok)
        {
            y[0]=y[1];
            fprintf(f,"%12.6g\t%12.6g\t%12.6g\n",j,y[0]-exp(j*lambda),y[0]);
            y[1]=y[0]+krok*lambda*y[0];
        }
}

void RK2(FILE *f, double krok)
{
    double lambda = -1, tmin = 0, tmax = 5, y0 =1.0,y1=1.0,k1,k2;
    for(double i=tmin; i<tmax; i+=krok)
    {
        y0=y1;
        fprintf(f,"%12.6g\t%12.6g\t%12.6g\n",i,y0-exp(i*lambda),y0);
        k1 = lambda*y0;
        k2 = lambda*(y0+krok*k1);
        y1 = y0 +(krok/2.0)*(k1+k2);
    }
}

void RK4(FILE *f, double krok)
{
    double lambda = -1, tmin = 0, tmax = 5, y0 =1.0,y1=1.0,k1,k2,k3,k4;

    for(double i=tmin; i<tmax; i+=krok)
    {
        y0=y1;
        fprintf(f,"%12.6g\t%12.6g\t%12.6g\n",i,y0-exp(i*lambda),y0);
        k1 = lambda*y0;
        k2 = lambda*(y0+krok*k1/2.0);
        k3 = lambda*(y0+krok*k2/2.0);
        k4 = lambda*(y0+krok*k3);
        y1 = y0+(krok/6.0)*(k1+2*k2+2*k3+k4);
    }

}

double V(double w, double t)
{
    return 10*sin(w*t);
}

void RRZ(FILE *f, double w_podane)
{
    double delta = pow(10,-4), R=100.0, L = 0.1, C = 0.001, w0 = 1./(sqrt(L*C));
    double T0 = 2*M_PI/w0, tmax = 4.0*T0;
    double Q=0, I=0,Q1=0,I1=0;
    double k1q, k1i,k2q,k2i,k3q,k3i,k4q,k4i;
    double wv = w_podane*w0;

    for(double i=0; i<=tmax; i+=delta)
    {
        Q=Q1;
        I=I1;
        fprintf(f,"%12.6g\t%12.6g\t%12.6g\n",i,Q,I);
        k1q=I;
        k1i=V(wv,i)/L - Q/(L*C) - R*I/L;
        k2q = I + k1i*delta/2;
        k2i = V(wv,i + delta/2)/L - (Q + delta*k1q/2)/(L*C) - (I+k1i*delta/2)*R/L;
        k3q = I + k2i*delta/2;
        k3i = V(wv, i + delta/2)/L - (Q - k2q*delta/2)/(L*C) - (I+k2i*delta/2)*R/L;
        k4q = I+ delta*k3i;
        k4i = V(wv, i + delta/2)/L - (Q - k3q*delta)/(L*C) - (I+k3i*delta)*R/L;

        Q1 = Q + (delta/6)*(k1q+2*k2q+2*k3q+k4q);
        I1 = I + (delta/6)*(k1i+2*k2i+2*k3i+k4i);

    }

}