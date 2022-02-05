#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double x0 = 0.01;
const double v0 = 0.;
const double dt0 = 1.;
const double S = 0.75;
const double p = 2.;
const double t_max = 40.;
const double alfa = 5.;
const double delta = pow(10, -10);
const double TOL_1 = pow(10, -2);
const double TOL_2 = pow(10, -5);

struct data
{
    double x;
    double v;
};

struct data trapezy(double x_n, double v_n, double dt)
{
    struct data d_n_1;
    double xn_1 = x_n;
    double vn_1 = v_n;
    double dx, dv, F, G;
    double a_11;
    double a_12; 
    double a_21; 
    double a_22;
    do
    {
        F = xn_1 - x_n - dt / 2. * (v_n + vn_1);
        G = vn_1 - v_n - dt / 2. * ((alfa * (1 - pow(x_n, 2)) * v_n - x_n) + (alfa * (1 - pow(xn_1, 2)) * vn_1 - xn_1));
        a_11 = 1.;
        a_12 = -(dt / 2.);
        a_21 = -(dt / 2.) * (-2 * alfa * xn_1 * vn_1 - 1);
        a_22 = 1 - (dt / 2.) * alfa * (1 - pow(xn_1, 2));
        
        dx = (-1. * F * a_22 - (-1. * G) * a_12) / (a_11 * a_22 - a_12 * a_21);
        dv = (a_11 * -1. * G - a_21 * -1. * F) / (a_11 * a_22 - a_12 * a_21);

        xn_1 += dx;
        vn_1 += dv;

    } while ((fabs(dv) > delta) && (fabs(dx) > delta));

    d_n_1.x = xn_1;
    d_n_1.v = vn_1;

    return d_n_1;
}

struct data RK2(double x_n, double v_n, double dt)
{
    struct data data_n_1;
    double k1x = v_n;
    double k1v = alfa * (1 - pow(x_n, 2)) * v_n - x_n;
    double k2x = v_n + dt * k1v;
    double k2v = alfa * (1- pow((x_n + dt * k1x), 2)) * (v_n + dt * k1v) - (x_n + dt * k1x);
    data_n_1.x = x_n + dt/2. * (k1x + k2x);
    data_n_1.v = v_n + dt/2. * (k1v + k2v);

    return data_n_1;
}

void kc(FILE* fp, struct data (* f)(double x_n, double v_n, double dt), double TOL)
{
    double t = 0.;
    double dt = dt0;
    double x_n = x0;
    double v_n = v0;
    do
    {
        struct data d2_n_1 = f(x_n, v_n, dt);
        d2_n_1 = f(d2_n_1.x, d2_n_1.v, dt);
        struct data d1_n_1 = f(x_n, v_n, dt * 2);

        double Ex = (d2_n_1.x- d1_n_1.x)/ (pow(2., p) - 1.);
        double Ev = (d2_n_1.v - d1_n_1.v) / (pow(2., p) - 1.);

        if(fmax(fabs(Ex), fabs(Ev)) < TOL)
        {
            t += 2. * dt;
            x_n = d2_n_1.x;
            v_n = d2_n_1.v;
            fprintf(fp, "%f \t %f \t %f \t %f \n", t, dt, x_n, v_n);
        }
        dt *= pow(S * TOL/(fmax(fabs(Ex), fabs(Ev))), (1./(p + 1.)));
    }while(t < t_max);
}

int main()
{
    FILE *f1 = fopen("RK2TOL1.txt", "w");
    kc(f1, RK2, TOL_1);
    fclose(f1);
    FILE *f2 = fopen("RK2TOL2.txt", "w");
    kc(f2, RK2, TOL_2);
    fclose(f2);
    FILE *f3 = fopen("trapezyTOL1.txt", "w");
    kc(f3, trapezy, TOL_1);
    fclose(f3);
    FILE *f4 = fopen("trapezyTOL2.txt", "w");
    kc(f4, trapezy, TOL_2);
    fclose(f4);

    return 0;
}
