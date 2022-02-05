#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double epsilon = 1.;
const double dt = 0.1;
const int nx = 150;
const int ny = 100;
const double V1 = 10.;
const double V2 = 0.;
const double xmax = 15.;
const double ymax = 10.;
const double sigx = 1.5;
const double sigy = 1.;
const double TOL = pow(10, -8);

double gestosc_1(double x, double y)
{
    return exp((-pow((x - 0.35 * xmax),2))/(pow(sigx, 2)) - (pow(y - 0.5 * ymax, 2))/(pow(sigy, 2)));
}

double gestosc_2(double x, double y)
{
    return -exp((-pow((x - 0.65 * xmax),2))/(pow(sigx, 2)) - (pow(y - 0.5 * ymax, 2))/(pow(sigy, 2)));
}

double gestosc(double x, double y)
{
    return gestosc_1(x,y)+gestosc_2(x,y);
}

void Relaksacja_globalna(double wg, FILE* fp_S, FILE* fp_Vn, FILE* fp_E)
{
    double V_n[nx + 1][ny + 1], V_s[nx + 1][ny + 1], error[nx + 1][ny + 1];
    double g[nx + 1][ny + 1];

    int iter = 0;

    for(int i = 0; i < nx + 1; i++)
    {
        for(int j = 0; j < ny + 1; j++)
        {
            V_n[i][j] = 0.;
            V_s[i][j] = 0.;
            error[i][j] = 0.;
            g[i][j] = gestosc(i * dt, j * dt);
        }
        V_n[i][0] = V1;
        V_s[i][0] = V1;
    }

    double S_it, S_it_1;
    S_it = 0.;
    S_it_1 = 0.;

    do
    {
        for(int i = 1; i < nx; i++)
        {
            for(int j = 1; j < ny; j++)
            {
                V_n[i][j] = 0.25 * (V_s[i + 1][j] + V_s[i - 1][j] + V_s[i][j + 1] + V_s[i][j - 1] + pow(dt, 2)/epsilon * g[i][j]);
            }
        }

        for(int j = 1;  j <= ny ; j++)
        {
            V_n[0][j] = V_n[1][j];
            V_n[nx][j] = V_n[nx - 1][j];
        }

        for(int i = 0; i <= nx ; i++)
        {
            for(int j = 0; j <= ny; j++)
            {
                V_s[i][j] = (1. - wg) * V_s[i][j] + wg * V_n[i][j];
            }
        }

        S_it_1 = S_it;
        S_it = 0.;


        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                S_it +=  pow(dt, 2) * (0.5 * pow(((V_n[i+1][j] - V_n[i][j])/dt), 2) + 0.5 * pow((V_n[i][j+1] - V_n[i][j])/dt, 2) - g[i][j] * V_n[i][j]);
            }
        }

        fprintf(fp_S, "%d \t %f \n", iter, S_it);
        iter++;
    }while(fabs((S_it - S_it_1) / S_it_1) >= TOL);

    for(int i = 1; i < nx; i++)
    {
        for(int j = 1; j < ny; j++)
        {
            fprintf(fp_Vn, "%f \t %f \t %f \n", dt * i, dt * j, V_n[i][j]);
        }
    }

    for(int i = 1; i < nx; i++)
    {
        for(int j = 1; j < ny; j++)
        {
            error[i][j] = (V_n[i+1][j] - 2 * V_n[i][j] + V_n[i - 1][j])/ (dt* dt) + (V_n[i][j + 1] - 2 * V_n[i][j] + V_n[i][j - 1])/ (dt * dt) + g[i][j]/epsilon;
            fprintf(fp_E, "%f \t %f \t %f \n", dt * i, dt * j, error[i][j]);
        }
    }

}

void Relaksacja_lokalna(double omega, FILE* fp_S)
{
    double V_n[nx + 1][ny + 1];
    double g[nx + 1][ny +1];

    int iter = 0;

    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++)
        {
            V_n[i][j] = 0.;
            g[i][j] = gestosc(i * dt, j * dt);
        }
        V_n[i][0] = V1;
    }

    double S_it, S_it_1;
    S_it = 0.;
    S_it_1 = 0.;

    do
    {
        for(int i = 1; i < nx; i++)
        {
            for(int j = 1; j < ny; j++)
            {
                V_n[i][j] = 0.25 * omega * (V_n[i + 1][j] + V_n[i - 1][j] + V_n[i][j + 1] + V_n[i][j - 1] + pow(dt, 2)/epsilon * g[i][j]) + (1. - omega) * V_n[i][j];
            }
        }

        for(int j = 1;  j <= ny; j++)
        {
            V_n[0][j] = V_n[1][j];
            V_n[nx][j] = V_n[nx - 1][j];
        }

        S_it_1 = S_it;
        S_it = 0.;

        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                S_it +=  pow(dt, 2) * (0.5 * pow(((V_n[i+1][j] - V_n[i][j])/dt), 2) + 0.5 * pow((V_n[i][j+1] - V_n[i][j])/dt, 2) - g[i][j] * V_n[i][j]);
            }
        }
        fprintf(fp_S, "%d \t %f \n", iter, S_it);
        iter++;
    }while(fabs((S_it - S_it_1)/S_it_1) >= TOL);

}

int main()
{

    FILE* global_S_06 = fopen("global_S_06.txt", "w");
    FILE* global_V_06 = fopen("global_V_06.txt", "w");
    FILE* global_E_06 = fopen("global_E_06.txt", "w");
    Relaksacja_globalna(0.6, global_S_06, global_V_06, global_E_06);
    fclose(global_S_06);
    fclose(global_V_06);
    fclose(global_E_06);

    FILE* global_S_1 = fopen("global_S_1.txt", "w");
    FILE* global_V_1 = fopen("global_V_1.txt", "w");
    FILE* global_E_1 = fopen("global_E_1.txt", "w");
    Relaksacja_globalna(1., global_S_1, global_V_1, global_E_1);
    fclose(global_S_1);
    fclose(global_V_1);
    fclose(global_E_1);

    FILE* local_1 = fopen("local_S_1.txt", "w");
    Relaksacja_lokalna(1., local_1);
    fclose(local_1);
    FILE* local_14 = fopen("local_S_14.txt", "w");
    Relaksacja_lokalna(1.4, local_14);
    fclose(local_14);
    FILE* local_18 = fopen("local_S_18.txt", "w");
    Relaksacja_lokalna(1.8, local_18);
    fclose(local_18);
    FILE* local_19 = fopen("local_S_19.txt", "w");
    Relaksacja_lokalna(1.9, local_19);
    fclose(local_19);


    return 0;
}