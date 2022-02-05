#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const int nx = 400;
const int ny = 90;
const int i1 = 200;
const int i2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 0.1;
const double xA = 0.45;
const double yA = 0.45;
const int IT_MAX = 10000;

void VFI(double v_x[nx + 1][ny + 1], double v_y[nx + 1][ny + 1], double psi[nx + 1][ny + 1], FILE* fp_1, FILE* fp_2)
{
    for(int i = 1; i <= nx - 1; i++)
    {
        for(int j = 1; j <= ny - 1; j++)
        {
            v_x[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2. * delta);
            v_y[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2. * delta);
        }
    }

    for(int i = i1; i <= i2; i++)
    {
        for(int j = 0; j <= j_1; j++)
        {
            v_x[i][j] = 0.;
            v_y[i][j] = 0.;
        }
    }

    for(int i = 1; i <= nx - 1; i++)
    {
        v_x[i][0] = 0.;
        v_y[i][ny] = 0.;
    }

    for(int j = 0; j <= ny; j++)
    {
        v_x[0][j] = v_x[1][j];
        v_x[nx][j] = v_x[nx - 1][j];
    }

    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++)
        {
            fprintf(fp_1, "%d \t %d \t %g \n", i, j, v_x[i][j]);
            fprintf(fp_2, "%d \t %d \t %g \n", i, j, v_y[i][j]);
        }
        fprintf(fp_1, "\n");
        fprintf(fp_2, "\n");
    }
}

void F_U(double u_0[nx+ 1][ny + 1])
{
    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++)
        {
            u_0[i][j] = 1./(2. * M_PI * sigma * sigma) * exp( -(pow((delta * i) - xA, 2) + pow((delta * j) - yA, 2))/(2. * sigma * sigma));
        }
    }
}

void AD(double u_0[nx + 1][ny + 1], double u_1[nx + 1][ny + 1], double v_x[nx + 1][ny + 1], double v_y[nx + 1][ny + 1], double psi[nx + 1][ny + 1], FILE*  fp_1, FILE* fp_2, double D, double delta_t)
{
    for(int it = 0; it <= IT_MAX; it++)
    {
        for(int i = 0; i <= nx; i++)
        {
            for(int j = 0; j <= ny; j++)
            {
                u_1[i][j] = u_0[i][j];
            }
        }

        for(int k = 1; k <= 20; k++)
        {
            for(int i = 0; i <= nx; i++)
            {
                for(int j = 1; j <= ny - 1; j++)
                {
                    if(i < i1 || i > i2 || j > j_1)
                    {
                        if(i == 0)
                        {
                            u_1[i][j] = (1. / (1. + ((2. * D * delta_t) / (delta * delta)))) * (u_0[i][j] - (delta_t / 2.) * v_x[i][j] * (((u_0[i + 1][j] - u_0[nx][j])/(delta * 2.)) + (u_1[i + 1][j] - u_1[nx][j])/ (2. * delta)) - (delta_t / 2.) * v_y[i][j] *((u_0[i][j + 1] - u_0[i][j - 1])/(2. * delta) + (u_1[i][j + 1] - u_1[i][j - 1])/(2. * delta)) + (delta_t / 2.) * D *((u_0[i + 1][j] + u_0[nx][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4. * u_0[i][j])/(delta * delta) + (u_1[i + 1][j] + u_1[nx][j] + u_1[i][j + 1] + u_1[i][j - 1])/(delta * delta)));
                        }
                        else if(i == nx)
                        {
                            u_1[i][j] = (1. / (1. + ((2. * D * delta_t) / (delta * delta)))) * (u_0[i][j] - (delta_t / 2.) * v_x[i][j] * (((u_0[0][j] - u_0[i - 1][j])/(delta * 2.)) + (u_1[0][j] - u_1[i - 1][j])/(2. * delta)) - (delta_t / 2.) * v_y[i][j] * ((u_0[i][j + 1] - u_0[i][j - 1])/(2. * delta) + (u_1[i][j + 1] - u_1[i][j - 1])/(2. * delta)) + (delta_t /2.) * D * ((u_0[0][j] + u_0[i - 1][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4. * u_0[i][j])/(delta * delta) + (u_1[0][j] + u_1[i - 1][j] + u_1[i][j + 1] + u_1[i][j - 1])/(delta * delta)));
                        }
                        else
                        {
                            u_1[i][j] = (1. / (1. +((2. * D * delta_t ) / (delta * delta)))) * (u_0[i][j] - (delta_t / 2.) * v_x[i][j] * (((u_0[i + 1][j] - u_0[i - 1][j])/(delta * 2.)) + (u_1[i + 1][j] - u_1[i - 1][j])/(2. * delta)) - (delta_t / 2.) * v_y[i][j] * ((u_0[i][j + 1] - u_0[i][j - 1])/(2. * delta) + (u_1[i][j + 1] - u_1[i][j - 1])/(2. * delta)) + (delta_t / 2.) * D * ((u_0[i + 1][j] + u_0[i - 1][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4. * u_0[i][j])/(delta * delta) + (u_1[i + 1][j] + u_1[i - 1][j] + u_1[i][j + 1] + u_1[i][j - 1])/(delta * delta)));
                        }
                    }
                }
            }
        }

        if(it == 4000)
        {
            for(int i = 0; i <= nx; i++)
            {
                for(int j = 0; j <= ny; j++)
                {
                    fprintf(fp_1, "%d \t %d \t %g \n", i, j, u_1[i][j]);
                }
                fprintf(fp_1, "\n");
            }
        }

        for(int i = 0; i <= nx; i++)
        {
            for(int j = 0; j <= ny; j++)
            {
                u_0[i][j] = u_1[i][j];
            }
        }

        double c = 0.;
        double x_sr = 0.;

        for(int i = 0; i <= nx; i++)
        {
            for(int j = 0; j <= ny; j++)
            {
                c += delta * delta * u_0[i][j];
                x_sr += delta * delta * delta * i * u_0[i][j];
            }
        }

        fprintf(fp_2, "%g \t %g \t %g \n", delta_t * it, c, x_sr);

    }
}


int main()
{
    FILE* vx = fopen("vx.txt", "w");
    FILE* vy = fopen("vy.txt", "w");
    FILE* integral = fopen("integral1.txt", "w");
    FILE* read = fopen("psi.dat", "r");

    double psi[nx + 1][ny + 1];
    double v_x[nx + 1][ny + 1];
    double v_y[nx + 1][ny + 1];
    double u_0[nx + 1][ny + 1];
    double u_1[nx + 1][ny + 1];

    double temp;
    int i, j;

    while(fscanf(read, "%d \t %d \t %lf", &i, &j, &temp) == 3)
    {
        psi[i][j] = temp;
    }

    fclose(read);

    VFI(v_x, v_y, psi, vx, vy);
    fclose(vx);
    fclose(vy);

    double maximum = 0.;
    for(int ii = 0; ii <= nx; ii++)
    {
        for(int jj = 0; jj <= ny; jj++)
        {
            if(maximum < sqrt(v_x[ii][jj] * v_x[ii][jj] + v_y[ii][jj] * v_y[ii][jj]))
            {
                maximum = sqrt(v_x[ii][jj] * v_x[ii][jj] + v_y[ii][jj] * v_y[ii][jj]);
            }
        }
    }

    double delta_t = delta / (4 * maximum);

    //FILE* map1_1 = fopen("map1_1.txt", "w");
    //FILE* map1_2 = fopen("map1_2.txt", "w");
    //FILE* map1_3 = fopen("map1_3.txt", "w");
    //FILE* map1_4 = fopen("map1_4.txt", "w");
    FILE* map1_5 = fopen("map1_5.txt", "w");

    F_U(u_0);

    AD(u_0, u_1, v_x, v_y, psi, map1_5, integral, 0., delta_t);

    //fclose(map1_1);
    //fclose(map1_2);
    //fclose(map1_3);
    //fclose(map1_4);
    fclose(map1_5);
    fclose(integral);

    FILE* integral1 = fopen("integral2.txt", "w");

    //FILE* map2_1 = fopen("map2_1.txt", "w");
    //FILE* map2_2 = fopen("map2_2.txt", "w");
    //FILE* map2_3 = fopen("map2_3.txt", "w");
    //FILE* map2_4 = fopen("map2_4.txt", "w");
    FILE* map2_5 = fopen("map2_5.txt", "w");


    F_U(u_0);

    AD(u_0, u_1, v_x, v_y, psi, map2_5, integral1, 0.1, delta_t);

    //fclose(map2_1);
    //fclose(map2_2);
    //fclose(map2_3);
    //fclose(map2_4);
    fclose(map2_5);
    fclose(integral1);

    return 0;
}
