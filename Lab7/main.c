#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

double delta = 0.01;
double ro = 1.;
double mi = 1.;
int nx = 200;
int ny = 90;
int i1 = 50;
int j_1 = 55;
int IT_MAX = 20000;


void WBpsi(double Q_we, double Q_wy, double psi[nx + 1][ny + 1], double y[ny + 1])
{
    for(int j = j_1; j <= ny; j++)
    {
        psi[0][j] = Q_we/(2 * mi) * (pow(y[j], 3)/3. - pow(y[j], 2)/2. * (y[j_1] + y[ny]) + y[j] * y[j_1] * y[ny]);
    }

    for(int j = 0; j <= ny; j++)
    {
        psi[nx][j] = Q_wy/ (2 * mi) * (pow(y[j], 3)/3. - pow(y[j], 2)/2. * y[ny]) + (Q_we * pow(y[j_1], 2) * (-y[j_1] + 3 * y[ny]))/ (12 * mi);
    }

    for(int i = 1; i <= nx - 1; i++)
    {
        psi[i][ny] = psi[0][ny];
    }

    for(int i = i1; i <= nx - 1; i++)
    {
        psi[i][0] = psi[0][j_1];
    }

    for(int j = 1; j <= j_1; j++)
    {
        psi[i1][j] = psi[0][j_1];
    }

    for(int i = 1; i <= i1; i++)
    {
        psi[i][j_1] = psi[0][j_1];
    }
}

void WBzeta(double Q_we, double Q_wy, double psi[nx + 1][ny + 1], double zeta[nx + 1][ny + 1], double y[ny + 1])
{
    for(int j = j_1; j <= ny; j++)
    {
        zeta[0][j] = Q_we/ (2 * mi) * (2 * y[j] - y[j_1] - y[ny]);
    }

    for(int j = 0; j <= ny; j++)
    {
        zeta[nx][j] = Q_wy / (2 * mi) * (2 * y[j] - y[ny]);
    }

    for(int i = 1; i <= nx - 1; i++)
    {
        zeta[i][ny] = 2. / (delta * delta) * (psi[i][ny - 1] - psi[i][ny]);
    }

    for(int i = i1 + 1; i <= nx -1; i++)
    {
        zeta[i][0] = 2. / (delta * delta) * (psi[i][1] - psi[i][0]);
    }

    for(int j = 1; j <= j_1 - 1; j++)
    {
        zeta[i1][j] = 2. / (delta * delta) * (psi[i1 + 1][j] - psi[i1][j]);
    }

    for(int i = 1; i <= i1; i++)
    {
        zeta[i][j_1] = 2. / (delta * delta) * (psi[i][j_1 + 1] - psi[i][j_1]);
    }

    zeta[i1][j_1] = (zeta[i1 - 1][j_1] + zeta[i1][j_1 - 1])/2.;
}


void nav_stokes(double Q_we, FILE* fp)
{
    double psi[nx + 1][ny + 1];
    double zeta[nx + 1][ny + 1];
    double y[ny + 1];

    for(int i = 0; i <= nx ; i++)
    {
        for(int j = 0; j <= ny; j++)
        {
            psi[i][j] = 0;
            zeta[i][j] = 0;
        }
    }

    for(int j = 0; j <= ny; j++)
        {
            y[j] = delta * j;
        }

    double Q_wy = Q_we * (pow(y[ny], 3) - pow(y[j_1], 3) - 3 * y[j_1] * pow(y[ny], 2) + 3 * pow(y[j_1], 2) * y[ny]) / pow(y[ny], 3);

    int omega;
    double error_control;
    int j_2 = j_1 + 2;

    WBpsi(Q_we, Q_wy, psi, y);

    for(int it = 1; it < IT_MAX; it++)
    {
        if(it < 2000)
        {
            omega = 0;
        }
        else
        {
            omega = 1;
        }

        for(int i = 1; i <= nx - 1; i++)
        {
            for(int j = 1; j <= ny - 1; j++)
            {
                if(!((i == 0) || (i == nx) || (j == 0) || (j == ny) || (i <= i1 && j == j_1) || (i == i1 && j <= j_1))&& (!(i <= i1 && j < j_1)))
                {
                    psi[i][j] = (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - zeta[i][j] * pow(delta, 2))/4.;
                    zeta[i][j] = (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1])/4. - omega * ro / (16.0 * mi) * ((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j]) - (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1]));
                }
            }
        }

        WBzeta(Q_we, Q_wy, psi, zeta, y);
        error_control = 0.;

        for(int i = 1; i <= nx - 1; i++)
        {
            error_control += (psi[i + 1][j_2] + psi[i - 1][j_2] + psi[i][j_2 + 1] + psi[i][j_2 - 1] - 4 * psi[i][j_2] - delta * delta * zeta[i][j_2]);
        }

        printf("It: %d \t Error_control: %f \n", it, error_control);

    }

    double u;
    double v;

    for(int i = 1; i <= nx - 1; i++)
    {
        for(int j = 1; j <= ny - 1; j++)
        {
            if(!((i == 0) || (i == nx) || (j == 0) || (j == ny) || (i <= i1 && j == j_1) || (i == i1 && j <= j_1))&& (!(i <= i1 && j < j_1)))
            {
                u = (psi[i][j + 1] - psi[i][j - 1])/ (2 * delta);
                v = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * delta);
            }
            else
            {
                u = 0;
                v = 0;
            }

            fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \n", i * delta, j* delta, psi[i][j], zeta[i][j], u, v);
        }
        fprintf(fp, "\n");
    }

}

int main(void)
{
    FILE* Q_minus_1000 = fopen("Q_minus_1000.txt", "w");
    nav_stokes(-1000, Q_minus_1000);
    fclose(Q_minus_1000);

    FILE* Q_minus_4000 = fopen("Q_minus_4000.txt", "w");
    nav_stokes(-4000, Q_minus_4000);
    fclose(Q_minus_4000);

    FILE* Q_4000 = fopen("Q_4000.txt", "w");
    nav_stokes(4000, Q_4000);
    fclose(Q_4000);
}