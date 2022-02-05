#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double dt = 0.2;
const int nx = 128;
const int ny = 128;
const double xmax = 128.*dt;
const double ymax = 128.*dt;
const double TOL = pow(10,-8);

void metoda_relaksacji_wielosiatkowej(int k, double V[nx+1][ny+1], FILE *fS, FILE *fV)
{
    double S_it=0., S_it_1;
    int it = 0;

    do
    {
        for(int i=k; i<=nx-k; i+=k)
        {
            for(int j=k; j<=ny-k; j+=k)
         {
            V[i][j] = 0.25 * (V[i + k][j] + V[i - k][j] + V[i][j + k] + V[i][j - k]);
         }
        }

        S_it_1 = S_it;
        S_it = 0.;
        for(int i=0; i<=nx-k; i+=k)
        {
            for(int j=0; j<=ny-k; j+=k)
         {
            S_it += (0.5 * k * k * dt * dt) * (pow(((V[i + k][j] - V[i][j])/(2 * k * dt) + (V[i + k][j + k] - V[i][j + k])/(2 * k * dt)), 2) + pow(((V[i][j + k] - V[i][j]) / (2 * k * dt) + (V[i + k][j + k] - V[i + k][j]) / (2 * k * dt)), 2));
         }
        }
        
        fprintf(fS,"%d \t %lf\n", it, S_it);
        it++;

    }while (fabs((S_it-S_it_1)/S_it_1)>=TOL);
    
    for (int i =0; i <= nx; i += k)
    {
        for (int j = 0; j <= ny; j += k)
        {
            fprintf(fV, "%f \t %f \t %f \n", i * dt, j * dt, V[i][j]);
        }
    }

    if(k!=1)
    {
      for(int i=0; i<=nx-k; i+=k)
        {
            for(int j=0; j<=ny-k; j+=k)
         {
            V[i + k / 2][j + k / 2] = 0.25*(V[i][j] + V[i + k][j] + V[i][j + k] + V[i + k][j + k]);
            if (i != nx - k)
            {
                V[i + k][j + k / 2] = 0.5*(V[i + k][j] + V[i + k][j + k]);
            }

            if (j != ny - k)
            {
                V[i + k / 2][j + k] = 0.5*(V[i][j + k] + V[i + k][j + k]);
            }

            if (j != 0)
            {
                V[i + k / 2][j] = 0.5*(V[i][j] + V[i + k][j]);
            }

            if (i != 0)
            {
                V[i][j + k / 2] = 0.5*(V[i][j] + V[i][j + k]);
            }
         }

        }   
    }

}

int main (void)
{
        double V[nx + 1][ny + 1];

    for(int i = 0; i <= nx; i++)
    {
        for(int j =0; j <= ny; j++)
        {
            V[i][j] = 0.;
        }
    }

    for(int i = 0; i <= ny; i++)
    {
        V[0][i] = sin(M_PI * i * dt / ymax);
        V[nx][i] = sin(M_PI * i * dt / ymax);
        V[i][0] = sin(2. * M_PI * i * dt/ xmax);
        V[i][ny] = -sin(2. * M_PI * i * dt / xmax);
    }

    FILE* multi_grid_S_16 = fopen("multi_grid_S_16.txt", "w");
    FILE* multi_grid_V_16 = fopen("multi_grid_V_16.txt", "w");
    metoda_relaksacji_wielosiatkowej(16, V, multi_grid_S_16, multi_grid_V_16);
    fclose(multi_grid_S_16);
    fclose(multi_grid_V_16);

    FILE* multi_grid_S_8 = fopen("multi_grid_S_8.txt", "w");
    FILE* multi_grid_V_8 = fopen("multi_grid_V_8.txt", "w");
    metoda_relaksacji_wielosiatkowej(8, V, multi_grid_S_8, multi_grid_V_8);
    fclose(multi_grid_S_8);
    fclose(multi_grid_V_8);

    FILE* multi_grid_S_4 = fopen("multi_grid_S_4.txt", "w");
    FILE* multi_grid_V_4 = fopen("multi_grid_V_4.txt", "w");
    metoda_relaksacji_wielosiatkowej(4, V, multi_grid_S_4, multi_grid_V_4);
    fclose(multi_grid_S_4);
    fclose(multi_grid_V_4);

    FILE* multi_grid_S_2 = fopen("multi_grid_S_2.txt", "w");
    FILE* multi_grid_V_2 = fopen("multi_grid_V_2.txt", "w");
    metoda_relaksacji_wielosiatkowej(2, V, multi_grid_S_2, multi_grid_V_2);
    fclose(multi_grid_S_2);
    fclose(multi_grid_V_2);

    FILE* multi_grid_S_1 = fopen("multi_grid_S_1.txt", "w");
    FILE* multi_grid_V_1 = fopen("multi_grid_V_1.txt", "w");
    metoda_relaksacji_wielosiatkowej(1, V, multi_grid_S_1, multi_grid_V_1);
    fclose(multi_grid_S_1);
    fclose(multi_grid_V_1);
}