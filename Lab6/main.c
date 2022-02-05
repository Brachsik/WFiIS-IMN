#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mgmres.h"
#include <stdbool.h>
#include "mgmres.c"

double delta = 0.1;
int e1=1, e2=1;
int v1=10,v3=10,v2=-10,v4=-10;

double ro(double x, double y, double x_max, double y_max)
{
    double sigma=x_max/10;
    
    return exp(-1. * pow(x - 0.25 * x_max, 2.)/pow(sigma, 2.) - pow(y - 0.5 * y_max, 2.)/pow(sigma, 2.)) - exp(-1. * pow(x - 0.75 * x_max, 2.)/pow(sigma, 2.) - pow(y - 0.5 * y_max, 2.)/pow(sigma, 2.));
}

int reindex_j(int l, int nx)
{
    return floor(l/(nx+1.));
}

int reindex_i(int l, int nx)
{
    return l - reindex_j(l, nx) * (nx + 1);
}

int El(int nx, int i, int e1, int e2)
{
    if(i <= nx/2)
    {
        return e1;
    }
    else
    {
        return e2;
    }
}

void poisson(double e1, double e2, double nx, double ny, double v1, double v2, double v3, double v4, bool g, FILE *f, bool fa, bool fb)
{
    int N =(nx+1)*(ny+1);
    double a[5*N];
    int ja[5*N];
    int ia[N+1];

    for(int i=0; i<N+1; i++)
    {
        ia[i]=-1;
    }


    double b[N];
    double V[N];

    //Algorytm wypełniania macierzy rzadkiej w formacie CSR + WB Dirichleta
    int k=-1;

    for(int l=0; l<N; l++)
    {
        int brzeg=0;
        double vb = 0.;
        if(reindex_i(l, nx) == 0)
        {
            brzeg = 1;
            vb = v1;
        }

        if(reindex_j(l, nx) == ny)
        {
            brzeg = 1;
            vb = v2;
        }

        if(reindex_i(l, nx) == nx)
        {
            brzeg = 1;
            vb = v3;
        }

        if(reindex_j(l, nx) == 0)
        {
            brzeg = 1;
            vb = v4;
        }

        if(!g)
        {
            b[l] = -ro(delta * reindex_i(l, nx), delta * reindex_j(l, nx), delta*nx, delta*ny);
        }
        else
        {
            b[l] = 0.;
        }

        if(brzeg == 1)
        {
            b[l] = vb;
        }

        ia[l]=-1;
        //lewa skrajna przekatna
        if(l - nx - 1 >= 0 && brzeg == 0)
        {
            k++;
            if(ia[l] < 0)
                ia[l] = k;

            a[k] = El(nx, reindex_i(l, nx), e1, e2)/ pow(delta, 2);
            ja[k] = l - nx - 1;
        }

        //poddiagonalna
         if(l - 1 >= 0 && brzeg == 0)
        {
            k++;
            if(ia[l] < 0)
                ia[l] = k;

            a[k] = El(nx, reindex_i(l, nx), e1, e2)/ pow(delta, 2);
            ja[k] = l - 1;
        }

        //diagonalna

        k++;
        if(ia[l] < 0)
            ia[l] = k;
        
        if(brzeg==0)
        a[k]= - (2*El(nx, reindex_i(l, nx), e1, e2)+El(nx, reindex_i(l+1, nx), e1, e2)+El(nx, reindex_i(l+nx+1, nx), e1, e2))/pow(delta,2);
        else
        a[k]=1;

        ja[k]=l;

        //naddiagonalna
        if(l<N && brzeg==0)
        {
            k++;
            a[k]=El(nx, reindex_i(l+1, nx), e1, e2)/pow(delta,2);
            ja[k]=l+1;
        }

        //prawa skrajna przekatna
        if(l<N-nx-1 && brzeg==0)
        {
            k++;
            a[k]=El(nx, reindex_i(l+nx+1, nx), e1, e2)/pow(delta,2);
            ja[k]=l+nx+1;
        }

        if(fb)
        fprintf(f, "%d \t %d \t %d \t %f \n", l, reindex_i(l, nx), reindex_j(l, nx), b[l]);

    }

    int nz_num=k+1;
    ia[N]=nz_num;

    if(fa)
    {
        for(int iter = 0; iter < nz_num; iter++)
        {
        fprintf(f, "%d \t %0.f \n", iter, a[iter]);
        }
    }
    

    if(!fa && !fb)
    {
        int itr_max = 500;
        int mr = 500;
        double tol_abs = 1.0e-8;
        double tol_rel = 1.0e-8;

        pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

        double war= 0;

        for(int iter = 0; iter < N; iter++)
        {
            if(delta * reindex_i(iter, nx) < war)
            {
                fprintf(f, "\n");
            }
        fprintf(f, "%f \t %f \t %lf \n", delta* reindex_i(iter, nx), reindex_j(iter, nx) * delta, V[iter]);
        war = delta *reindex_i(iter, nx);
        }
    }
    

}

int main()
{
    FILE* vect = fopen("vector.txt", "w");

    poisson(1, 1, 4, 4, 10, -10, 10, -10, true, vect, false , true);
    fclose(vect);

    FILE* matr = fopen("matrix.txt", "w");
    poisson(1, 1, 4, 4, 10, -10, 10, -10, true, matr,true,false);
    fclose(matr);

    FILE* no_ro_50 = fopen("map50_50.txt", "w");
    poisson(1, 1, 50, 50, 10, -10, 10, -10, true, no_ro_50,false,false);
    fclose(no_ro_50);

    FILE* no_ro_100 = fopen("map100_100.txt", "w");
    poisson(1, 1, 100, 100, 10, -10, 10, -10, true, no_ro_100,false,false);
    fclose(no_ro_100);

    FILE* no_ro_200 = fopen("map200_200.txt", "w");
    poisson(1, 1, 200, 200, 10, -10, 10, -10, true, no_ro_200,false,false);
    fclose(no_ro_200);

    FILE* ro_1_1 = fopen("map_e_1_1.txt", "w");
    poisson(1, 1, 100, 100, 0, 0, 0, 0, false, ro_1_1,false,false);
    fclose(ro_1_1);

    FILE* ro_1_2 = fopen("map_e_1_2.txt", "w");
    poisson(1, 2, 100, 100, 0, 0, 0, 0, false, ro_1_2,false,false);
    fclose(ro_1_2);

    FILE* ro_1_10 = fopen("map_e_1_10.txt", "w");
    poisson(1, 10, 100, 100, 0, 0, 0, 0, false, ro_1_10,false,false);
    fclose(ro_1_10);

    return 0;
}