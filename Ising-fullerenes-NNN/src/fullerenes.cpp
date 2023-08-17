#include <fstream>
#include <cstring>
#include <iostream>
#include "fullerenes.h"
fullerenes::fullerenes(std::string filename)
{
    std::ifstream ifile(filename, std::ios::in);
    ifile>>Np;
    std::cout<<"# of C atoms: "<<Np<<std::endl;
    std::string s;
    ifile>>s;
    adj = (int**) malloc(sizeof(int*)*3);
    adj[0] = (int *) malloc(sizeof(int) * Np);
    adj[1] = (int *) malloc(sizeof(int) * Np);
    adj[2] = (int *) malloc(sizeof(int) * Np);
    x = (double *) malloc(sizeof(double) * Np);
    y = (double *) malloc(sizeof(double) * Np);
    z = (double *) malloc(sizeof(double) * Np);
    for (int i = 0; i < Np; i++)
    {

        char tmpc;
        int tmpi;
        ifile>>tmpc>>x[i]>>y[i]>>z[i]>>tmpi>>adj[0][i]>>adj[1][i]>>adj[2][i];
        if (tmpi != i+1)
        {
            std::cout<<"not sorted!"<<std::endl;
        }
        adj[0][i]--;
        adj[1][i]--;
        adj[2][i]--;
    }
    ifile.close();

    
    int * bond_idx_matrix = (int *)malloc(Np * Np * sizeof (int));
    for (int i = 0; i < Np * Np; i++)
    {
        bond_idx_matrix[i] = -1;
    }
    dual_Np = 0;
    for (int i = 0; i < Np; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (bond_idx_matrix[i * Np + adj[j][i]] == -1)
            {
                bond_idx_matrix[i * Np + adj[j][i]] = dual_Np ;
                bond_idx_matrix[adj[j][i] * Np + i] = dual_Np ;
                dual_Np ++;
            }
        }
    }
    std::cout<<"# of dual lattice sites: "<< dual_Np<<std::endl;

    dual_adj = (int**) malloc(sizeof(int*)*4);
    dual_adj[0] = (int *) malloc(sizeof(int) * dual_Np);
    dual_adj[1] = (int *) malloc(sizeof(int) * dual_Np);
    dual_adj[2] = (int *) malloc(sizeof(int) * dual_Np);
    dual_adj[3] = (int *) malloc(sizeof(int) * dual_Np);
    dual_x = (double *) malloc(sizeof(double) * dual_Np);
    dual_y = (double *) malloc(sizeof(double) * dual_Np);
    dual_z = (double *) malloc(sizeof(double) * dual_Np);

    for (int i = 0; i < Np; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (bond_idx_matrix[i * Np + adj[j][i]] == -1)
            {
                bond_idx_matrix[i * Np + adj[j][i]] = dual_Np ;
                
            }
        }
    }
}