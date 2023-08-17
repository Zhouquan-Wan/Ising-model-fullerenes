#ifndef FULLERENES_H
#define FULLERENES_H

class fullerenes
{
public:
    int Np,dual_Np; // # of points 
    int** adj,**dual_adj;
    double *x, *y, *z, *dual_x, *dual_y, *dual_z;
    fullerenes(std::string filename);
};

#endif