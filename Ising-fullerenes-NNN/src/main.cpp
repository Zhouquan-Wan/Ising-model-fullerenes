#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <cstring>
#include <cstdlib>
#include "ap_int.h"
#include "gsf_tensor.h"
#include "fullerenes.h"
#include <bitset>
#include <omp.h>
int main(int argc, char *argv[])
{
    double itime, ftime, exec_time;
    itime = omp_get_wtime();
    // fullerenes("C60-Ih.xyz");

    double J1 = std::atof(argv[1]);
    double J2 = std::atof(argv[2]);
    std::cout << "J1 = " << J1 << std::endl
              << "J2 = " << J2 << std::endl;
    GSF_Tensor tri(3, 3);
    tri.poly_init(J1, 1);
    GSF_Tensor penta(5, 5);
    penta.poly_init(J2, 2);
    GSF_Tensor hex(6, 6);
    hex.poly_init(J2, 2);

    GSF_Tensor A(5, 21);
    GSF_Tensor B(5, 21);
    GSF_Tensor *pA = &A;
    GSF_Tensor *pB = &B;
    GSF_Tensor *tmp;
    A.poly_init(J2, 2);

    int bond_type[6] = {0, 1, 0, 0, 0, 0};
    int bond_idx[6] = {0, 0, 1, 0, 0, 0};

    for (int i = 0; i < 5; i++)
    {
        tri.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);
        pA->lroll(1);
    }

    bond_type[0] = 0;
    bond_type[1] = 1;
    bond_type[2] = 1;
    bond_type[3] = 1;
    bond_type[4] = 0;
    bond_type[5] = -1;
    bond_idx[0] = 0;
    bond_idx[1] = 0;
    bond_idx[2] = 0;
    bond_idx[3] = 0;
    bond_idx[4] = 2;
    bond_idx[5] = 1;
    for (int i = 0; i < 5; i++)
    {
        hex.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);
        pA->lroll(2);
    }

    pA->rroll(1);
    bond_type[0] = 0;
    bond_type[1] = -1;
    bond_type[2] = 0;
    bond_type[3] = 0;
    bond_type[4] = 0;
    bond_type[5] = 0;
    bond_idx[0] = 0;
    bond_idx[1] = 1;
    bond_idx[2] = 2;
    bond_idx[3] = 0;
    bond_idx[4] = 0;
    bond_idx[5] = 0;
    for (int i = 0; i < 5; i++)
    {
        tri.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);
        pA->lroll(4);
    }

    for (int i = 0; i < 5; i++)
    {
        pA->lroll(1);
        bond_type[0] = 0;
        bond_type[1] = 1;
        bond_type[2] = 0;
        bond_type[3] = 0;
        bond_type[4] = 0;
        bond_type[5] = 0;
        bond_idx[0] = 0;
        bond_idx[1] = 0;
        bond_idx[2] = 1;
        bond_idx[3] = 0;
        bond_idx[4] = 0;
        bond_idx[5] = 0;
        tri.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);

        pA->lroll(1);
        bond_type[0] = 0;
        bond_type[1] = 1;
        bond_type[2] = 1;
        bond_type[3] = 0;
        bond_type[4] = -1;
        bond_type[5] = 0;
        bond_idx[0] = 0;
        bond_idx[1] = 0;
        bond_idx[2] = 0;
        bond_idx[3] = 2;
        bond_idx[4] = 1;
        bond_idx[5] = 0;
        penta.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);

        pA->lroll(1);
        bond_type[0] = 0;
        bond_type[1] = 0;
        bond_type[2] = -1;
        bond_type[3] = 0;
        bond_type[4] = 0;
        bond_type[5] = 0;
        bond_idx[0] = 0;
        bond_idx[1] = 2;
        bond_idx[2] = 1;
        bond_idx[3] = 0;
        bond_idx[4] = 0;
        bond_idx[5] = 0;
        tri.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);
    }

    pA->rroll(2);
    for (int i = 0; i < 5; i++)
    {

        pA->lroll(1);
        bond_type[0] = 0;
        bond_type[1] = 1;
        bond_type[2] = 0;
        bond_type[3] = 0;
        bond_type[4] = 0;
        bond_type[5] = 0;
        bond_idx[0] = 0;
        bond_idx[1] = 0;
        bond_idx[2] = 1;
        bond_idx[3] = 0;
        bond_idx[4] = 0;
        bond_idx[5] = 0;
        tri.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);

        pA->lroll(2);
        bond_type[0] = 0;
        bond_type[1] = 1;
        bond_type[2] = 1;
        bond_type[3] = 0;
        bond_type[4] = -1;
        bond_type[5] = -1;
        bond_idx[0] = 0;
        bond_idx[1] = 0;
        bond_idx[2] = 0;
        bond_idx[3] = 3;
        bond_idx[4] = 2;
        bond_idx[5] = 1;
        hex.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);

        pA->lroll(1);
        bond_type[0] = 0;
        bond_type[1] = 0;
        bond_type[2] = -1;
        bond_type[3] = 0;
        bond_type[4] = 0;
        bond_type[5] = 0;
        bond_idx[0] = 0;
        bond_idx[1] = 2;
        bond_idx[2] = 1;
        bond_idx[3] = 0;
        bond_idx[4] = 0;
        bond_idx[5] = 0;
        tri.contraction(pA, pB, bond_type, bond_idx);
        swap<GSF_Tensor *>(pA, pB);
    }

    double e_min = __DBL_MAX__;
    ap_int<128> degeneracy = 0;
    uint64_t gs[3] = {0, 0, 0};
    for (int i = 0; i < 1ULL << pA->rank; i++)
    {
        int j = rbitroll(i, 2, pA->rank);
        double e_min_try = pA->e_min[i] + pA->e_min[j];
        if (std::abs(e_min_try - e_min) < 1e-11)
        {
            degeneracy += pA->degeneracy[i] * pA->degeneracy[j];
        }
        else if (e_min_try < e_min)
        {
            e_min = e_min_try;
            degeneracy = pA->degeneracy[i] * pA->degeneracy[j];
            gs[0] = pA->gs[i];
            gs[1] = i;
            gs[2] = pA->gs[j];
        }
    }
    std::cout << "Ground state energy is " << e_min << std::endl
              << "Degeneracy is " << degeneracy << std::endl;
    std::cout << "Configuration is " << std::bitset<35>(gs[0]) << " " << std::bitset<20>(gs[1]) << "  " << std::bitset<35>(gs[2]) << std::endl;

    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    std::cout << std::endl
              << "Time taken is " << exec_time << std::endl;
}
