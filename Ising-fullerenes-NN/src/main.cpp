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
#include <omp.h>
int main(int argc, char *argv[])
{
    double itime, ftime, exec_time;
    itime = omp_get_wtime();
    std::ifstream ifile(argv[1], std::ios::in);
    int N;
    ifile >> N;
    int *roll = (int *)malloc(sizeof(int) * N);
    int *con = (int *)malloc(sizeof(int) * N);
    int _rank = 3, max_rank = 0;
    for (int i = 0; i < N - 1; i++)
    {
        ifile >> roll[i] >> con[i];
        _rank += 3 - 2 * con[i];
        if (_rank > max_rank)
        {
            max_rank = _rank;
        }
    }
    ifile.close();
    GSF_Tensor A(3, max_rank);
    GSF_Tensor *pA = &A;
    A.tri_init();
    GSF_Tensor tri(3, 4);
    tri.tri_init();
    uint64_t bond_type[3] = {0, 1, 0};
    uint64_t bond_idx[3] = {0, 0, 1};
    for (int i = 0; i < N - 1; i++)
    {
        if (roll[i] != 0)
        {
            pA->lroll(roll[i]);
        }
        for (int j = 0; j < 3; j++)
        {
            bond_type[j] = 1;
            bond_idx[j] = 0;
        }
        for (int j = 0; j < con[i]; j++)
        {
            bond_type[j] = -1;
            bond_idx[j] = j;
        }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        std::cout << "Number of points: " << i + 1 << "  Tensor rank: " << pA->rank << "   Time taken: " << exec_time << std::endl;
        tri.contraction(pA, bond_type, bond_idx);
    }

    std::cout << "Ground state degeneracy is " << pA->degeneracy[0] << std::endl;
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    std::cout << "Time taken is " << exec_time << std::endl;
}
