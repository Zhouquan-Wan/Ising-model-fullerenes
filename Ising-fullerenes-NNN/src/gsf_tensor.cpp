#include <cstring>
#include <stdint.h>
#include <float.h>
#include <ap_int.h>
#include <omp.h>
#include "gsf_tensor.h"


void GSF_Tensor::init()
{
    #pragma omp parallel for
    for (uint64_t i = 0; i < (1ULL << rank); i++)
    {
        degeneracy[i] = 0;
        // degeneracy_tmp[i] = 0;
        e_min[i] = DBL_MAX;
        // e_min_tmp[i] = DBL_MAX;
        gs[i] = 0;
        // gs_tmp[i] = 0;
    }
}

GSF_Tensor::GSF_Tensor(int rank, int max_rank)
{
    GSF_Tensor::rank = rank;
    GSF_Tensor::n_inset = 0;
    GSF_Tensor::D = 1ULL << rank;
    GSF_Tensor::D_MAX = (1ULL << max_rank);
    e_min = (double *)malloc((1ULL << max_rank) * sizeof(double));
    e_min_tmp = (double *)malloc((1ULL << max_rank) * sizeof(double));
    degeneracy = (UINT_t *)malloc((1ULL << max_rank) * sizeof(UINT_t));
    degeneracy_tmp = (UINT_t *)malloc((1ULL << max_rank) * sizeof(UINT_t));
    gs = (uint64_t *)malloc((1ULL << max_rank) * sizeof(uint64_t));
    gs_tmp = (uint64_t *)malloc((1ULL << max_rank) * sizeof(uint64_t));

    init();
}

void GSF_Tensor::rebuild(int rank)
{
    GSF_Tensor::rank = rank;
    GSF_Tensor::D = 1ULL << rank;
    n_inset = 0;
    if (D_MAX < D)
    {
        std::cout << "realloc!" << std::endl;
        free(degeneracy_tmp);
        free(degeneracy);
        free(e_min);
        free(e_min_tmp);
        free(gs);
        free(gs_tmp);
        D_MAX = D;
        e_min = (double *)malloc(D*sizeof(double));
        e_min_tmp = (double *)malloc(D*sizeof(double));
        degeneracy = (UINT_t *)malloc(D*sizeof(UINT_t));
        degeneracy_tmp = (UINT_t *)malloc(D*sizeof(UINT_t));
        gs = (uint64_t *)malloc(D * sizeof(uint64_t));
        gs_tmp = (uint64_t *)malloc(D * sizeof(uint64_t));
    }
    init();
}



void GSF_Tensor::lroll(int d)
{
    swap<UINT_t *>(degeneracy, degeneracy_tmp);
    swap<double *>(e_min, e_min_tmp);
    swap<uint64_t *>(gs, gs_tmp);
    #pragma omp parallel for
    for (uint64_t i = 0; i < (1ULL << rank); i++)
    {
        degeneracy[i] = degeneracy_tmp[lbitroll(i, d, rank)];
        e_min[i] = e_min_tmp[lbitroll(i, d, rank)];
        gs[i] = gs_tmp[lbitroll(i, d, rank)];
    }
}

void GSF_Tensor::rroll(int d)
{
    swap(degeneracy, degeneracy_tmp);
    swap(e_min, e_min_tmp);
    swap<uint64_t *>(gs, gs_tmp);
    #pragma omp parallel for
    for (uint64_t i = 0; i < (1ULL << rank); i++)
    {
        degeneracy[i] = degeneracy_tmp[rbitroll(i, d, rank)];
        e_min[i] = e_min_tmp[rbitroll(i, d, rank)];
        gs[i] = gs_tmp[rbitroll(i, d, rank)];
    }
}

void GSF_Tensor::contraction(GSF_Tensor *itensor, GSF_Tensor *otensor, int *bond_type, int *bond_idx)
{
    int rank_change = 0;
    int n_bond_in = 0;
    int rank_rest = itensor->rank;
    int o_n_inset = itensor->n_inset;
    for (int i = 0; i < rank; i++)
    {
        rank_change += bond_type[i];
        if (bond_type[i] != 1)
        {
            n_bond_in++;
            rank_rest--;
        }
        if (bond_type[i] == -1)
        {
            o_n_inset++;
        }
    }
    otensor->rebuild(itensor->rank + rank_change);
    otensor->n_inset = o_n_inset;
    for (uint64_t i = 0; i < (1ULL << rank); i++)
    {
        uint64_t oidx_pre = 0;
        uint64_t iidx_pre = 0;
        int shift_res = 0;
        for (int j = rank - 1; j > -1; j--)
        {
            if (bond_type[j] != -1)
            {
                oidx_pre += si(i, rank - j - 1) << shift_res;
                shift_res++;
            }
            if (bond_type[j] != 1)
            {
                iidx_pre += si(i, rank - 1 - j) << (n_bond_in - bond_idx[j] - 1);
            }
        }
        
        oidx_pre = oidx_pre << (rank_rest);
        iidx_pre = iidx_pre << (rank_rest);
        #pragma omp parallel for
        for (uint64_t rest_idx = 0; rest_idx < (1ULL << (rank_rest)); rest_idx++)
        {
            uint64_t iidx = iidx_pre + rest_idx;
            uint64_t oidx = oidx_pre + rest_idx;
            double E_min_try = itensor->e_min[iidx] + e_min[i];
            if (std::abs(E_min_try - otensor->e_min[oidx]) < 1e-11)
            {
                // otensor->degeneracy[oidx] += itensor->degeneracy[iidx] * degeneracy[i];
                otensor->degeneracy[oidx] += itensor->degeneracy[iidx];
            }
            else if (E_min_try < otensor->e_min[oidx])
            {
                otensor->e_min[oidx] = E_min_try;
                // otensor->degeneracy[oidx] = itensor->degeneracy[iidx] * degeneracy[i];
                otensor->degeneracy[oidx] = itensor->degeneracy[iidx];
                int shift = itensor->n_inset;
                otensor->gs[oidx] = itensor->gs[iidx];
                for (int xx= 0; xx < rank; xx++)
                {
                    if (bond_type[xx]==-1)
                    {
                        otensor->gs[oidx] += (si(i,rank-1-xx)<<shift);
                        shift++;
                    }
                }
            }
        }
    }
}

void GSF_Tensor::poly_init(double J, int d)
{
    for (int i = 0; i < 1 << rank; i++)
    {
        e_min[i] = 0;
        for (int j = 0; j < rank; j++)
        {
            if (si(i, j) == si(i, (j + d) % rank))
            {
                e_min[i] += J;
            }
            else
            {
                e_min[i] -= J;
            }
        }
        degeneracy[i] = 1;
    }
}