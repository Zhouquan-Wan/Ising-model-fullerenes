#include <cstring>
#include <stdint.h>
#include <float.h>
#include <ap_int.h>
#include <omp.h>
#include "gsf_tensor.h"

void GSF_Tensor::init()
{
    if (rank > 0)
    {
#pragma omp parallel for
        for (uint64_t i = 0; i < (1ULL << (rank - 1)); i++)
        {
            degeneracy[i] = 0;
        }
    }
    else
    {
        degeneracy[0] = 0;
    }
}

GSF_Tensor::GSF_Tensor(uint64_t rank, uint64_t max_rank)
{
    GSF_Tensor::rank = rank;
    GSF_Tensor::D = 1ULL << rank;
    GSF_Tensor::D_MAX = (1ULL << max_rank);
    degeneracy = (UINT_t *)malloc((1ULL << (max_rank - 1)) * sizeof(UINT_t));
    degeneracy_tmp = (UINT_t *)malloc((1ULL << (max_rank - 1)) * sizeof(UINT_t));
    init();
}

void GSF_Tensor::rebuild(uint64_t rank)
{
    GSF_Tensor::rank = rank;
    if (rank != 0)
    {
        GSF_Tensor::D = 1ULL << (rank - 1);
    }
    if (D_MAX < D)
    {
        std::cout << "realloc!" << std::endl;
        free(degeneracy_tmp);
        free(degeneracy);
        D_MAX = D;
        degeneracy = (UINT_t *)malloc(D / 2 * sizeof(UINT_t));
        degeneracy_tmp = (UINT_t *)malloc(D / 2 * sizeof(UINT_t));
    }
    init();
}

void GSF_Tensor::lroll(uint64_t d)
{
    swap<UINT_t *>(degeneracy, degeneracy_tmp);
#pragma omp parallel for
    for (uint64_t i = 0; i < (1ULL << (rank - 1)); i++)
    {
        degeneracy[i] = degeneracy_tmp[sbar(lbitroll(i, d, rank), rank)];
    }
}

void GSF_Tensor::rroll(uint64_t d)
{
    swap<UINT_t *>(degeneracy, degeneracy_tmp);
#pragma omp parallel for
    for (uint64_t i = 0; i < (1ULL << (rank - 1)); i++)
    {
        degeneracy[i] = degeneracy_tmp[sbar(rbitroll(i, d, rank), rank)];
    }
}

void GSF_Tensor::contraction(GSF_Tensor *itensor, uint64_t *bond_type, uint64_t *bond_idx)
{
    int rank_change = 0;
    int n_bond_in = 0;
    int rank_rest = itensor->rank;
    for (uint64_t i = 0; i < rank; i++)
    {
        rank_change += bond_type[i];
        if (bond_type[i] != 1)
        {
            n_bond_in++;
            rank_rest--;
        }
    }

    int o_rank = itensor->rank + rank_change;
    if (o_rank == 0)
    {
        o_rank = 1;
    }
#pragma omp parallel for
    for (uint64_t i = 0; i < (1ULL << (o_rank - 1)); i++)
    {
        itensor->degeneracy_tmp[i] = 0;
    }
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
            if ((degeneracy[i] == 1) & (si(oidx, o_rank - 1) == 0))
            {
                itensor->degeneracy_tmp[oidx] += itensor->degeneracy[sbar(iidx, itensor->rank)];
            }
        }
    }
    swap<UINT_t *>(itensor->degeneracy, itensor->degeneracy_tmp);
    itensor->rank = o_rank;
    itensor->D = 1ULL << o_rank;
}

void GSF_Tensor::tri_init()
{
    degeneracy[0b111] = 0;
    degeneracy[0b001] = 1;
    degeneracy[0b010] = 1;
    degeneracy[0b100] = 1;
    degeneracy[0b110] = 1;
    degeneracy[0b101] = 1;
    degeneracy[0b011] = 1;
    degeneracy[0b000] = 0;
}