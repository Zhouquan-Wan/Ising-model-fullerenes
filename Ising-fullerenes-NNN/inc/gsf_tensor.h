#ifndef GSF_TENSOR_H
#define GSF_TENSOR_H

#include <stdint.h>
#include <ap_int.h>
using UINT_t = ap_uint<128> ;
template <class T>
void swap(T &a, T &b)
{
    T c(a);
    a = b;
    b = c;
}

inline uint64_t si(uint64_t s, int i)
{
    return (s & (1ULL << i)) >> i;
}

inline uint64_t lbitroll(uint64_t s, int d, int rank)
{
    return ((s << d) | (s >> (rank - d))) & ((1ULL << rank) - 1);
}

inline uint64_t rbitroll(uint64_t s, int d, int rank)
{
    return ((s >> d) | (s << (rank - d))) & ((1ULL << rank) - 1);
}

class GSF_Tensor
{
public:
    int rank;
    int n_inset;
    UINT_t *degeneracy;
    UINT_t *degeneracy_tmp;
    double *e_min;
    double *e_min_tmp;
    uint64_t* gs,*gs_tmp;
    GSF_Tensor(int rank, int max_rank);
    GSF_Tensor(){};
    void rebuild(int rank);
    void init();
    void lroll(int ind);
    void rroll(int ind);
    void contraction(GSF_Tensor *itensor, GSF_Tensor *otensor, int *bond_type, int *bond_idx);
    void poly_init(double J, int d);

private:

    uint64_t D_MAX;
    uint64_t D;
};

#endif