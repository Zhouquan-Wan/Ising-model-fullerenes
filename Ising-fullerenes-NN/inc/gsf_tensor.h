#ifndef GSF_TENSOR_H
#define GSF_TENSOR_H

using UINT_t = ap_uint<800>;
// using UINT_t = long double ;
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

inline uint64_t sbar(uint64_t s, int rank)
{
    if (si(s, rank - 1) == 1)
    {
        return ((1ULL << rank) - 1) - s;
    }
    else
    {
        return s;
    }
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
    UINT_t *degeneracy;
    UINT_t *degeneracy_tmp;
    GSF_Tensor(uint64_t rank, uint64_t max_rank);
    GSF_Tensor(){};
    void rebuild(uint64_t rank);
    void init();
    void lroll(uint64_t ind);
    void rroll(uint64_t ind);
    void contraction(GSF_Tensor *itensor, uint64_t *bond_type, uint64_t *bond_idx);
    void tri_init();

private:
    uint64_t D_MAX;
    uint64_t D;
};

#endif