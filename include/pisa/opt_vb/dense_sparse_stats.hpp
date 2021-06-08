#pragma once

namespace pvb {
    class stats{
        public:
        stats(bool b, bool interpolative)
            : dense_sparse(b), interpolative(interpolative),
            cantidad_integers_con_interpolative(0),
            cantidad_integers_sin_interpolative(0),
            dense_short(0), dense_medium(0), dense_large(0), 
            sparse_short(0), sparse_medium(0), sparse_large(0),
            dense_short_cost(0), dense_medium_cost(0), dense_large_cost(0), 
            sparse_short_cost(0), sparse_medium_cost(0), sparse_large_cost(0)
        {}
            uint64_t cantidad_integers_con_interpolative,
            cantidad_integers_sin_interpolative,
            dense_short, dense_medium, dense_large, 
            sparse_short, sparse_medium, sparse_large,
            dense_short_cost, dense_medium_cost, dense_large_cost, 
            sparse_short_cost, sparse_medium_cost, sparse_large_cost;
            bool dense_sparse, interpolative;
    };
}