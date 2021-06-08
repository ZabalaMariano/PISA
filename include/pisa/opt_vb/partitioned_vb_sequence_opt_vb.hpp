#pragma once

#include "./partitioned_sequence_enumerator_opt_vb.hpp"
#include "configuration.hpp"//#include "./configuration_opt_vb.hpp"
#include "./global_parameters_opt_vb.hpp"
#include "../codec/integer_codes.hpp"
#include "./util_opt_vb.hpp"
#include "./indexed_sequence_opt_vb.hpp"
#include "./compact_ranked_bitvector_opt_vb.hpp"
#include "./optimizer_opt_vb.hpp"
#include "./opt_vb/dense_sparse_stats.hpp"

#include <limits>
#include <cmath>
#include <type_traits>

namespace pvb {

    template<typename VByteBlockType, typename VByteBlockType2, typename Enumerator>
    struct partitioned_vb_sequence_opt_vb {

        typedef Enumerator enumerator;

        typedef VByteBlockType  VBBlock;
        typedef VByteBlockType2 RBBlock;

        template<typename Iterator>
        static void write(pisa::bit_vector_builder& bvb,
                          Iterator begin,
                          uint64_t universe, uint64_t n,
                          global_parameters_opt_vb const& params, stats& stats)
        {
            assert(n > 0);
            auto partition = optimizer_opt_vb<VByteBlockType, VByteBlockType2>::compute_partition(begin, n, stats);
            size_t partitions = partition.size();
            assert(partitions > 0);

            pisa::write_gamma_nonzero(bvb, partitions);

            if (partitions == 1) {                
                auto const& singleton = partition.front();
                uint64_t base = *begin;
                uint64_t universe_bits = ceil_log2(universe);
                bvb.append_bits(base, universe_bits);
                uint64_t relative_universe = *(begin + n - 1) - base;
                // write universe only if non-singleton and not tight
                if (n > 1) {
                    if (base + relative_universe + 1 == universe) {
                        pisa::write_delta(bvb, 0); // tight universe
                    } else {
                        pisa::write_delta(bvb, relative_universe);
                    }
                }

                push_pad(bvb);
                write_block(bvb, begin, singleton.type, base,
                            *(begin + n - 1) + 1,
                            n, params, n, stats);
                
            } else {
                pisa::bit_vector_builder bv_sequences;
                std::vector<uint64_t> sizes;
                sizes.reserve(partitions);
                std::vector<uint64_t> endpoints;
                endpoints.reserve(partitions);
                std::vector<uint64_t> upper_bounds;
                upper_bounds.reserve(partitions + 1);
                upper_bounds.push_back(*begin);

                auto b = begin;
                for (uint64_t prev_size = 0, i = 0; i < partitions; ++i) {
                    uint64_t curr_base = i == 0
                           ? *(begin)
                           : *(partition[i - 1].begin - 1) + 1;
                    uint64_t curr_n = std::distance(b, partition[i].begin);
                    uint64_t curr_universe = *(b + curr_n - 1);

                    write_block(bv_sequences,
                                b,
                                partition[i].type,
                                curr_base,
                                curr_universe + 1,
                                curr_n, params, n, stats);

                    sizes.push_back(prev_size + curr_n);
                    endpoints.push_back(bv_sequences.size());
                    upper_bounds.push_back(curr_universe);
                    prev_size = sizes.back();
                    b = partition[i].begin;
                }

                assert(upper_bounds.back() == *(begin + n - 1));
                assert(sizes.front() != 0);
                assert(sizes.back() == n);
                assert(sizes.size() == partitions);
                assert(endpoints.size() == partitions);
                assert(upper_bounds.size() == partitions + 1);

                pisa::bit_vector_builder bv_sizes;
                pisa::bit_vector_builder bv_upper_bounds;

                compact_elias_fano_opt_vb::write(bv_sizes, sizes.begin(),
                                          n, partitions - 1, params);
                compact_elias_fano_opt_vb::write(bv_upper_bounds, upper_bounds.begin(),
                                          universe, partitions + 1, params);
                uint64_t endpoint_bits = ceil_log2(bv_sequences.size() + 1);
                pisa::write_gamma(bvb, endpoint_bits);

                bvb.append(bv_sizes);
                bvb.append(bv_upper_bounds);

                for (uint64_t p = 0; p < endpoints.size() - 1; ++p) {
                    bvb.append_bits(endpoints[p], endpoint_bits);
                }

                push_pad(bvb);
                bvb.append(bv_sequences);
            }
        }

    private:

        static const uint64_t type_bits = indexed_sequence_opt_vb<
            compact_elias_fano_opt_vb,
            block_sequence_opt_vb<VByteBlockType2>
        >::type_bits;

        template<typename Iterator>
        static void write_block(pisa::bit_vector_builder& bvb,
                                Iterator begin, int type,
                                uint64_t base, uint64_t universe, uint64_t n,
                                global_parameters_opt_vb const& params, uint64_t size_posting_list, stats& stats)
        {
            assert(n > 0);
            switch (type) {
                case VBBlock::type:
                {
                    if(stats.dense_sparse){
                        uint64_t cost_encoder1 = VBBlock::bitsize(begin, params, universe, n);
                        stats.cantidad_integers_sin_interpolative += n;
                        if (size_posting_list < 10000) {
                            stats.dense_short+=n;
                            stats.dense_short_cost+=cost_encoder1;
                        } else if (size_posting_list < 7000000) {
                            stats.dense_medium+=n;
                            stats.dense_medium_cost+=cost_encoder1;
                        } else {
                            stats.dense_large+=n;
                            stats.dense_large_cost+=cost_encoder1;
                        }
                    }

                    bvb.append_bits(type, type_bits);
                    push_pad(bvb);
                    VBBlock::write(bvb, begin, base, universe, n, params, size_posting_list);
                    break;
                }
                case RBBlock::type:
                {
                    if (!std::is_same<RBBlock, pvb::compact_ranked_bitvector_opt_vb>::value){
                        
                        if(stats.dense_sparse){
                            uint64_t cost_encoder2 = RBBlock::bitsize(begin, params, universe, n);
                            stats.cantidad_integers_sin_interpolative += n;
                            if (size_posting_list < 10000) {
                                stats.sparse_short+=n;
                                stats.sparse_short_cost+=cost_encoder2;
                            } else if (size_posting_list < 7000000) {
                                stats.sparse_medium+=n;
                                stats.sparse_medium_cost+=cost_encoder2;
                            } else {
                                stats.sparse_large+=n;
                                stats.sparse_large_cost+=cost_encoder2;
                            }
                        }

                        bvb.append_bits(type, type_bits);
                        push_pad(bvb);
                        RBBlock::write(bvb, begin, base, universe, n, params, size_posting_list);
                    } else {
                        bvb.append_bits(type, type_bits);
                        RBBlock::write(bvb, begin, base, universe, n, params, size_posting_list);   
                    }
                    break;
                }
                case 3://*
                {
                    if(stats.dense_sparse){
                        stats.cantidad_integers_con_interpolative += n;
                    }

                    bvb.append_bits(type, type_bits);
                    push_pad(bvb);
                    interpolative_opt_vb::write(bvb, begin, base, universe, n, params, size_posting_list);
                }
                default:
                {
                    assert(false);
                }
            }
        }
    };
}
