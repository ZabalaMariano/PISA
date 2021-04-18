#pragma once

#include "./partitioned_sequence_enumerator_opt_vb.hpp"
#include "configuration.hpp"//#include "./configuration_opt_vb.hpp"
#include "./global_parameters_opt_vb.hpp"
#include "../codec/integer_codes.hpp"
#include "./util_opt_vb.hpp"
#include "./indexed_sequence_opt_vb.hpp"
#include "./compact_ranked_bitvector_opt_vb.hpp"
#include "./optimizer_opt_vb.hpp"

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
                          global_parameters_opt_vb const& params,
                          uint64_t& dense_short, uint64_t& dense_medium, uint64_t& dense_large,
                          uint64_t& sparse_short, uint64_t& sparse_medium, uint64_t& sparse_large,
                          uint64_t& dense_short_cost, uint64_t& dense_medium_cost, uint64_t& dense_large_cost,
                          uint64_t& sparse_short_cost, uint64_t& sparse_medium_cost, uint64_t& sparse_large_cost,
                          uint64_t& cantidad_integers_con_interpolative,
                          uint64_t& cantidad_integers_con_varintg8iu, bool dense_sparse)
        {
            assert(n > 0);
            auto partition = optimizer_opt_vb<VByteBlockType, VByteBlockType2>::compute_partition(begin, n);
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
                            n, params, n,
                            dense_short, dense_medium, dense_large,
                            sparse_short, sparse_medium, sparse_large,
                            dense_short_cost, dense_medium_cost, dense_large_cost,
                            sparse_short_cost, sparse_medium_cost, sparse_large_cost,
                            cantidad_integers_con_interpolative,
                            cantidad_integers_con_varintg8iu, dense_sparse);
                
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
                                curr_n, params, n,
                                dense_short, dense_medium, dense_large,
                                sparse_short, sparse_medium, sparse_large,
                                dense_short_cost, dense_medium_cost, dense_large_cost,
                                sparse_short_cost, sparse_medium_cost, sparse_large_cost,
                                cantidad_integers_con_interpolative,
                                cantidad_integers_con_varintg8iu, dense_sparse);

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
                                global_parameters_opt_vb const& params, uint64_t size_posting_list,
                                uint64_t& dense_short, uint64_t& dense_medium, uint64_t& dense_large,
                                uint64_t& sparse_short, uint64_t& sparse_medium, uint64_t& sparse_large,
                                uint64_t& dense_short_cost, uint64_t& dense_medium_cost, uint64_t& dense_large_cost,
                                uint64_t& sparse_short_cost, uint64_t& sparse_medium_cost, uint64_t& sparse_large_cost,
                                uint64_t& cantidad_integers_con_interpolative,
                                uint64_t& cantidad_integers_con_varintg8iu, bool dense_sparse)
        {
            assert(n > 0);
            switch (type) {
                case VBBlock::type:
                {
                    uint64_t cost_encoder1 = VBBlock::bitsize(begin, params, universe, n);

                    if(dense_sparse){
                        if (size_posting_list < 10000) {
                            dense_short+=n;
                            dense_short_cost+=cost_encoder1;
                        } else if (size_posting_list < 7000000) {
                            dense_medium+=n;
                            dense_medium_cost+=cost_encoder1;
                        } else {
                            dense_large+=n;
                            dense_large_cost+=cost_encoder1;
                        }
                    }

                    bvb.append_bits(type, type_bits);
                    push_pad(bvb);
                    VBBlock::write(bvb, begin, base, universe, n, params);
                    break;
                    }
                case RBBlock::type:
                {
                    if (!std::is_same<RBBlock, pvb::compact_ranked_bitvector_opt_vb>::value){
                        uint64_t cost_encoder2 = RBBlock::bitsize(begin, params, universe, n);

                        if(dense_sparse){
                            if (size_posting_list < 10000) {
                                sparse_short+=n;
                                sparse_short_cost+=cost_encoder2;
                            } else if (size_posting_list < 7000000) {
                                sparse_medium+=n;
                                sparse_medium_cost+=cost_encoder2;
                            } else {
                                sparse_large+=n;
                                sparse_large_cost+=cost_encoder2;
                            }

                            if(n<8){
                                cantidad_integers_con_interpolative += n;
                            } else {
                                cantidad_integers_con_varintg8iu += n;
                            }
                        }

                        bvb.append_bits(type, type_bits);
                        push_pad(bvb);
                        RBBlock::write(bvb, begin, base, universe, n, params);
                    } else {
                        bvb.append_bits(type, type_bits);
                        RBBlock::write(bvb, begin, base, universe, n, params);   
                    }
                    break;
                }
                default:
                {
                    assert(false);
                }
            }
        }
    };
}
