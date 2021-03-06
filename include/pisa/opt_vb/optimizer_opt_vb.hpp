#pragma once

#include "configuration.hpp"
#include "./opt_vb/dense_sparse_stats.hpp"

namespace pvb {

    template<typename Encoder, typename Encoder2>
    struct optimizer_opt_vb {

        typedef Encoder     encoder1;
        typedef Encoder2    encoder2;

        static const int invalid_block_type = -1;

        template<typename Iterator>
        struct block {
            block(Iterator b)
                : type(invalid_block_type)
                , begin(b)
            {}

            int type;
            Iterator begin;
        };

        template<typename Iterator>
        static std::vector<block<Iterator>>
        compute_partition(Iterator begin, uint64_t n, stats& stats)//configuration_opt_vb const& conf)
        {
            assert(encoder1::type == 0 and
                   encoder2::type == 1);

            auto curr_begin = begin;   // begin of current block
            auto curr_bv_end = begin;  // end of bin. vect. block
            auto curr_vb_end = begin;  // end of VByte block
            auto end = begin + n;      // end of list
            auto it = begin;           // current position

            std::vector<block<Iterator>> partition;
            block<Iterator> curr_block(curr_begin);
            int last_block_type = invalid_block_type;

            int64_t best_bv_gain = 0;
            int64_t best_vb_gain = 0;
            int64_t curr_gain = -1;
            int64_t F = int64_t(pisa::configuration::get().fix_cost);
            int64_t T = F;

            posting_type curr;
            posting_type last = *begin;

            auto push_block = [&]()
            {
                last_block_type = curr_block.type;
                partition.push_back(curr_block);
                curr_block = block<Iterator>(curr_begin);
                T = 2 * F;
            };

            auto encode_encoder2 = [&]()
            {
                curr_block.type = encoder2::type;
                curr_begin = curr_bv_end;
                curr_block.begin = curr_begin;
                push_block();
                curr_gain -= best_bv_gain;
                best_bv_gain = 0;
                best_vb_gain = curr_gain;
                curr_vb_end = it + 1;
                assert(curr_gain <= 0);
            };

            auto encode_encoder1 = [&]()
            {
                curr_block.type = encoder1::type;
                curr_begin = curr_vb_end;
                curr_block.begin = curr_begin;
                push_block();
                curr_gain -= best_vb_gain;
                best_vb_gain = 0;
                best_bv_gain = curr_gain;
                curr_bv_end = it + 1;
                assert(curr_gain >= 0);
            };

            auto close = [&]()
            {
                if (best_bv_gain > F and best_bv_gain - curr_gain > F) {
                    encode_encoder2();
                }

                if (best_vb_gain < -F and best_vb_gain - curr_gain < -F) {
                    encode_encoder1();
                }

                if (curr_gain > 0) {
                    curr_bv_end = end;
                    encode_encoder2();
                } else {
                    curr_vb_end = end;
                    encode_encoder1();
                }
            };
        
            if(stats.interpolative && n<128){
                curr_block.type = 3;//*//encoder1::type;
                curr_block.begin = end;
                partition.push_back(curr_block);
            }else{
                for (int64_t last_gain = 0; it != end; ++it,
                    last = curr, last_gain = curr_gain)
                {
                    curr = *it;
                    curr_gain += encoder1::posting_cost(curr, last)
                            - encoder2::posting_cost(curr, last);

                    if (curr_gain >= last_gain) { // gain is not decreasing

                        if (curr_gain > best_bv_gain) {
                            best_bv_gain = curr_gain;
                            curr_bv_end = it + 1;
                        }

                        assert(best_vb_gain <= curr_gain);
                        if (best_vb_gain - curr_gain < -2 * F and best_vb_gain < -T) {
                            encode_encoder1();
                        }

                    } else { // gain is decreasing

                        if (curr_gain <= best_vb_gain) {
                            best_vb_gain = curr_gain;
                            curr_vb_end = it + 1;
                        }

                        assert(best_bv_gain >= curr_gain);
                        if (best_bv_gain - curr_gain > 2 * F and best_bv_gain > T) {
                            encode_encoder2();
                        }
                    }
                }

                assert(it == end);
                close();
            }
            return partition;
        }
    };
}