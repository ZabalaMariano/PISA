#pragma once

#include "configuration.hpp"

namespace pvb {

    template<typename Encoder, typename Encoder2>
    struct optimizer_opt_vb {

        typedef Encoder     VBBlock;
        typedef Encoder2    RBBlock;

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
        compute_partition(Iterator begin, uint64_t n)//configuration_opt_vb const& conf)
        {
            assert(VBBlock::type == 0 and
                   RBBlock::type == 1);

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

            auto encode_bin_vector = [&]()
            {
                curr_block.type = RBBlock::type;
                curr_begin = curr_bv_end;
                curr_block.begin = curr_begin;
                push_block();
                curr_gain -= best_bv_gain;
                best_bv_gain = 0;
                best_vb_gain = curr_gain;
                curr_vb_end = it + 1;
                assert(curr_gain <= 0);
            };

            auto encode_VByte = [&]()
            {
                curr_block.type = VBBlock::type;
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
                    encode_bin_vector();
                }

                if (best_vb_gain < -F and best_vb_gain - curr_gain < -F) {
                    encode_VByte();
                }

                if (curr_gain > 0) {
                    curr_bv_end = end;
                    encode_bin_vector();
                } else {
                    curr_vb_end = end;
                    encode_VByte();
                }
            };

            for (int64_t last_gain = 0; it != end; ++it,
                 last = curr, last_gain = curr_gain)
            {
                curr = *it;
                curr_gain += VBBlock::posting_cost(curr, last)
                           - RBBlock::posting_cost(curr, last);

                if (curr_gain >= last_gain) { // gain is not decreasing

                    if (curr_gain > best_bv_gain) {
                        best_bv_gain = curr_gain;
                        curr_bv_end = it + 1;
                    }

                    assert(best_vb_gain <= curr_gain);
                    if (best_vb_gain - curr_gain < -2 * F and best_vb_gain < -T) {
                        encode_VByte();
                    }

                } else { // gain is decreasing

                    if (curr_gain <= best_vb_gain) {
                        best_vb_gain = curr_gain;
                        curr_vb_end = it + 1;
                    }

                    assert(best_bv_gain >= curr_gain);
                    if (best_bv_gain - curr_gain > 2 * F and best_bv_gain > T) {
                        encode_bin_vector();
                    }
                }
            }

            assert(it == end);
            close();

            return partition;
        }
    };
}
/*
#pragma once

#include "configuration.hpp"

namespace pvb {

    template<typename Encoder, typename Encoder2>
    struct optimizer_opt_vb {

        typedef Encoder     VBBlock;
        typedef Encoder2    RBBlock;

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
        compute_partition(Iterator begin, uint64_t n)
        {
            std::cout << "\n\nCompute partition" << std::endl;
            assert(VBBlock::type == 0 and
                   RBBlock::type == 1);

            auto curr_begin = begin;   // begin of current block
            auto curr_bv_end = begin;  // end of bin. vect. block
            auto curr_vb_end = begin;  // end of VByte block
            auto end = begin + n;      // end of list
            auto it = begin;           // current position

            int64_t curr_begin_index = 1;
            int64_t curr_vb_end_index = 1;
            int64_t curr_bv_end_index = 1;
            int64_t it_index = 1;

            std::vector<block<Iterator>> partition;
            block<Iterator> curr_block(curr_begin);
            int last_block_type = invalid_block_type;

            int64_t best_bv_gain = 0;
            int64_t best_vb_gain = 0;
            int64_t best_bv_unfix_gain = 0;
            int64_t best_vb_unfix_gain = 0;
            int64_t curr_gain = -1;
            int64_t F = int64_t(pisa::configuration::get().fix_cost);
            int64_t T = F;
            int64_t curr_gain_padding = 0;

            posting_type curr;
            posting_type last = *begin;

            auto push_block = [&]()
            {
                last_block_type = curr_block.type;
                partition.push_back(curr_block);
                curr_block = block<Iterator>(curr_begin);
                T = 2 * F;
            };

            auto encode_bin_vector = [&]()
            {
                std::cout << "\nencode_bin_vector" << std::endl;
                curr_block.type = RBBlock::type;
                std::cout << "curr_begin viejo: " << curr_begin << std::endl;
                std::cout << "curr_bv_end: " << curr_bv_end << std::endl;
                curr_begin = curr_bv_end;
                std::cout << "curr_begin nuevo: " << curr_begin << std::endl;
                std::cout << "curr_begin_index viejo: " << curr_begin_index << std::endl;
                std::cout << "curr_bv_end_index: " << curr_bv_end_index << std::endl;
                curr_begin_index = curr_bv_end_index;
                std::cout << "curr_begin_index nuevo: " << curr_begin_index << std::endl;
                curr_block.begin = curr_begin;
                push_block();
                std::cout << "curr_gain_padding viejo: " << curr_gain_padding << std::endl;
                std::cout << "best_bv_gain: " << best_bv_gain << std::endl;
                curr_gain_padding -= best_bv_gain;
                std::cout << "curr_gain_padding nuevo: " << curr_gain_padding << std::endl;
                std::cout << "curr_gain viejo: " << curr_gain << std::endl;
                std::cout << "best_bv_unfix_gain: " << best_bv_unfix_gain << std::endl;
                curr_gain -= best_bv_unfix_gain;
                std::cout << "curr_gain nuevo: " << curr_gain << std::endl;
                best_bv_gain = 0;
                best_bv_unfix_gain = 0;
                best_vb_gain = curr_gain_padding;
                std::cout << "best_vb_gain: " << best_vb_gain << std::endl;
                curr_vb_end = it + 1;
                curr_vb_end_index = it_index + 1;
                std::cout << "curr_vb_end_index: " << curr_vb_end_index << std::endl;
                assert(curr_gain_padding <= 0);
                std::cout << "\n" << std::endl;
            };

            auto encode_VByte = [&]()
            {
                std::cout << "\nencode_VByte" << std::endl;
                curr_block.type = VBBlock::type;
                std::cout << "curr_begin viejo: " << curr_begin << std::endl;
                std::cout << "curr_vb_end: " << curr_vb_end << std::endl;
                curr_begin = curr_vb_end;
                std::cout << "curr_begin nuevo: " << curr_begin << std::endl;
                std::cout << "curr_begin_index viejo: " << curr_begin_index << std::endl;
                std::cout << "curr_vb_end_index: " << curr_vb_end_index << std::endl;
                curr_begin_index = curr_vb_end_index;
                std::cout << "curr_begin_index nuevo: " << curr_begin_index << std::endl;
                curr_block.begin = curr_begin;
                push_block();
                std::cout << "curr_gain_padding viejo: " << curr_gain_padding << std::endl;
                std::cout << "best_vb_gain: " << best_vb_gain << std::endl;
                curr_gain_padding -= best_vb_gain;
                std::cout << "curr_gain_padding nuevo: " << curr_gain_padding << std::endl;
                std::cout << "curr_gain viejo: " << curr_gain << std::endl;
                std::cout << "best_vb_unfix_gain: " << best_vb_unfix_gain << std::endl;
                curr_gain -= best_vb_unfix_gain;
                std::cout << "curr_gain nuevo: " << curr_gain << std::endl;
                best_vb_gain = 0;
                best_vb_unfix_gain = 0;
                best_bv_gain = curr_gain_padding;
                std::cout << "best_bv_gain: " << best_vb_gain << std::endl;
                curr_bv_end = it + 1;
                curr_bv_end_index = it_index + 1;
                std::cout << "curr_bv_end_index: " << curr_vb_end_index << std::endl;
                assert(curr_gain_padding >= 0);
                std::cout << "\n" << std::endl;
            };

            auto close = [&]()
            {
                if (best_bv_gain > F and best_bv_gain - curr_gain_padding > F) {
                    std::cout << "close: encode_bin_vector" << std::endl;
                    encode_bin_vector();
                }

                if (best_vb_gain < -F and best_vb_gain - curr_gain_padding < -F) {
                    std::cout << "close: encode_VByte" << std::endl;
                    encode_VByte();
                }

                if (curr_gain_padding > 0) {
                    std::cout << "close: curr_gain_padding > 0" << std::endl;
                    curr_bv_end = end;
                    curr_bv_end_index = n+1;
                    encode_bin_vector();
                } else {
                    std::cout << "close: curr_gain_padding <= 0" << std::endl;
                    curr_vb_end = end;
                    curr_vb_end_index = n+1;
                    encode_VByte();
                }
            };

            for (int64_t last_gain = 0; it != end; ++it,
                 last = curr, last_gain = curr_gain_padding)
            {
                curr = *it;
                curr_gain += VBBlock::posting_cost(curr, last)
                           - RBBlock::posting_cost(curr, last);
                std::cout << "it_index: " << it_index << std::endl;
                std::cout << "curr_gain: " << curr_gain << std::endl;
                curr_gain_padding = curr_gain;
                while(curr_gain_padding % 8 != 0){
                    curr_gain_padding -= 1;
                }
                std::cout << "curr_gain_padding: " << curr_gain_padding << std::endl;
                if (curr_vb_end_index < curr_bv_end_index) {
                    if (curr_vb_end_index-curr_begin_index % 4 != 0) {
                        std::cout << "curr_gain_padding: le resto 8" << std::endl;
                        curr_gain_padding-=8;
                    }
                }

                if (curr_bv_end_index < curr_vb_end_index) {
                    if (curr_bv_end_index-curr_begin_index % 4 != 0) {
                        std::cout << "curr_gain_padding: le resto 8" << std::endl;
                        curr_gain_padding-=8;
                    }
                }
                */
                /* 
                curr_begin_index  67
                curr_bv_end_index 31
                curr_vb_end_index 67
                /*
                it empezando en 1.
                17 doc_id igual a 1:
                1)  curr_gain -2  | best_vb_gain -8  | curr_vb_end 2  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -8
                2)  curr_gain -4  | best_vb_gain -8  | curr_vb_end 3  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -8
                3)  curr_gain -6  | best_vb_gain -8  | curr_vb_end 4  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -8
                4)  curr_gain -8  | best_vb_gain -8  | curr_vb_end 5  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -8
                5)  curr_gain -10 | best_vb_gain -16 | curr_vb_end 6  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -16
                6)  curr_gain -12 | best_vb_gain -16 | curr_vb_end 7  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -16
                7)  curr_gain -14 | best_vb_gain -16 | curr_vb_end 8  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -16
                8)  curr_gain -16 | best_vb_gain -16 | curr_vb_end 9  | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -16
                9)  curr_gain -18 | best_vb_gain -24 | curr_vb_end 10 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -24
                10) curr_gain -20 | best_vb_gain -24 | curr_vb_end 11 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -24
                11) curr_gain -22 | best_vb_gain -24 | curr_vb_end 12 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -24
                12) curr_gain -24 | best_vb_gain -24 | curr_vb_end 13 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -24
                13) curr_gain -26 | best_vb_gain -32 | curr_vb_end 14 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -32
                14) curr_gain -28 | best_vb_gain -32 | curr_vb_end 15 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -32
                15) curr_gain -30 | best_vb_gain -32 | curr_vb_end 16 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -32
                16) curr_gain -32 | best_vb_gain -32 | curr_vb_end 17 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -32
                17) curr_gain -34 | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 0 | curr_bv_end 1 | curr_gain_padding -40
                
                13 doc_id igual a 128:
                18) curr_gain -28 | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 0  | curr_bv_end 1  | curr_gain_padding -32
                19) curr_gain -22 | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 0  | curr_bv_end 1  | curr_gain_padding -24
                20) curr_gain -16 | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 0  | curr_bv_end 1  | curr_gain_padding -16
                21) curr_gain -10 | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 0  | curr_bv_end 1  | curr_gain_padding -16
                22) curr_gain -4  | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 0  | curr_bv_end 1  | curr_gain_padding -8
                23) curr_gain 2   | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 0  | curr_bv_end 1  | curr_gain_padding 0
                24) curr_gain 8   | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 8  | curr_bv_end 25 | curr_gain_padding 8
                25) curr_gain 14  | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 8  | curr_bv_end 25 | curr_gain_padding 0
                26) curr_gain 20  | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 8  | curr_bv_end 25 | curr_gain_padding 8
                27) curr_gain 26  | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 16 | curr_bv_end 28 | curr_gain_padding 16
                28) curr_gain 32  | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 24 | curr_bv_end 29 | curr_gain_padding 24
                29) curr_gain 38  | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 24 | curr_bv_end 29 | curr_gain_padding 24
                30) curr_gain 44  | best_vb_gain -40 | curr_vb_end 18 | best_bv_gain 32 | curr_bv_end 31 | curr_gain_padding 32

                encode_VByte: curr_gain 78  | best_vb_gain 0 | curr_vb_end 18 | best_bv_gain 72 | curr_bv_end 30 | curr_gain_padding 72

                masked: 26bytes-->208 bits
                varint: 13bytes+2*13bits-->130 bits
                varint: 13bytes+4bytes-->136 bits
                208-130=78
                208-136=72

                36 doc_id igual a 1:
                31) curr_gain 76  | best_vb_gain 0 | curr_vb_end 18 | best_bv_gain 72 | curr_bv_end 31 | curr_gain_padding 72
                32) curr_gain 74  | best_vb_gain 0 | curr_vb_end 18 | best_bv_gain 72 | curr_bv_end 31 | curr_gain_padding 72
                33) curr_gain 72  | best_vb_gain 0 | curr_vb_end 18 | best_bv_gain 72 | curr_bv_end 31 | curr_gain_padding 72
                34) curr_gain 70  | best_vb_gain 0 | curr_vb_end 18 | best_bv_gain 72 | curr_bv_end 31 | curr_gain_padding 64
                ...
                64) curr_gain 10  | best_vb_gain 0 | curr_vb_end 18 | best_bv_gain 72 | curr_bv_end 31 | curr_gain_padding 8
                65) curr_gain 8   | best_vb_gain 0 | curr_vb_end 18 | best_bv_gain 72 | curr_bv_end 31 | curr_gain_padding 8
                66) curr_gain 6   | best_vb_gain 0 | curr_vb_end 67 | best_bv_gain 72 | curr_bv_end 31 | curr_gain_padding 0

                encode_bin_vector: curr_gain -72 | best_vb_gain -72 | curr_vb_end 67 | best_bv_gain 0 | curr_bv_end 31 | curr_gain_padding -72

                masked: 36bytes--> 288 bits
                varint: 36bytes+2*36bits--> 360 bits
                varint: 36bytes+9bytes--> 360 bits
                288-360=-72

                encode_VByte: curr_gain 0 | best_vb_gain 0 | curr_vb_end 67 | best_bv_gain 0 | curr_bv_end 31 | curr_gain_padding 0

                -------------

                7 doc_id igual a 128:
                1) curr_gain 6  | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  0  | curr_bv_end 1 | curr_gain_padding 0
                2) curr_gain 12 | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  8  | curr_bv_end 3 | curr_gain_padding 8
                3) curr_gain 18 | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  16 | curr_bv_end 4 | curr_gain_padding 16
                4) curr_gain 24 | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  24 | curr_bv_end 5 | curr_gain_padding 24
                5) curr_gain 30 | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  24 | curr_bv_end 5 | curr_gain_padding 24
                6) curr_gain 36 | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  32 | curr_bv_end 7 | curr_gain_padding 32
                7) curr_gain 42 | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  40 | curr_bv_end 8 | curr_gain_padding 40

                33 doc_id igual a 1:
                8)  curr_gain 40 | best_vb_gain 0 | curr_vb_end 1 | best_bv_gain  40 | curr_bv_end 8 | curr_gain_padding 40
                ...
                28) curr_gain 0  | best_vb_gain -8 | curr_vb_end 29 | best_bv_gain  40 | curr_bv_end 8 | curr_gain_padding -8
                ...
                40) curr_gain -24| best_vb_gain -32 | curr_vb_end 41 | best_bv_gain  40 | curr_bv_end 8 | curr_gain_padding -32
            
                encode_bin_vector: curr_gain -66 | best_vb_gain -72 | curr_vb_end 41 | best_bv_gain 0 | curr_bv_end 8 | curr_gain_padding -72

                masked: 33bytes-->264 bits
                varint: 33bytes+2*33bits-->330 bits
                varint: 33bytes+9bytes-->336 bits
                264-330= -66
                264-336= -72

                ARRELGAR LOS ADD linea 72 y 96, Funcionan de pedo creo.

                ***Probar si el archivo final pesa menos con esto o simplemente sumando 2 al header.***
                */
/*
                std::cout << "curr_gain >= last_gain: " << curr_gain << " >= " << last_gain << std::endl;
                if (curr_gain >= last_gain) { // gain is not decreasing
                    std::cout << "curr_gain > best_bv_unfix_gain: " << curr_gain << " > " << best_bv_unfix_gain << std::endl;
                    if (curr_gain > best_bv_unfix_gain){
                        best_bv_unfix_gain = curr_gain;
                    }
                    std::cout << "curr_gain_padding > best_bv_gain: " << curr_gain_padding << " > " << best_bv_gain << std::endl;
                    if (curr_gain_padding > best_bv_gain) {
                        best_bv_gain = curr_gain_padding;
                        curr_bv_end = it + 1;
                        curr_bv_end_index = it_index + 1;
                    }
                    std::cout << "best_bv_gain: " << best_bv_gain << std::endl;
                    std::cout << "curr_bv_end: " << curr_bv_end << std::endl;
                    std::cout << "curr_bv_end_index: " << curr_bv_end_index << std::endl;
                    assert(best_vb_gain <= curr_gain_padding);
                    if (best_vb_gain - curr_gain_padding < -2 * F and best_vb_gain < -T) {
                        std::cout << "---encode_VByte---" << std::endl;
                        encode_VByte();
                    }

                } else { // gain is decreasing
                    std::cout << "curr_gain < best_vb_unfix_gain: " << curr_gain << " < " << best_vb_unfix_gain << std::endl;
                    if (curr_gain < best_vb_unfix_gain){
                        best_vb_unfix_gain = curr_gain;
                    }
                    std::cout << "curr_gain_padding <= best_vb_gain: " << curr_gain_padding << " <= " << best_vb_gain << std::endl;
                    if (curr_gain_padding <= best_vb_gain) {
                        best_vb_gain = curr_gain_padding;
                        curr_vb_end = it + 1;
                        curr_vb_end_index = it_index + 1;
                    }
                    std::cout << "best_vb_gain: " << best_vb_gain << std::endl;
                    std::cout << "curr_vb_end: " << curr_vb_end << std::endl;
                    std::cout << "curr_vb_end_index: " << curr_vb_end_index << std::endl;
                    assert(best_bv_gain >= curr_gain_padding);
                    if (best_bv_gain - curr_gain_padding > 2 * F and best_bv_gain > T) {
                        std::cout << "---encode_bin_vector---" << std::endl;
                        encode_bin_vector();
                    }
                }
                it_index++;
            }

            assert(it == end);
            close();

            return partition;
        }
    };
}
*/