#pragma once

#include "boost/preprocessor/cat.hpp"
#include "boost/preprocessor/seq/for_each.hpp"
#include "boost/preprocessor/stringize.hpp"

#include "codec/block_codecs.hpp"
#include "codec/maskedvbyte.hpp"
#include "codec/qmx.hpp"
#include "codec/simdbp.hpp"
#include "codec/simple16.hpp"
#include "codec/simple8b.hpp"
#include "codec/streamvbyte.hpp"
#include "codec/varintgb.hpp"

#include "binary_freq_collection.hpp"
#include "block_freq_index.hpp"

#include "freq_index.hpp"
#include "sequence/partitioned_sequence.hpp"
#include "sequence/positive_sequence.hpp"
#include "sequence/uniform_partitioned_sequence.hpp"

//opt_vb
#include "opt_vb/block_codecs_opt_vb.hpp"
#include "opt_vb/freq_index_opt_vb.hpp"
#include "opt_vb/partitioned_vb_sequence_opt_vb.hpp"
#include "opt_vb/positive_sequence_opt_vb.hpp"
#include "opt_vb/partitioned_sequence_enumerator_opt_vb.hpp"
#include "opt_vb/indexed_sequence_opt_vb.hpp"
#include "opt_vb/compact_ranked_bitvector_opt_vb.hpp"

namespace pisa {
using ef_index = freq_index<compact_elias_fano, positive_sequence<strict_elias_fano>>;

using single_index = freq_index<indexed_sequence, positive_sequence<>>;

using pefuniform_index =
    freq_index<uniform_partitioned_sequence<>, positive_sequence<uniform_partitioned_sequence<strict_sequence>>>;

using pefopt_index =
    freq_index<partitioned_sequence<>, positive_sequence<partitioned_sequence<strict_sequence>>>;

using block_optpfor_index = block_freq_index<pisa::optpfor_block>;
using block_varintg8iu_index = block_freq_index<pisa::varint_G8IU_block>;
using block_streamvbyte_index = block_freq_index<pisa::streamvbyte_block>;
using block_maskedvbyte_index = block_freq_index<pisa::maskedvbyte_block>;
using block_varintgb_index = block_freq_index<pisa::varintgb_block>;
using block_interpolative_index = block_freq_index<pisa::interpolative_block>;
using block_qmx_index = block_freq_index<pisa::qmx_block>;
using block_simple8b_index = block_freq_index<pisa::simple8b_block>;
using block_simple16_index = block_freq_index<pisa::simple16_block>;
using block_simdbp_index = block_freq_index<pisa::simdbp_block>;

using masked_bitvector_index =
    freq_index_opt_vb<
        pvb::partitioned_vb_sequence_opt_vb<
            pvb::maskedvbyte_block_opt_vb, // Encoder1
            pvb::compact_ranked_bitvector_opt_vb, // Encoder2
            pvb::partitioned_sequence_enumerator_opt_vb<
                pvb::indexed_sequence_opt_vb<
                    pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                    pvb::compact_ranked_bitvector_opt_vb // Encoder2 Enumerator
                >
            >
        >,
        pvb::positive_sequence_opt_vb<
            pvb::partitioned_vb_sequence_opt_vb<
                pvb::maskedvbyte_block_opt_vb, // Encoder1
                pvb::compact_ranked_bitvector_opt_vb, // Encoder2
                pvb::partitioned_sequence_enumerator_opt_vb<
                    pvb::indexed_sequence_opt_vb<
                        pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                        pvb::compact_ranked_bitvector_opt_vb // Encoder2 Enumerator
                    >
                >
            >
        >
    >;

using masked_vgb_index =
    freq_index_opt_vb<
        pvb::partitioned_vb_sequence_opt_vb<
            pvb::maskedvbyte_block_opt_vb, // Encoder1
            pvb::varintgb_block, // Encoder2
            pvb::partitioned_sequence_enumerator_opt_vb<
                pvb::indexed_sequence_opt_vb<
                    pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                    pvb::block_sequence_opt_vb<pvb::varintgb_block> // Encoder2 Enumerator
                >
            >
        >,
        pvb::positive_sequence_opt_vb<
            pvb::partitioned_vb_sequence_opt_vb<
                pvb::maskedvbyte_block_opt_vb, // Encoder1
                pvb::varintgb_block, // Encoder2
                pvb::partitioned_sequence_enumerator_opt_vb<
                    pvb::indexed_sequence_opt_vb<
                        pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                        pvb::block_sequence_opt_vb<pvb::varintgb_block> // Encoder2 Enumerator
                    >
                >
            >
        >
    >;

using masked_vg8iu_index =
    freq_index_opt_vb<
        pvb::partitioned_vb_sequence_opt_vb<
            pvb::maskedvbyte_block_opt_vb, // Encoder1
            pvb::varintg8iu_block, // Encoder2
            pvb::partitioned_sequence_enumerator_opt_vb<
                pvb::indexed_sequence_opt_vb<
                    pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                    pvb::block_sequence_opt_vb<pvb::varintg8iu_block> // Encoder2 Enumerator
                >
            >
        >,
        pvb::positive_sequence_opt_vb<
            pvb::partitioned_vb_sequence_opt_vb<
                pvb::maskedvbyte_block_opt_vb, // Encoder1
                pvb::varintg8iu_block, // Encoder2
                pvb::partitioned_sequence_enumerator_opt_vb<
                    pvb::indexed_sequence_opt_vb<
                        pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                        pvb::block_sequence_opt_vb<pvb::varintg8iu_block> // Encoder2 Enumerator
                    >
                >
            >
        >
    >;

using masked_stream_index =
    freq_index_opt_vb<
        pvb::partitioned_vb_sequence_opt_vb<
            pvb::maskedvbyte_block_opt_vb, // Encoder1
            pvb::streamvbyte_block, // Encoder2
            pvb::partitioned_sequence_enumerator_opt_vb<
                pvb::indexed_sequence_opt_vb<
                    pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                    pvb::block_sequence_opt_vb<pvb::streamvbyte_block> // Encoder2 Enumerator
                >
            >
        >,
        pvb::positive_sequence_opt_vb<
            pvb::partitioned_vb_sequence_opt_vb<
                pvb::maskedvbyte_block_opt_vb, // Encoder1
                pvb::streamvbyte_block, // Encoder2
                pvb::partitioned_sequence_enumerator_opt_vb<
                    pvb::indexed_sequence_opt_vb<
                        pvb::block_sequence_opt_vb<pvb::maskedvbyte_block_opt_vb>, // Encoder1 Enumerator
                        pvb::block_sequence_opt_vb<pvb::streamvbyte_block> // Encoder2 Enumerator
                    >
                >
            >
        >
    >;

}  // namespace pisa

#define PISA_INDEX_TYPES                                                                    \
    (ef)(single)(pefuniform)(pefopt)(block_optpfor)(block_varintg8iu)(block_streamvbyte)(   \
        block_maskedvbyte)(block_interpolative)(block_qmx)(block_varintgb)(block_simple8b)( \
        block_simple16)(block_simdbp)(masked_bitvector)(masked_vgb)(masked_vg8iu)(masked_stream)
 #define PISA_BLOCK_INDEX_TYPES                                                                    \
    (block_optpfor)(block_varintg8iu)(block_streamvbyte)(block_maskedvbyte)(block_interpolative)( \
        block_qmx)(block_varintgb)(block_simple8b)(block_simple16)(block_simdbp)