#pragma once

#include "global_parameters_opt_vb.hpp"
#include "configuration.hpp"//#include "configuration_opt_vb.hpp"
#include "strict_sequence_opt_vb.hpp"
#include "util_opt_vb.hpp"

namespace pvb {

    template<typename BaseSequence = strict_sequence_opt_vb>
    struct positive_sequence_opt_vb {

        typedef BaseSequence base_sequence_type;
        typedef typename base_sequence_type::enumerator base_sequence_enumerator;

        template <typename Iterator>
        static void write(pisa::bit_vector_builder& bvb,
                          Iterator begin,
                          uint64_t universe, uint64_t n,
                          global_parameters_opt_vb const& params,
                        uint64_t& dense_short_freq, uint64_t& dense_medium_freq, uint64_t& dense_large_freq,
                        uint64_t& sparse_short_freq, uint64_t& sparse_medium_freq, uint64_t& sparse_large_freq,
                        uint64_t& dense_short_cost_freq, uint64_t& dense_medium_cost_freq, uint64_t& dense_large_cost_freq,
                        uint64_t& sparse_short_cost_freq, uint64_t& sparse_medium_cost_freq, uint64_t& sparse_large_cost_freq,
                        uint64_t& cantidad_integers_con_interpolative_freq,
                        uint64_t& cantidad_integers_con_varintg8iu_freq, bool dense_sparse)
        {

            assert(n > 0);
            std::vector<uint64_t> prefixes;
            prefixes.reserve(n);
            auto it = begin;
            prefixes.push_back(*it);
            ++it;
            for (uint64_t i = 1; i < n; ++i, ++it) {//Reverse dgaps for when dgaps is apply to freqs in "write" from "block_codecs"
                prefixes.push_back(prefixes.back() + *it);
            }
            base_sequence_type::write(bvb,
                                      prefixes.begin(),
                                      universe, n,
                                      params,
                                    dense_short_freq, dense_medium_freq, dense_large_freq,
                                    sparse_short_freq, sparse_medium_freq, sparse_large_freq,
                                    dense_short_cost_freq, dense_medium_cost_freq, dense_large_cost_freq,
                                    sparse_short_cost_freq, sparse_medium_cost_freq, sparse_large_cost_freq,
                                    cantidad_integers_con_interpolative_freq,
                                    cantidad_integers_con_varintg8iu_freq, dense_sparse);
        }

        static void decode(pisa::bit_vector const& bv,
                           uint32_t* out, uint64_t offset,
                           uint64_t universe, uint64_t n)
        {
            BaseSequence::decode(bv, out, offset, universe, n);
        }

        class enumerator {
        public:

            typedef std::pair<uint64_t, uint64_t> value_type; // (position, value)

            enumerator()
            {}

            enumerator(pisa::bit_vector const& bv, uint64_t offset,
                       uint64_t universe, uint64_t n,
                       global_parameters_opt_vb const& params, bool queries)
                : m_base_enum(bv, offset, universe, n, params, queries)
                , m_position(m_base_enum.size())
            {}

            value_type move(uint64_t position)
            {
                // we cache m_position and m_cur to avoid the call overhead in
                // the most common cases
                uint64_t prev = m_cur;
                if (position != m_position + 1) {
                    if (DS2I_UNLIKELY(position == 0)) {
                        // we need to special-case position 0
                        m_cur = m_base_enum.move(0).second;
                        m_position = 0;
                        return value_type(m_position, m_cur);
                    }
                    prev = m_base_enum.move(position - 1).second;
                }

                m_cur = m_base_enum.next().second;
                m_position = position;
                return value_type(position, m_cur - prev);
            }

            base_sequence_enumerator const& base() const {
                return m_base_enum;
            }

        private:
            base_sequence_enumerator m_base_enum;
            uint64_t m_position;
            uint64_t m_cur;
        };
    };
}
