#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <iomanip>
#include <locale>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <sys/time.h>
#include <sys/resource.h>

#include "../bit_vector.hpp"

namespace pvb {

    static const uint64_t alignment = 8;

    inline uint64_t push_pad_v2(pisa::bit_vector_builder& bvb, uint64_t alignment = 8) {
        uint64_t mod = bvb.size() % alignment;
        if (mod) {
            uint64_t pad = alignment - mod;
            bvb.append_bits(0, pad);
            assert(bvb.size() % alignment == 0);
            return pad;
        }
        return 0;
    }

    inline uint64_t eat_pad_v2(pisa::bit_vector::enumerator& it, uint64_t alignment = 8) {
        uint64_t mod = it.position() % alignment;
        if (mod) {
            uint64_t pad = alignment - mod;
            it.take(pad);
            assert(it.position() % alignment == 0);
            return pad;
        }
        return 0;
    }

}  // namespace pvb