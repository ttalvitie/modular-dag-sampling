#pragma once

#include "common.h"

#include <immintrin.h> // _pext_u32

template <typename T>
class SubTable {
public:
    SubTable(uint32_t n) {
        data.resize(1 << n);
        for(uint32_t U = 0; U < data.size(); ++U) {
            data[U].resize(1 << __builtin_popcount(U));
        }
    }
    SubTable() : SubTable(0) {}

    T& operator()(uint32_t R, uint32_t U) {
        return data[U][_pext_u32(R, U)];
    }
    const T& operator()(uint32_t R, uint32_t U) const {
        return data[U][_pext_u32(R, U)];
    }

private:
    std::vector<std::vector<T>> data;
};
