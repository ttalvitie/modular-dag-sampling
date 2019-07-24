#include "lognum.h"

static std::array<double, 16384> init_log_factorial_table() {
    std::array<double, 16384> ret;
    for(int i = 0; i < 16384; ++i) {
        ret[i] = std::lgamma((double)i + 1.0);
    }
    return ret;
}

const std::array<double, 16384> Lognum::log_factorial_table = init_log_factorial_table();