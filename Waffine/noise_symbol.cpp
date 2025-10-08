#include "noise_symbol.h"
//
// Created by will on 9/24/25.
//

// Prime noise symbols with value 0.
static noise_symbol_t max_noise_symbol = 0;

noise_symbol_t new_noise_symbol() {
    return max_noise_symbol++;
}

bool valid_noise_symbol(noise_symbol_t symbol) {
    return symbol < max_noise_symbol;
}
