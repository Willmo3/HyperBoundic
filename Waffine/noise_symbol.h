//
// Created by will on 9/24/25.
//

#ifndef WAFFINE_NOISESYMBOL_H
#define WAFFINE_NOISESYMBOL_H
#include <cstdint>

/*
 * Note: we use 32-bit integers here to avoid class overhead, specifically with hashing.
 */

/**
 * Noise symbols are integers that refer to unique sources of error in affine forms.
 */
typedef uint32_t noise_symbol_t;

/**
 * @return A new unique noise symbol.
 */
noise_symbol_t new_noise_symbol();

/**
 * @return whether the provided symbol is <= the last symbol allocated.
 */
bool valid_noise_symbol(noise_symbol_t symbol);

#endif //WAFFINE_NOISESYMBOL_H