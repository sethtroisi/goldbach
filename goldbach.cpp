// Copyright 2025 Seth Troisi
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <atomic>
#include <bitset>
#include <cassert>
#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <getopt.h>

#include <omp.h>
#include <primesieve.hpp>

using std::atomic;
using std::cout;
using std::endl;
using std::bitset;
using std::vector;
using namespace std::chrono;

void goldbach_search(const struct Config &config);
void goldbach_search_single(const Config &config);

const uint32_t THREAD_SPACING = 256'000;
const uint32_t OVERLAP = 1536;

struct Config {
    // start
    uint64_t n_0 = 0;
    // end
    uint64_t n_1 = 0;

    // K = number of primes
    int K = 200;

    int thread_i = -1;
    int threads = 1;
};

Config argparse(int argc, char* argv[]) {
    // NOTE: Remember to add to getopt_long(argc, argv, OPTIONS_STRING, ...) below
    static struct option long_options[] = {
        {"start",            required_argument, 0,  's' },
        {"end",              required_argument, 0,  'e' },
        {"threads",          required_argument, 0,  't' },
        {0,                  0,                 0,   0  }
    };

    Config config;

    int option_index = 0;
    char c;
    while ((c = getopt_long(argc, argv, "s:e:K:t:", long_options, &option_index)) >= 0) {
        switch (c) {
            case 's':
                config.n_0 = atoll(optarg);
                break;
            case 'e':
                config.n_1 = atoll(optarg);
                break;
            case 'K':
                config.K = atoi(optarg);
                break;
            case 't':
                config.threads = atoi(optarg);
                break;
            case 0: {
                int option_index = 0;
                printf("option %s arg %s\n", long_options[option_index].name, optarg);
                break;
            }

            default:
                printf("getopt returned \"%c\"(%d)\n", c, c);
                break;
        }
    }

    return config;
}

int main(int argc, char* argv[]) {
    // Display %'d with commas i.e. 12,345
    setlocale(LC_NUMERIC, "");


    Config config = argparse(argc, argv);

    goldbach_search(config);
}

atomic<uint16_t> max_needed_P(0);
bool set_max(uint16_t P_i) {
    // Use atomic::fetch_max after c++26

    uint16_t was = max_needed_P.load();
    while (was < P_i) {
        // Try to store
        if (max_needed_P.compare_exchange_strong(was, P_i)) {
            return true;
        }
        was = max_needed_P.load();
    }
    return false;
}

int brute_check(uint64_t N, const vector<uint16_t> &P) {
    primesieve::iterator it(N);
    size_t P_i = 0;
    uint64_t q = P[P_i];
    for (uint64_t p = it.prev_prime(); P_i < P.size(); p = it.prev_prime()) {
        while (p + q < N && (++P_i) < P.size()) {
            q = P[P_i];
        }
        // Check if p + q = N
        if (p + q == N) {
            if (set_max(P_i)) {
                printf("\t\t%lu = %lu + %lu (%lu)\n", N, p, q, P_i);
            }
            return P_i;
        }
    }
    return -1;
}

const size_t MASK_BITS = 64;
void generate_P_masks(vector<uint16_t> P, vector<uint64_t> *maskss) {
    for (size_t i = 0; i < MASK_BITS; i++) {
        auto &masks = maskss[i];

        uint64_t mask = 0;
        uint32_t offset = 0;
        for (const auto p : P) {
            uint32_t bit = (i - offset) + (p+1) / 2;
            if (bit >= MASK_BITS) {
                masks.push_back(mask);
                mask = 0;
                offset += MASK_BITS;
                bit -= MASK_BITS;
                assert( bit < MASK_BITS );
                mask = 0;
            }
            mask |= 1UL << bit;
        }
        if (mask > 0) {
            masks.push_back(mask);
        }

        if (false) {
            printf("mask(%lu) =", i);
            for (auto mask : masks) {
                printf(" %16lx", mask);
            }
            printf("\n");
        }
    }
}

void goldbach_search(const Config &config) {
    if (config.threads == 1) {
        goldbach_search_single(config);
    } else {
        // Not the best split but what ever
        uint64_t spaces = (config.n_1 - config.n_0) / THREAD_SPACING + 1;

        #pragma omp parallel for num_threads(config.threads)
        for (int t = 0; t < config.threads; t++) {
            Config config_t(config);
            config_t.n_0 = config.n_0 + (t * spaces / config.threads) * THREAD_SPACING;
            config_t.n_1 = config.n_0 + ((t+1) * spaces / config.threads) * THREAD_SPACING;
            config_t.thread_i = t;
            printf("\t[%lu, %lu]\n", config_t.n_0, config_t.n_1);
            goldbach_search_single(config_t);
        }
    }
}

void goldbach_search_single(const Config &config) {
    vector<uint16_t> P;
    primesieve::generate_primes(3, config.K, &P);

    // 64 copies of P mask depending on what bit we start at
    vector<uint64_t> P_masks[MASK_BITS];
    generate_P_masks(P, P_masks);

    // See https://oeis.org/A025019 and
    // "EMPIRICAL VERIFICATION OF THE EVEN GOLDBACH CONJECTURE" by TOMAS OLIVEIRA e SILVA
    vector<uint16_t> BRUTE_P;
    primesieve::generate_primes(3, 10000, &BRUTE_P);

    printf("Searching [%lu, %lu] P <= %u, spacing = %u (overlap: %u)\n",
            config.n_0, config.n_1, config.K, THREAD_SPACING, OVERLAP);
    printf("\t|P| = %lu (%lu masks)\n", P.size(), P_masks[31].size());
    printf("\n");
    assert( (unsigned) config.K <= OVERLAP );
    assert( config.n_1 >= config.n_0 );

    assert( THREAD_SPACING % 64 == 0);
    assert( OVERLAP % 64 == 0);

    size_t THREAD_BYTES = THREAD_SPACING / 8;
    size_t THREAD_UINT64_PER = THREAD_BYTES / 8;
    size_t OVERLAP_BYTES = OVERLAP / (2 * 8);

    uint64_t *cleared = (uint64_t*) calloc(THREAD_BYTES + OVERLAP_BYTES, 1);
    uint64_t *overlap = cleared + THREAD_UINT64_PER;
    uint64_t *cleared_after_overlap = cleared + (OVERLAP_BYTES / 8);

    size_t NON_OVERLAPPED_BYTES = THREAD_BYTES - OVERLAP_BYTES;
    uint64_t FULL_MASK = std::numeric_limits<uint64_t>::max();

    // Break up [n_0, n_1] into ranges of THREAD_SPACING
    uint64_t M = 2 * THREAD_SPACING;
    uint64_t M_COUNT = (config.n_1 - config.n_0 + 1 + M-1) / M;
    atomic<size_t> brute_searches = 0;

    // Handle overflow from previous interval
    if (config.n_0 > 0) {
        uint64_t n_0 = config.n_0;
        assert( n_0 >= OVERLAP );
        primesieve::iterator it(std::max(n_0 - OVERLAP, uint64_t(3)));
        for (uint64_t prime = it.next_prime(); prime < n_0; prime = it.next_prime()) {
            for (const uint32_t q : P) {
                int64_t i = prime + q - n_0;
                if (i >= 0) {
                    uint64_t t = i >> 1;
                    cleared[t >> 6] |= 1UL << (t & (MASK_BITS-1));
                    //printf("underflow %lu = %lu + %u -> (%lu, %lu)\n",
                    //       prime+q, prime, q, t>>6, t & (MASK_BITS-1));
                }
            }
        }
    }

    for (uint64_t m_i = 0 ; m_i < M_COUNT; m_i ++) {
        uint64_t M_0 = config.n_0 + M * m_i;
        uint64_t M_1 = M_0 + M-1;
        // [M_0, M_1]

        primesieve::iterator it(std::max(M_0, uint64_t(3)));
        uint64_t prime;

        for (prime = it.next_prime(); prime < M_1; prime = it.next_prime()) {
            uint64_t i = prime - M_0;
            assert( i <= M );
            if (false) {
                for (const uint32_t q : P) {
                    uint32_t t = (i + q) >> 1;
                    cleared[t >> 6] |= 1UL << (t & (MASK_BITS-1));
                }
            } else {
                // Faster mask strategy
                uint32_t t = i >> 1;
                uint32_t mask_index = t & (MASK_BITS-1);
                uint32_t c_i = t >> 6;
                for (auto mask : P_masks[mask_index]) {
                    cleared[c_i++] |= mask;
                }

                // Verify mask strategy
                if (false) {
                    for (const uint32_t q : P) {
                        uint32_t t = (i + q) >> 1;
                        bool valid = cleared[t >> 6] & (1UL << (t & (MASK_BITS-1)));

                        if ( !valid ) {
                            printf("%u -> (%u, %u) from %lu + %u\n", t, t >> 6, t & 64, prime, q);

                            uint32_t s = i >> 1;
                            uint32_t c_i = s >> 6;
                            printf("maskss[%u]\n", mask_index);
                            for (auto mask : P_masks[mask_index]) {
                                printf("\t%u |= %lx ([%u] = %lu)\n", c_i++, mask, t & 63, mask & (1UL << (t & 63)));
                            }
                        }
                        assert( valid );
                    }
                }
            }
        }

        // Verify all cleared
        for (size_t i = 0; i < THREAD_UINT64_PER; i++) {
            uint64_t bits = cleared[i];
            if (bits != FULL_MASK) {
                for (size_t j = 0; j < MASK_BITS; j++) {
                    if (!(bits & (1UL << j))) {
                        uint64_t t = M_0 + 2*(MASK_BITS*i+j);
                        brute_searches++;
                        if (t > 4) {
                            int q_i = brute_check(t, BRUTE_P);
                            // Need to ignore
                            if (q_i < 0 || q_i < (signed) P.size() ) {
                                printf("\t%lu not found [%lu, %lu) = %lx\n",
                                        t, M_0 + 2*MASK_BITS*i, M_0 + 2*MASK_BITS*(i+1), bits);
                            }
                            assert( q_i >= 0 && q_i >= (signed) P.size() );
                        }
                    }
                }
            }
        }

        // Move bits down and clear stuff out
        memcpy(cleared, overlap, OVERLAP_BYTES);
        memset(cleared_after_overlap, 0, NON_OVERLAPPED_BYTES);

        if (m_i < 10 || ((m_i & 511) == 0) || (m_i + 5 > M_COUNT)) {
            printf("\t[%lu, %lu] (brute: %lu)\n", M_0, M_1, brute_searches.load());
        }
    }
}
