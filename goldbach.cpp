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

using std::cout;
using std::endl;
using std::bitset;
//using std::mutex;
using std::vector;
using namespace std::chrono;

void goldbach_search(const struct Config& config);

const uint32_t THREAD_SPACING = 512'000;
const uint32_t OVERLAP = 512;

struct Config {
    // start
    uint64_t n_0 = 0;
    // end
    uint64_t n_1 = 0;

    // K = number of primes
    int K = 200;

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

void goldbach_search(const Config &config) {

    vector<uint16_t> P;
    primesieve::generate_primes(3, config.K, &P);

    printf("Searching [%lu, %lu] P <= %u, spacing = %u (overlap: %u)\n",
            config.n_0, config.n_1, config.K, THREAD_SPACING, OVERLAP);
    printf("\t|P| = %lu\n", P.size());
    assert( (unsigned) config.K <= OVERLAP );

    assert( THREAD_SPACING % 32 == 0);
    assert( OVERLAP % 64 == 0);

    size_t THREAD_BYTES = THREAD_SPACING / 8;
    size_t THREAD_UINT32_PER = THREAD_BYTES / 4;
    size_t OVERLAP_BYTES = OVERLAP / (2 * 8);

    uint32_t *cleared = (uint32_t*) malloc(THREAD_BYTES + OVERLAP_BYTES);
    uint32_t *overlap = cleared + (THREAD_BYTES / 4);
    uint32_t *cleared_after_overlap = cleared + (OVERLAP_BYTES / 4);

    size_t NON_OVERLAPPED_BYTES = THREAD_BYTES - OVERLAP_BYTES;
    uint32_t FULL_UINT32 = std::numeric_limits<uint32_t>::max();

    // #pragma omp parallel for num_threads(THREADS)

    // Break up [n_0, n_1] into ranges of THREAD_SPACING
    uint64_t M = 2 * THREAD_SPACING;
    uint64_t M_COUNT = (config.n_1 - config.n_0 + 1 + M-1) / M;

    for (uint64_t m_i = 0 ; m_i < M_COUNT; m_i ++) {
        uint64_t M_0 = config.n_0 + M * m_i;
        uint64_t M_1 = M_0 + M-1;
        // [M_0, M_1]
        printf("\t[%lu, %lu]\n", M_0, M_1);

        primesieve::iterator it(M_0);
        uint64_t prime;
        for (prime = it.next_prime(); prime < M_1; prime = it.next_prime()) {
            uint64_t i = prime - M_0;
            assert( i <= M );
            // TODO this can be converted to masks that get ORd with cleared.
            for (const uint32_t q : P) {
                uint32_t t = (i + q) >> 1;
                cleared[t >> 5] |= 1 << (t & 31);
            }
        }

        // Verify all cleared
        for (size_t i = 0; i < THREAD_UINT32_PER; i++) {
            uint32_t bits = cleared[i];
            if (bits != FULL_UINT32) {
                // Ignores 2 and 4.
                if (M_0 == 0 and i == 0 && bits == (FULL_UINT32 - 3))
                    continue;
                for (size_t j = 0; j < 32; j++) {
                    if (!(bits & (1 << j))) {
                        printf("\t%lu not found [%lu, %lu) = %x\n",
                                M_0 + 2*(32*i+j), M_0 + 2*32*i, M_0 + 2*32*(i+1), bits);
                    }
                }
            }
        }

        // Move bits down and clear stuff out
        memcpy(cleared, overlap, OVERLAP_BYTES);
        memset(cleared_after_overlap, 0, NON_OVERLAPPED_BYTES);
    }
}
