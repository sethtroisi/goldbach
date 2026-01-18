// Copyright 2026 Seth Troisi
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
#include <cassert>
#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <getopt.h>

#include <primesieve.hpp>

using std::array;
using std::atomic;
using std::vector;
using namespace std::chrono;

void goldbach_gpu_search(const struct Config &config);

// How much each bytes GPU thread processes
constexpr uint32_t GPU_BYTES = 320'000;

// See https://oeis.org/A025019 and
// "EMPIRICAL VERIFICATION OF THE EVEN GOLDBACH CONJECTURE" by TOMAS OLIVEIRA e SILVA
constexpr uint32_t MAX_Q = 10'000;


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

    goldbach_gpu_search(config);
}


// Depends on K and must be handled at compile time.
const size_t MASK_LENGTH = (64 + 1024-1) / (2*64);
const size_t MASK_BITS = 64;
void generate_P_masks(vector<uint16_t> P, uint64_t maskss[MASK_BITS][MASK_LENGTH]) {
    for (size_t i = 0; i < MASK_BITS; i++) {
        for (size_t j = 0; j < MASK_LENGTH; j++) {
            maskss[i][j] = 0;
        }

        auto &masks = maskss[i];

        uint64_t mask = 0;
        uint64_t m_i = 0;
        uint32_t offset = 0;
        for (const auto p : P) {
            uint32_t bit = (i - offset) + (p+1) / 2;
            if (bit >= MASK_BITS) {
                masks[m_i++] = mask;
                mask = 0;
                offset += MASK_BITS;
                bit -= MASK_BITS;
                assert( bit < MASK_BITS );
            }
            mask |= 1UL << bit;
        }
        if (mask > 0) {
            masks[m_i] = 0;
        }

        if (false) {
            printf("mask(%lu) =", i);
            for (size_t j = 0; j < MASK_LENGTH; j++) {
                printf(" %16lx", maskss[i][j]);
            }
            printf("\n");
        }
    }
}

/*
void goldbach_search_single(const Config &config) {
    vector<uint16_t> P;
    primesieve::generate_primes(3, config.K, &P);

    // 64 copies of P mask depending on what bit we start at
    // Indexed as (p & 63) *
    MASK_MAX = (MASK_BITS - 1) + (P.back() + 1) / 2;
    MASK_MAX_LENGTH = (MASK_MAX + MASK_BITS-1) / MASK_BITS;
    uint64_t P_masks[MASK_BITS][MASK_MAX_LENGTH];
    generate_P_masks(P, P_masks);

    // See https://oeis.org/A025019 and
    // "EMPIRICAL VERIFICATION OF THE EVEN GOLDBACH CONJECTURE" by TOMAS OLIVEIRA e SILVA
    vector<uint16_t> Q;
    primesieve::generate_primes(3, MAX_Q, &Q);

    uint64_t prime_start = 0;
    uint64_t prime_end = 0;

    uint64_t FULL_MASK = std::numeric_limits<uint64_t>::max();

    // Break up [n_0, n_1] into ranges of GPU_BYTES
    uint64_t M = 2 * THREAD_SPACING;
    uint64_t M_COUNT = (n_1 - n_0 + 1 + M-1) / M;
    size_t brute_searches = 0;
    double setup_time = 0.0;
    double search_time = 0.0;
    double brute_time = 0.0;

    primesieve::iterator it(n_0);

    uint64_t prime = it.next_prime();
    for (uint64_t m_i = 0 ; m_i < M_COUNT; m_i ++) {
        uint64_t M_0 = n_0 + M * m_i;
        uint64_t M_1 = M_0 + M-1 + OVERLAP_AFTER;
        // Process all numbers
        // [M_0, M_1 + OVERLAP_AFTER)
        // First interval verifies [0, M_1 + OVERLAP)
        // Other intervals verifiy [M_0 + OVERLAP, M_1 + OVERLAP]

        // Handle updating prime_start to last prime <= M_0 - MAX_Q
        if (m_i > 0) {
            auto setup_t = high_resolution_clock::now();

            min_needed_prime = std::max<int64_t>((int64_t) M_0 - MAX_Q, int64_t(3));
            //printf("\tAdvancing up to %lu, prime_start[%lu] = %lu to ", min_needed_prime, prime_start, circular_primes[prime_start]);
            uint32_t next_i = prime_start;
            while (circular_primes[next_i] <= min_needed_prime) {
                assert( prime_start != prime_end ); // Can't advance past prime_end
                prime_start = next_i;
                next_i = next_i == (PRIMES_SIZE - 1) ? 0 : next_i + 1;
            }
            //printf("[%lu] = %lu\n", prime_start, circular_primes[prime_start]);
            assert( circular_primes[prime_start] <= min_needed_prime );
            assert( circular_primes[prime_start] + MAX_PRIME_GAP > min_needed_prime );
            assert( circular_primes[next_i] > min_needed_prime );

            setup_time += duration<double>(high_resolution_clock::now() - setup_t).count();
        }

        // Should/Could check prev_prime
        assert(M_0 <= prime && prime < M_0 + MAX_PRIME_GAP);
        auto search_t = high_resolution_clock::now();

        for (; prime <= M_1; prime = it.next_prime()) {
            circular_primes[prime_end++] = prime;
            if (prime_end == PRIMES_SIZE) prime_end = 0;

            uint64_t i = prime - M_0;
            assert( i <= M );
            if (true) {
                // Fast mask strategy
                uint32_t t = i >> 1;
                uint32_t mask_index = t & (MASK_BITS-1);
                uint32_t c_i = i >> 7; // t >> 6;
                for (auto mask : P_masks[mask_index]) {
                    cleared[c_i++] |= mask;
                }
            } else {
                // Slower prime-by-prime method of masking of bits.
                for (const uint32_t q : P) {
                    uint32_t t = (i + q) >> 1;
                    cleared[t >> 6] |= 1UL << (t & (MASK_BITS-1));
                }
            }
        }

        search_time += duration<double>(high_resolution_clock::now() - search_t).count();

        // Verify all cleared
        {
            // Add final prime at prime_end without increment.
            circular_primes[prime_end] = prime;

            auto brute_t = high_resolution_clock::now();

            //printf("Verifying [%lu, %lu] primes(indexes %lu %lu) [%lu, %lu]\n",
            //       M_0, M_1, prime_start, prime_end, circular_primes[prime_start], circular_primes[prime_end]);

            assert( (circular_primes[prime_start] == 3) || (circular_primes[prime_start] + MAX_Q < M_0) );
            assert( circular_primes[prime_end] >= M_1 );

            // Search forward till prime just less than M_0.
            uint32_t prime_i = prime_start;
            uint32_t next_i = prime_start;
            while (circular_primes[next_i] <= M_0) {
                assert( prime_i != prime_end ); // Can't advance past prime_end
                prime_i = next_i;
                next_i = next_i == (PRIMES_SIZE - 1) ? 0 : next_i + 1;
            }

            for (size_t i = 0; i < THREAD_UINT64_PER; i++) {
                uint64_t bits = cleared[i];
                if (bits != FULL_MASK) {
                    for (size_t j = 0; j < MASK_BITS; j++) {
                        if (!(bits & (1UL << j))) {
                            brute_searches++;
                            uint64_t t = M_0 + 2*(MASK_BITS*i+j);
                            if (t > 4) {
                                // Could be (t - K) but that verifies less
                                while (circular_primes[prime_i] < t) {
                                    assert( prime_i != prime_end ); // Can't advance past prime_end
                                    prime_i = next_i;
                                    next_i = next_i == (PRIMES_SIZE - 1) ? 0 : next_i + 1;
                                }

                                size_t q_i;
                                {
                                    // Walk backwards in circular_primes and forwards in Q
                                    uint64_t p_i = prime_i;
                                    uint64_t p = circular_primes[p_i];
                                    uint64_t p_i_0 = prime_i;
                                    uint64_t p_0 = p;
                                    assert( p >= t );

                                    bool valid = false;
                                    q_i = 0;
                                    for (uint64_t q = Q[q_i]; q_i < Q.size(); q = Q[++q_i]) {
                                        while (p + q > t) {
                                            // Can't advance backwards to less than prime_start
                                            assert( p_i != prime_start );
                                            p_i = (p_i == 0) ? (PRIMES_SIZE - 1) : (p_i - 1);
                                            p = circular_primes[p_i];
                                        }
                                        // Check if p + q = N = t
                                        if (p + q == t) {
                                            if (set_max(q_i)) {
                                                printf("\t\t%lu = %lu + %lu (%lu)\n", t, p, q, q_i);
                                            }
                                            valid = true;
                                            break;
                                        }
                                    }
                                    if (!valid) {
                                        printf("Didn't find p + q = %lu (searched p = %lu to %lu (indexes %lu to %lu), q = %u to %u)\n",
                                                t, p_0, p, p_i_0, p_i, Q.front(), Q.back());
                                        assert( false );
                                    }
                                }

                                if (q_i < P.size() ) {
                                    printf("\t%lu not found [%lu, %lu) = %lx\n",
                                            t, M_0 + 2*MASK_BITS*i, M_0 + 2*MASK_BITS*(i+1), bits);
                                    assert( false );
                                }
                                assert( Q[q_i] > config.K );
                            }
                        }
                    }
                }
            }
            brute_time += duration<double>(high_resolution_clock::now() - brute_t).count();
        }


        // Move bits down and clear stuff out
        memcpy(cleared, overlap, OVERLAP_BYTES);
        memset(cleared_after_overlap, 0, NON_OVERLAPPED_BYTES);

        if (m_i < 10 || (m_i + 5 > M_COUNT)
                || (m_i < 1024 && ((m_i & 127) == 0))
                || (m_i < 10240 && ((m_i & 2047) == 0))
                || ((m_i & 0xFFFFF) == 0)) {
            printf("\t[%lu, %lu] (%.1f%%) (brute: %lu)\n",
                   M_0, M_1, 100.0 * m_i / M_COUNT, brute_searches);
        }
    }

    printf("\t[%lu, %lu) took %.1f = setup %.1f, search %.1f, verify %.1f\n",
            n_0, n_1, setup_time + search_time + brute_time,
            setup_time, search_time, brute_time);
}
*/

#ifndef __CUDACC__
#error "This file should only be compiled with nvcc"
#endif

#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <cuda.h>
#include "CUDASieve/cudasieve.hpp"


#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)

inline void cuda_check(cudaError_t status, const char *action=NULL, const char *file=NULL, int32_t line=0) {
  if (status != cudaSuccess) {
    fprintf (stderr, "CUDA error (%d) occurred: %s\n", status, cudaGetErrorString(status));
    if (action!=NULL)
      fprintf (stderr, "While running %s   (file %s, line %d)\n", action, file, line);
    exit(EXIT_FAILURE);
  }
}


__global__ void kernel_goldbach_search(
        uint64_t found_start,
        uint64_t found_end,
        uint64_t *found,
        uint64_t num_instances,
        uint64_t *primes,
        uint64_t *q_primes,
        uint64_t masks[MASK_BITS][MASK_LENGTH]
        ) {

  // decode an instance_i number from the blockIdx and threadIdx
  //int32_t instance_i = threadIdx.x; //blockIdx.x*blockDim.x + threadIdx.x;


  // PROCESS ALL PRIMESS
  // XXX: PASS OVERLAP FLAG

  // CHECK FOR GAPS AND DO THE THING (HERE OR CPU?)
  // XXX: return max_needed_q

  return;

}

void goldbach_gpu_search(const Config &config)
{
    assert( config.n_0 % 2 == 0 );

    vector<uint16_t> P;
    primesieve::generate_primes(3, config.K, &P);

    // 64 copies of P mask depending on what bit we start at
    uint64_t test = (MASK_BITS + P.back()) / (2 * MASK_BITS);
    printf("%lu vs %lu\n", test, MASK_LENGTH);
    assert( test == MASK_LENGTH );
    uint64_t P_masks[MASK_BITS][MASK_LENGTH];
    generate_P_masks(P, P_masks);

    // See https://oeis.org/A025019 and
    // "EMPIRICAL VERIFICATION OF THE EVEN GOLDBACH CONJECTURE" by TOMAS OLIVEIRA e SILVA
    vector<uint16_t> Q;
    primesieve::generate_primes(3, MAX_Q, &Q);

    //cudaEvent_t global_start, batch_start, stop;
    //CUDA_CHECK(cudaEventCreate (&global_start));
    //CUDA_CHECK(cudaEventCreate (&batch_start));
    //CUDA_CHECK(cudaEventCreate (&stop));
    //CUDA_CHECK(cudaEventRecord (global_start));

    // Handle this many more at the end
    uint64_t OVERLAP = ((config.K + 63) / 64) * 64;
    uint64_t OVERLAP_BYTES = OVERLAP / 2 / 8;
    assert( OVERLAP >= config.K );
    assert( (unsigned) config.K <= OVERLAP );

    /**
     * GPU handles overlap by processing the first few p
     * then marking the bottom of it's range as "safe" to modify by other threads
     */
    size_t INTERVAL_BYTES = GPU_BYTES * config.threads;
    size_t TOTAL_BYTES = INTERVAL_BYTES + OVERLAP_BYTES;
    assert( GPU_BYTES % 64 == 0 );
    assert( TOTAL_BYTES % 64 == 0 );
    assert( INTERVAL_BYTES > OVERLAP_BYTES );

    const auto n_0 = config.n_0 > OVERLAP ? (config.n_0 - OVERLAP) : 0;
    const auto n_1 = config.n_1;

    printf("Testing [%lu, %lu) P <= %u, GPU_BYTES = %u, %lu primes\n\n",
            n_0, n_1, config.K, GPU_BYTES, P.size());
    assert( n_1 >= n_0 );

    //uint64_t *gpu_cleared;
    //CUDA_CHECK(cudaMalloc((void **)&gpu_cleared, TOTAL_BYTES));
    //CUDA_CHECK(cudaMemset((void **)&gpu_cleared, 0, TOTAL_BYTES));

    //uint64_t *gpu_primes;
    //CUDA_CHECK(cudaMalloc((void **)&gpu_primes, PRIME_BYTES);

    size_t interval_size = 2 * 8 * INTERVAL_BYTES;
    uint64_t interval_start = n_0;
    uint64_t interval_end = std::min(interval_start + interval_size, n_1);

    uint64_t total_primes = 0;
    // XXX: atomic<uint16_t> max_needed_Q(0);

    CudaSieve sieve;

    while (interval_start < n_1) {
        {
            sieve.top = interval_end - 1;
            sieve.bottom = interval_start;
            // TODO use sieve.setFlagOn(<N>)
            // silent
            sieve.flags[30] = 0;
            // flags 29 and 20 are set but seem unused.

            // Init d_primeOut
            sieve.flags[0] = (interval_start == n_0);
            //
            sieve.flags[1] = (interval_start > n_0);

            sieve.launchCtl();

            uint64_t primes_in_range = sieve.kerneldata.getCount();

            safeCudaFree(sieve.d_primeList);

            printf("\tProcessing [%lu, %lu) with %lu primes\n\n",
                    interval_start, interval_end, primes_in_range);

            total_primes += primes_in_range;
        }

        /* Call CUDA Kernel. */
        // kernel_goldbach_search<<<config.threads>>>(
        //    gpu_cleared,
        //    interval_start);


        // Copy overlap from the end back to the start for the next set of threads.
        //CUDA_CHECK(cudaMemcpy((void *)gpu_cleared, gpu_cleared + INTERVAL_BYTES, OVERLAP_BYTES, cudaMemcpyDeviceToDevice));

        interval_start += interval_size;
        interval_end = std::min(interval_start + interval_size, n_1);
    }

    // clean up
    //CUDA_CHECK(cudaFree(gpu_cleared));
    //CUDA_CHECK(cudaEventDestroy (global_start));
    //CUDA_CHECK(cudaEventDestroy (batch_start));
    //CUDA_CHECK(cudaEventDestroy (stop));
    printf("\tTotal Primes: %lu\n", total_primes);
}
