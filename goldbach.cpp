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

void goldbach_search(const struct Config &config);
void goldbach_search_single(const Config &config);

constexpr uint32_t THREAD_SPACING = 256'000;
constexpr uint32_t OVERLAP = 1024 + 512;

// See https://oeis.org/A025019 and
// "EMPIRICAL VERIFICATION OF THE EVEN GOLDBACH CONJECTURE" by TOMAS OLIVEIRA e SILVA
constexpr uint32_t MAX_Q = 10'000;
// must be greater than primepi(2 * THREAD_SPACING + MAX_Q)
constexpr uint32_t PRIMES_SIZE = []() constexpr -> uint32_t {
    uint32_t RANGE = 2 * THREAD_SPACING + MAX_Q;
    uint32_t PRIME_PI = 1.1 * RANGE / log(RANGE);
    return 1024 * ((PRIME_PI + 1024 - 1) / 1024);
}();
constexpr uint32_t MAX_PRIME_GAP = 772; // Good up to 10^20

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

atomic<uint16_t> max_needed_Q(0);
bool set_max(uint16_t q_i) {
    // Use atomic::fetch_max after c++26

    uint16_t was = max_needed_Q.load();
    while (was < q_i) {
        // Try to store
        if (max_needed_Q.compare_exchange_strong(was, q_i)) {
            return true;
        }
        // Was updated in the middle, try again with new value from max_needed_Q
        was = max_needed_Q.load();
    }
    return false;
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
    vector<uint16_t> Q;
    primesieve::generate_primes(3, MAX_Q, &Q);

    const auto n_0 = config.n_0;
    const auto n_1 = config.n_1;
    printf("Searching [%lu, %lu] P <= %u, spacing = %u (overlap: %u, <= %u primes)\n",
            n_0, n_1, config.K, THREAD_SPACING, OVERLAP, PRIMES_SIZE);
    printf("\t|P| = %lu (%lu masks)\n", P.size(), P_masks[31].size());
    printf("\n");
    assert( n_1 >= n_0 );
    assert( (unsigned) config.K <= OVERLAP );
    assert( OVERLAP <= MAX_Q );

    assert( THREAD_SPACING % 64 == 0);
    assert( OVERLAP % 64 == 0);

    size_t THREAD_BYTES = THREAD_SPACING / 8;
    size_t THREAD_UINT64_PER = THREAD_BYTES / 8;
    size_t OVERLAP_BYTES = OVERLAP / (2 * 8);

    uint64_t *cleared = (uint64_t*) calloc(THREAD_BYTES + OVERLAP_BYTES, 1);
    uint64_t *overlap = cleared + THREAD_UINT64_PER;
    uint64_t *cleared_after_overlap = cleared + (OVERLAP_BYTES / 8);

    uint64_t prime_start = 0;
    uint64_t prime_end = 0;

    array<uint64_t, PRIMES_SIZE> circular_primes;

    size_t NON_OVERLAPPED_BYTES = THREAD_BYTES - OVERLAP_BYTES;
    uint64_t FULL_MASK = std::numeric_limits<uint64_t>::max();

    // Break up [n_0, n_1] into ranges of THREAD_SPACING
    uint64_t M = 2 * THREAD_SPACING;
    uint64_t M_COUNT = (n_1 - n_0 + 1 + M-1) / M;
    size_t brute_searches = 0;
    double setup_time = 0.0;
    double search_time = 0.0;
    double brute_time = 0.0;

    // Need primes up to MAX_Q prior to n_0
    uint64_t min_needed_prime = std::max<int64_t>((int64_t)n_0 - MAX_Q - MAX_PRIME_GAP, int64_t(3));
    assert( min_needed_prime >= 3);
    primesieve::iterator it(min_needed_prime);

    // Handle overflow from previous interval and primes needed for brute check
    if (n_0 > 0) {
        assert( n_0 >= OVERLAP );
        uint64_t prime;
        for (prime = it.next_prime(); prime < n_0 - OVERLAP; prime = it.next_prime()) {
            circular_primes[prime_end++] = prime;
        }
        for (; prime < n_0; prime = it.next_prime()) {
            circular_primes[prime_end++] = prime;
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
        assert(prime_end <= PRIMES_SIZE);
    }

    uint64_t prime = it.next_prime();
    for (uint64_t m_i = 0 ; m_i < M_COUNT; m_i ++) {
        uint64_t M_0 = n_0 + M * m_i;
        uint64_t M_1 = M_0 + M-1;
        // [M_0, M_1]

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
