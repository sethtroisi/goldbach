# Goldbach

Reasonable implementation of naive algorithm for verifying Goldbach conjector.

## Algorithm

Uses a naive algorithm with one or two clever bithacks.

The naive algorithm is two nested loops.
```

for p in primes:
  for q in small primes:
    found[p + q] = 1

for i < target:
  if found[i] == 0:
    print("Unknown!")
```


### Replacing the inner loop with a bitmask.

The inner loop (`q`) iterates over 1229 primes less than 10000.

A sizable (10x) speed up is achieved by replacing iteration over all small
primes with a bitmask of length 10000 representing small primes.
This bitmask is then shifted over by all prime p and OR'd with `found`.

The implementation needs to pre-compute a 157 element bitmask
(64 * 157 > 10000) representing the values of q. Then 64 copies of that are
needed to account for `p mod 64`.

### Threading

This problem is trivially parallelized by splitting found to ranges.

A small amount of care is needed for ranges near the end of the interval.

### Combining p

the gap between `p_i` and `p_{i+1}` is small by pre-computing the bitmask
array for different gaps (2, 4, 6, 8), multiple p could be handled at once.

TODO: Implement this. Maybe also handling more than 2 at a time.


## Benchmark

```
time ./goldbach -t8 -K 1024 -e 100''000''000''000
real	0m6.670s
user	0m52.846s
sys	0m0.088s

time ./goldbach -t12 -K 1500 -e 1''000''000''000''000
  267468893942 = 267468891139 + 2803 (407)
  926868341768 = 926868338701 + 3067 (437)
  599533546358 = 599533542901 + 3457 (481)
real	1m6.552s
user	12m20.815s
sys	    0m0.109s
```
