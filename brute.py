import primesieve
import sys


if len(sys.argv) < 2:
    print("requires N to test")
    exit(1)

N = int(sys.argv[1])
print("Searching for", N)
assert N >= 4

# small prime
q_it = primesieve.Iterator()
q = q_it.next_prime()

p_it = primesieve.Iterator()
p_it.skipto(N)
p = p_it.prev_prime()

while q <= N:
    t = p + q
    if t == N:
        print(f"Found {N} = {p} + {q}")
        exit(0)
    if t < N:
        q = q_it.next_prime()
    else:
        p = p_it.prev_prime()

print("NOT FOUND!?!")
