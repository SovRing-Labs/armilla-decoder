// JS-ARMILLA-0 probe F1 (J 2026-10-16: "verify what you've said, our armilla should be able to do the arithmetic and the routing").
// The Armilla doing JS int32 arithmetic at speed: the value lives on the ring Z/2^64 (the rounding helper — the CPU register's own
// arithmetic, so the op is exact by construction), two guard rings 251 and 241 carry residues. add / sub / mul: the guards update
// homomorphically (r = (ra op rb) mod p); xor / and / or / shifts are not ring maps: the guard is re-derived from the value (constant
// modulo = multiply-shift). After every op the guard checks value mod p == r. Workload: YouTube-VM-like chains of int32 ops on a small
// register file, ops drawn at random; the same chain run natively. Faults: one bit of the value flipped after a random op.
// Build: cc -O2 -march=native fast_probe.c -o fast_probe ; ./fast_probe [N]
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
static double now() { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
#define R 16
typedef struct { uint8_t op, a, b, d; } Ins;
static inline int32_t js32(uint64_t v) { return (int32_t)(uint32_t)v; }   // ToInt32 of the ring value
// residue of the int32 value v (its two's-complement integer) mod p: shift into [0, 2^31·2) first so the residue is of the signed value
static inline uint32_t res(int32_t v, uint32_t p) { int64_t x = v; int64_t r = x % (int64_t)p; return (uint32_t)(r < 0 ? r + p : r); }
int main(int argc, char **argv) {
  long N = argc > 1 ? atol(argv[1]) : 20000000; srand(7);
  Ins *P = malloc(sizeof(Ins) * N);
  for (long i = 0; i < N; i++) { P[i].op = rand() % 7; P[i].a = rand() % R; P[i].b = rand() % R; P[i].d = rand() % R; }
  // native
  int32_t n[R]; for (int i = 0; i < R; i++) n[i] = 0x9e3779b9u * (i + 1);
  double t0 = now();
  for (long i = 0; i < N; i++) { Ins c = P[i]; int32_t x = n[c.a], y = n[c.b], z;
    switch (c.op) { case 0: z = (int32_t)((uint32_t)x + (uint32_t)y); break; case 1: z = (int32_t)((uint32_t)x - (uint32_t)y); break;
      case 2: z = (int32_t)((uint32_t)x * (uint32_t)y); break; case 3: z = x ^ y; break; case 4: z = x & y | (x >> 3); break;
      case 5: z = (int32_t)((uint32_t)x << (y & 31)); break; default: z = (int32_t)((uint32_t)x >> (y & 31)); }
    n[c.d] = z; }
  double tn = now() - t0;
  // Armilla: value on Z/2^64 (low 32 bits = the JS int32), guards mod 251 and 241
  uint64_t v[R]; uint32_t g1[R], g2[R]; long alarms = 0;
  for (int i = 0; i < R; i++) { v[i] = (uint32_t)(0x9e3779b9u * (i + 1)); g1[i] = res(js32(v[i]), 251); g2[i] = res(js32(v[i]), 241); }
  t0 = now();
  for (long i = 0; i < N; i++) { Ins c = P[i]; int32_t x = js32(v[c.a]), y = js32(v[c.b]), z; uint32_t r1, r2;
    switch (c.op) {
      case 0: z = (int32_t)((uint32_t)x + (uint32_t)y); r1 = (g1[c.a] + g1[c.b]) % 251; r2 = (g2[c.a] + g2[c.b]) % 241; break;
      case 1: z = (int32_t)((uint32_t)x - (uint32_t)y); r1 = (g1[c.a] + 251 - g1[c.b]) % 251; r2 = (g2[c.a] + 241 - g2[c.b]) % 241; break;
      case 2: z = (int32_t)((uint32_t)x * (uint32_t)y); r1 = g1[c.a] * g1[c.b] % 251; r2 = g2[c.a] * g2[c.b] % 241; break;
      case 3: z = x ^ y; r1 = res(z, 251); r2 = res(z, 241); break;
      case 4: z = x & y | (x >> 3); r1 = res(z, 251); r2 = res(z, 241); break;
      case 5: z = (int32_t)((uint32_t)x << (y & 31)); r1 = res(z, 251); r2 = res(z, 241); break;
      default: z = (int32_t)((uint32_t)x >> (y & 31)); r1 = res(z, 251); r2 = res(z, 241); }
    // the ring result: add / sub / mul wrap mod 2^32 in JS — the guard residues are of the unwrapped integer, so re-base them on the
    // wrapped value only when it wrapped (the 2^32 carry = the routing between the value ring and the guards)
    if (c.op <= 2) { int64_t wide = c.op == 0 ? (int64_t)x + y : c.op == 1 ? (int64_t)x - y : (int64_t)x * y;
      if (wide != z) { r1 = res(z, 251); r2 = res(z, 241); } }
    v[c.d] = (uint32_t)z; g1[c.d] = r1; g2[c.d] = r2;
    if (res(js32(v[c.d]), 251) != r1 || res(js32(v[c.d]), 241) != r2) alarms++; }
  double ta = now() - t0;
  int same = 1; for (int i = 0; i < R; i++) same &= js32(v[i]) == n[i];
  // faults: flip one bit of the value after the op, before the guard check
  long injected = 0, caught = 0;
  for (int i = 0; i < R; i++) { v[i] = (uint32_t)(0x9e3779b9u * (i + 1)); g1[i] = res(js32(v[i]), 251); g2[i] = res(js32(v[i]), 241); }
  for (long i = 0; i < N / 10; i++) { Ins c = P[i]; int32_t x = js32(v[c.a]), y = js32(v[c.b]), z = x ^ y;
    if (c.op == 0) z = (int32_t)((uint32_t)x + (uint32_t)y); uint32_t r1 = res(z, 251), r2 = res(z, 241); uint64_t val = (uint32_t)z;
    if (i % 7 == 3) { val ^= 1ull << (rand() % 32); injected++; if (res(js32(val), 251) != r1 || res(js32(val), 241) != r2) caught++; }
    v[c.d] = (uint32_t)z; g1[c.d] = r1; g2[c.d] = r2; }
  printf("{\"ops\": %ld, \"native_ns_per_op\": %.2f, \"armilla_ns_per_op\": %.2f, \"ratio\": %.2f, \"results_equal_native\": %s, \"guard_alarms_clean\": %ld, \"faults_injected\": %ld, \"faults_caught\": %ld}\n",
         N, tn / N * 1e9, ta / N * 1e9, ta / tn, same ? "true" : "false", alarms, injected, caught);
  return 0;
}
