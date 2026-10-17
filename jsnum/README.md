# Armilla JS-number unit (J 10-15: "the Armilla should be able to do any kind of arithmetic … it's a matter of organizing the rings, coprimes, how they interact within each unit and group, and the routing")
`src/lib.rs` (the unit) · `src/bin/verify.rs` (`cargo build --release && ./target/release/verify 2000000`). Semantics: ECMA-262 Number =
IEEE 754 binary64 roundTiesToEven (−0, ±∞, NaN, subnormals, overflow), ToInt32 + & | ^ << >> >>>, <. Reference = the CPU's own f64.
## Ring organisation of one unit
| Ring(s) | Role | Why there |
|---|---|---|
| Z/2^128 | carries the exact value: rounding split, leading bit, bitwise, shifts | low bits are direct reads — no bridge |
| 3·43·127 (2^14 ≡ +1), 5·29·113 (2^14 ≡ −1) | routing rings: ×/÷ 2^k = rotation of a 14-cycle; guard of each op's integer core | alignment, shifts and the rounding split are routing |
| M61 = 2^61 − 1 | second guard; ×2^k = 61-bit rotation | routing rings are blind when an operand is a multiple of 2^28 − 1 (repeating binary fractions like 4/17): 13 faults slipped through without it |
| math core | division = long division on Z/2^128, verified q·d + r = n on the guards | division is not a ring map |
## Results (2026-10-15, this laptop, 2,000,000 cases per op: random bit patterns over the whole range + near pairs + 25×25 edge values)
| Op | Mismatches vs IEEE | Armilla ns/op | Native ns/op |
|---|---|---|---|
| add / sub | 0 / 0 | 311–335 | 2.3–2.5 |
| mul | 0 | 411–449 | 2.6–2.8 |
| div | 0 | 690–787 | 2.4–2.7 |
| and / or / xor / shl / sar / shr / lt | 0 each | — | — |
Guard: 24.9 M checks on clean runs, 0 false alarms. Faults (one bit of the 2^128 value flipped after an op's integer core, random op / bit):
358,210 outputs changed, **358,210 caught** (routing rings alone: 358,197 — the 13 misses were divisions by repeating fractions, fixed by M61).
Speed: scalar u128 `%` per ring — 130–310× native; the lanes / batching lever from the VP9 grid applies next.
## F1 — the fast organisation (2026-10-16, JS-ARMILLA-0; J: "verify what you've said, our armilla should be able to do the arithmetic and the routing")
`fast_probe.c`: the value on the ring Z/2^64 (the rounding helper = the register's own arithmetic, exact by construction), guard rings
251 and 241 (add / sub / mul homomorphic, re-based when the 2^32 wrap carries; xor / and / or / shifts re-derived by constant modulo),
checked after every op. YouTube-VM-like int32 chains (7 op kinds, 16 registers), 20 M ops × 3 runs (`result-fast-probe.jsonl`):
**11.1–11.4 ns/op vs native 7.9–8.1 (1.4×), results identical, 0 false alarms, 285,714 / 285,714 single-bit faults caught.**
The 130–310× above is the scalar u128-`%`-per-ring build, not a property of the Armilla; our JS interpreter spends ~100 ns per step.
