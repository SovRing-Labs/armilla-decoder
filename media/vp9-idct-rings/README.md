# vp9-idct-rings — VP9 inverse DCT on Armilla rings (MEDIA-PLAYER-0 Armilla lead, J 10-15 readback confirmed)
Spec: VP9 Bitstream & Decoding Process Specification v0.6 §8.7.1 (research/specs/google.vp9-bitstream). Run: `cargo build --release && ./target/release/vp9-idct-rings 1000000`.
Plain = §8.7.1 on i32. Rings = residues on 251·241·239·233·229 (≈2^39.5) + a 2^16 ring; Round2(x,14): low 14 bits from the 2^16 ring,
exact division on the odd rings, then the 2^16 residue rebuilt by mixed-radix base extension (the bridge).
| Date | Vectors | Mismatches | Plain ns/vec | Rings ns/vec | Ratio | Note |
|---|---|---|---|---|---|---|
| 2026-10-15 | 1,000,000 (8-point, inputs ±4096) | 0 | 43.2 | 959.6 | 22.2× | scalar u64 `%`, base extension at every Round2 (13 ops → 16 rounds per vector) |
| 2026-10-15 | 1,000,000 | 0 | 42–45 | 529–550 | 11.7–13.0× | **mode 2, routed rings** (J: "an organization array in the rings, routed so the rounding is removed"): moduli 3·43·127 (2^14 ≡ +1) and 5·29·113 (2^14 ≡ −1), the divisors of 2^28−1 → Round2 = subtract low bits + route straight / negated; low bits from a 2^64 ring; NO base extension |
| 2026-10-15 | 200,000 faults (1 bit of one 2^64-ring value, random op / slot / bit 0–31) | — | — | — | — | 199,895 outputs changed, **199,895 caught** by the routed rings disagreeing (100 %) |
Next lever: residues ≤ 127 fit u8, products < 2^14 → per-ring lookup tables + AVX2 u8/u16 lanes (32 / 16 per op) instead of scalar u64 `%`.
| 2026-10-15 | 1,000,000 | 0 (0 false alarms) | 47 | 8,437 | 176× | **mode 3, C5 three-group format on T50 primes** (J: "the 3x3x3 with math core, output all routed"): A = routing rings 3·43·127·5·29·113, B = T50 7·11·13·17·19·23 (Round2 = × inv 2^14, no bridge), C = 2^64 guard; centre decodes each group (u128 CRT, scalar — the 176× is mostly this decode) and votes 2-of-3 |
| 2026-10-15 | 200,000 faults (⅓ per group: one residue +1 in A or B, one bit 0–63 in C) | — | — | — | — | detected 197,630; **A faults 0 wrong / 66,775, B faults 0 wrong / 66,447 (vote corrects)**; **C faults 17,134 wrong / 66,778** — A and B take their Round2 low bits from C, so a C fault enters all three groups: the groups are not independent |

## J's layout (10-15, readback confirmed): 3 groups × 3 diamonds × (1-3-5-3-1) = 117 T50-class systems → centre = math core
`src/bin/armilla_grid.rs` — `./target/release/armilla_grid FAULT_RATE THREADS`. Each system: own 50 disjoint odd primes (system 0 = T50 3…233,
largest prime 57,751) + own 2^64 low-bit ring; the full 8×8 2D inverse DCT (§8.7.2) stays in residues; centre checks all 50 residues of
every output against the carrier, applies the final Round2(x,5), recomputes a failing block. Frame: 1080p luma, 32,400 blocks.
| Run | Blocks wrong | Faults injected | Repaired by centre | Armilla s/frame | Plain s/frame (1 thread) | Ratio |
|---|---|---|---|---|---|---|
| clean, 8 threads | 0 | 0 | 0 | 0.944 | 0.023 | 41× |
| 10 % blocks faulted, 8 threads | 0 | 3,250 | 3,249 (1 fault had no effect on any output) | 1.066 | 0.025 | 43× |
| clean, 1 thread | 0 | 0 | 0 | 3.991 | 0.021 | 193× |
Layers 1-3-5-3-1 are a share assignment here (each system does whole blocks); they do not yet have different roles.

## Per-diamond ring placement + sub-cores (J 10-15, readbacks confirmed): layout A = 9 diamond cores + final (10), layout B = 3 group cores + final (4)
`src/bin/armilla_cores.rs` — `./target/release/armilla_cores A|B FAULT_RATE WORKERS`. Diamond 0 = routing primes 127·113·43 (data) + 29 (check);
diamonds 1–8 = 2 data + 1 check primes 1031…2903, disjoint; each system keeps its own 2^64 ring. Sub-core: residue decode vs 2^64 copy,
check ring decides, in-place correction; recompute only when neither copy is confirmed; final core verifies a checksum per block.
| Layout | Cores | Faults (10 % of blocks) | Corrected in place | Recomputed | Wrong blocks | Frames/s (8 workers) | vs plain 1 thread |
|---|---|---|---|---|---|---|---|
| A | 10 | 0 | — | — | 0 | 7.4–8.4 | 4.2–5.1× |
| A | 10 | 3,250 | 2,467 | 783 | 0 | 7.2–7.4 | 5.8–6.6× |
| B | 4 | 0 | — | — | 0 | 7.6–8.2 | 5.8–6.1× |
| B | 4 | 3,250 | 2,467 | 783 | 0 | 7.9–8.4 | 4.1–5.7× |
Recomputed ≈ the faults that hit a 2^64 ring (expected share 1/(rings+1) ≈ 795): a carrier fault also bends the rings' later roundings, so neither copy is confirmed.
Laptop load varied between runs (ranges are two runs each). 50 primes → 3–4 rings per system: ~1 → ~8 frames/s.

## Router (J 10-15 "yes add the router"): `src/bin/armilla_router.rs` — `./target/release/armilla_router FAULT_RATE`
Tile column → group; block class (DC-only / sparse ≤ 3 rows / dense) → exact path (zero rows skipped; DC-only = one row + one column,
replicated); within the group the diamond with the lowest measured cost (startup calibration, C3-style) + current load. Group cores (B) + final.
Block mix is SYNTHETIC (40 % DC, 33 % sparse, 28 % dense; detail clustered in the right-hand tile). Two runs each, clean:
| Mode | Frames/s | Wrong | Slowest group busy (s) |
|---|---|---|---|
| fixed (band split, full transform) | 8.6–8.7 | 0 | 0.107–0.108 |
| fixed + class paths | 10.7–10.9 | 0 | 0.065–0.073 |
| routed (tile→group, class path, cost+load diamond) | 10.5–11.5 | 0 | 0.068–0.079 (the busy tile's group) |
10 % faults, routed: 0 wrong, 2,493 corrected in place, 733 recomputed, 0 checksum failures.
Reading: the gain is the class paths; routing balances diamonds inside a group, but tile → group pins the busy tile to one group, which
becomes the bottleneck. Next: let blocks that have no intra dependency (inter blocks) spill to the least-loaded group, or size tiles by cost.

## Simplex core + cyclone counting (J 10-15, readbacks "correct on both"): `src/bin/armilla_simplex.rs` — `./target/release/armilla_simplex FAULT_RATE`
117 systems = vertices of a regular 116-simplex (all pairwise inner products −1/116, checked); superblock (64×64) routing; 30 % intra superblocks
pinned to their tile's group (SYNTHETIC); 20 outer T50 moduli 131…233 count traffic per link (router→group→diamond→system); bottleneck =
hottest link / mean link per level (1.00 = balanced). Ternary query = Σ desirability·key, quantised; sift = nearest vertex.
| Router | Frames/s (3 runs) | Wrong | Bottleneck group / diamond level |
|---|---|---|---|
| fixed split | 8.7–9.1 | 0 | 1.25–1.33 / 1.42–1.45 |
| direct (cost + outer-ring counters) | 11.0–11.6 | 0 | 1.00 / 1.01–1.02 |
| **simplex keys, ternary query** | **11.1–11.8** | 0 | 1.01 / 1.03–1.04 |
| random ternary keys (control) | 6.4–6.7 | 0 | 2.1–2.3 / 5.6–6.1 |
10 % faults: every router 0 wrong (≈2,450 corrected in place, ≈770 recomputed, 0 checksum failures).
Reading: the simplex geometry routes as well as the exact chooser through a ternary query; random keys under the same query pile load onto
one diamond (bottleneck 6×) — law R1 holds here. Agreement with the direct pick is low (4 % system, ~20 % diamond) because many systems
tie in cost; balance is what matters. The router is compiled Rust (not yet JIT-emitted by the math core).

## Step 2 — SIMD lanes (J 10-15 "implement all 3"): `src/bin/armilla_simd.rs` — `[CORE_INLINE=1] ./target/release/armilla_simd FAULT_RATE`
16 blocks per batch as lane arrays per ring; mod m by a float reciprocal + exact integer correction (vectorised: vpmulld / vpmuludq / vcvtdq2ps);
DC-only blocks keep the scalar shortcut; sub-core fast path = compare the 2^64 copy with each ring (CRT only on a mismatch).
CORE_INLINE=1 = layout A (math core inside each diamond); default = layout B (group cores). 10 % faults, 2 runs each:
| | Fixed split | Direct router | Simplex router |
|---|---|---|---|
| Layout B (group cores) | 15.2–17.2 | 13.3–20.9 | 17.6–21.1 |
| Layout A (diamond cores) | 17.6–19.4 | 16.4–20.7 | 15.8–18.6 (one outlier 7.6) |
All runs: 0 wrong blocks, ≈2,470 corrected in place, ≈755 recomputed, 0 checksum failures. Before step 2: 8.7–11.8.
Not reached: 30 frames/s. Lane work halved diamond busy time (≈0.065 → 0.031–0.047 s); the rest is per-block copies / allocations between
the lanes, the cores and the channels, and 13 threads on 8 hardware threads. Next levers: keep outputs in lane form to the core, one
allocation per batch, cores per diamond without a channel hop.

## Step 1 — JIT router (J 10-15 "implement all 3"): `src/bin/armilla_jit.rs` — `./target/release/armilla_jit FAULT_RATE`
The centre's router step is emitted by the engine's own JIT emitter (omni-render `src/jit/ternary.rs` JITW-0, included by path: +1 → add,
−1 → sub, 0 → nothing) from the ternarised key matrix at startup. Simplex keys ternarised = the identity (117 non-zeros, emitted in 14 µs):
on a regular simplex the nearest vertex is the largest query coordinate. Random ternary keys: 9,125 non-zeros, 166 µs.
| Router | Route µs / superblock | Frames/s | Wrong | Bottleneck (diamond level) | JIT decisions = interpreted |
|---|---|---|---|---|---|
| direct (cost + counters) | — | 16.7 | 0 | 1.01 | — |
| simplex, interpreted | 10.26 | 21.8 | 0 | 1.06 | — |
| **simplex, JIT** | **0.56** | 21.1 | 0 | 1.06 | **identical** |
| random keys, interpreted | 9.88 | 15.1 | 0 | 4.35 | — |
| random keys, JIT | 1.74 | 14.2 | 0 | 4.35 | identical |
(10 % faults; one run; frames/s noisy ±3 on this laptop.)

## 1-3-7-3-1 with primes by layer size (J 10-15 19:40, readback confirmed, "yes, explore") — reading V (vertical): `src/bin/armilla_layers.rs`
`./target/release/armilla_layers base|v137 FAULT_RATE K [WORKERS]`. v137: one value spans the diamond's 15 systems, one ring each —
middle 7 = information rings 3·5·7·11·13·17·19 (≈ 2^22.2, all byte-shuffle sized), upper 3 = routing primes 29·43·127, lower 3 = 113·241·251,
top = 2^31 − 1, bottom = 2^64 carrier. 7 redundant odd moduli, each larger than every information modulus (the RRNS condition) → up to 3
wrong rings per value corrected in place. base = the 1-3-5-3-1 per-diamond placement of armilla_cores.rs (one check ring) in the same binary.
10 % of blocks get K wrong rings (K distinct rings of one value; the carrier counts as a ring). Layout B cores. CELT search running on all
8 cores during these runs (load ≈ 14): counts are exact, seconds are only relative.
| K wrong rings | Layout | Wrong blocks | Faulted | Hit the carrier | Corrected in place | Recomputed | s/frame (loaded) |
|---|---|---|---|---|---|---|---|
| 1 | base 1-3-5-3-1 | 0 | 3,250 | 829 | 2,421 | 829 | 0.78 |
| 1 | v137 | 0 | 3,250 | 217 | 2,868 | 217 | 2.10 |
| 2 | base | 0 | 3,217 | 1,582 | **0** | 3,217 | 0.71 |
| 2 | v137 | 0 | 3,234 | 405 | **2,824** | 405 | 2.00 |
| 3 | base | 0 | 3,189 | 2,335 | **0** | 3,189 | 0.64 |
| 3 | v137 | 0 | 3,219 | 633 | **2,586** | 633 | 1.64 |
Clean (no faults): both 0 wrong; base 0.59 s, v137 2.32 s, armilla_cores B 0.56 s (same load).
Reading: the size-tiered diamond corrects 2 and 3 wrong rings in place where the one-check diamond can only recompute; v137's recomputes
are EXACTLY the faults that hit the 2^64 carrier (a carrier fault bends every later rounding, so no ring set can repair it — next: a second
low-bit source). Scalar speed: 14 rings cost ~3–4× the 4-ring diamond; the middle layer's rings are all ≤ 19, i.e. byte-shuffle sized
(shuffle_probe.c: constant multiply = 1 instruction / 32 values), which is the lever not used yet.

### D3 calibration — the butterfly step on one ring, by prime size class (`projects/armilla/exec/butterfly_probe.c`, AVX2, one core)
Result: `projects/armilla/results/d3-butterfly-calibration-loaded.json` — measured WHILE the CELT search held all 8 cores (load ≈ 14):
±50 % noise; rerun on a quiet machine before any placement is fixed. All classes 0 wrong vs the scalar reference.
| Class | Primes | ns per value × ring | Bits per ring | ns per bit of range |
|---|---|---|---|---|
| S: byte lanes, every constant multiply = 1 shuffle | 3 · 5 · 7 · 11 · 13 · 17 | 0.51–1.34 | 1.6–4.1 | 0.15 (17) … 0.45 (3) |
| M16: 16-bit lanes, Barrett | 19 · 29 · 43 · 113 · 127 · 241 · 251 | 1.28–2.73 | 4.2–8.0 | 0.17 (241) … 0.56 (29) |
| M32: 32-bit lanes, float reciprocal | 1031 · 4093 · 46337 | 4.5–5.7 | 10–15.5 | 0.30–0.57 |
| L: 2^31 − 1, 64-bit lanes, Mersenne fold (shifts only) | 2^31 − 1 | 5.69 | 31 | 0.18 |
| L%: same prime, scalar u64 `%` | 2^31 − 1 | 104 | 31 | 3.36 |
Per value and butterfly, all rings: v137 ≈ 22.5 ns (corrects 3 wrong rings) · 1-3-5-3-1 diamonds 1–8 (1031…4093) ≈ 15.4 ns · diamond 0
(127·113·43·29) ≈ 7.2 ns (both correct 1). Reading: per bit of range the cheapest rings are 17 (shuffle), 241/251 (top of the 16-bit lane)
and 2^31 − 1 (Mersenne fold); the tiny primes 3, 5 and the 1031…4093 range are the most expensive per bit — a byte lane holding mod 3
carries 1.6 of its 8 bits.

## Mixed systems (J 10-15 20:05, readbacks confirmed, "let's see what we get"): `armilla_layers mix1|mix2`
Every system of the 1-3-7-3-1 diamond = its main size class (information rings: the value range) + one reserved ring of EACH other size
class (checks) + 1 (mix1) or 2 (mix2) rounding helpers; each system computes whole blocks; blocks go to systems by weight 1/(rings + helpers).
middle 7 = small main (5–7 primes from 3…31, product ≥ 2^22) + 1 medium (223…251) + 1 large (2^31 − c) · upper/lower 3 = 3 medium main
(113…211) + 1 small + 1 large · top/bottom = 1 large main + 1 medium (107/109) + 1 small (3/5). With 2 helpers each feeds its OWN rings:
helper A the main rings' rounding, helper B the reserved rings' — a helper fault bends only its own side, the other side is an
independent witness. Decoder: legal range |v| ≤ 2^16 (inputs ±512, ≤ 8× growth per 8-point pass → |v| ≤ 2^15).
`result-layers-mix.jsonl` (same frame, 10 % of blocks faulted with K wrong slots; CELT search on all cores → seconds are relative):
| K | Layout | Wrong | Fixed in place | Recomputed | s/frame scalar (loaded) | est. ns/value vectorised (D3 table) |
|---|---|---|---|---|---|---|
| 0 | base 1-3-5-3-1 | 0 | — | — | 0.66 | 15.8 (diamonds 1–8) / 7.7 (diamond 0) |
| 0 | v137 (one value spans the diamond) | 0 | — | — | 1.97 | 23.0 |
| 0 | mix1 | 0 | — | — | 1.11 | 12.3 |
| 0 | mix2 | 0 | — | — | 1.07 | 13.0 |
| 1 | base / v137 / mix1 / **mix2** | 0 / 0 / 0 / 0 | 2,421 / 2,868 / 2,649 / **3,168** | 829 / 217 / 520 / **5** | | |
| 2 | base / v137 / mix1 / mix2 | 0 / 0 / 0 / 0 | 0 / **2,824** / 153 / 2,727 | 3,217 / 405 / 3,075 / 519 | | |
| 3 | base / v137 / mix1 / mix2 | 0 / 0 / 0 / 0 | 0 / **2,586** / 9 / 1,893 | 3,189 / 633 / 3,182 / 1,317 | | |
Reading: mix2 fixes 97.5 % of single faults in place (helper faults included — the split-fed helpers remove the carrier problem) at the
lowest vectorised cost of the self-checking layouts (13.0 ns vs 15.8 for today's diamonds); v137 (7 redundant rings per value) is the
strongest at 2–3 wrong rings but costs ~1.8× mix2. Way there (kept for the record): the first mix decoder trusted a system's own ring span
as the legal range (a single 2^31 ring accepts almost anything) and let helper A, which fed every ring, be "confirmed" by one main ring →
37–1,341 wrong blocks; fixed by the input-derived legal range, split-fed helpers, majority-of-main-rings confirmation.
