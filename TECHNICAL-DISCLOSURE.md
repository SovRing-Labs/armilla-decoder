# armilla-decoder — technical disclosure

Published as defensive prior art: the author claims no patent on what is described here and asks that none be granted on it.

Author: Juan Carlos De Santiago, with the Sovereign Raptor fleet (AI-assisted)
Date: 2026-10-16

## Problem and method
H.264 and VP9 inverse-transform probes plus JavaScript number arithmetic expressed on rings. The shipped references compare synthetic inputs with their own plain integer or CPU twins. This is a reproducible arithmetic workload, not a full bitstream decoder or playback system.

The H.264 probe encodes each coefficient on twelve small primes, divided into three groups. Additions and subtractions operate per ring. At arithmetic shifts each group reconstructs a signed integer with Garner digits; a two-of-three vote selects a value, applies the shift, and re-encodes. A row butterfly pass, column pass and final (x+32)>>6 produce each 4x4 block. Multiple faults can defeat a majority, and the current probe continues with a selected value; this is not fail-closed correction.

The VP9 probe generates the 8-point butterfly schedule from specification prose, stores residues on odd primes plus a power-of-two rounding helper, and uses exact base extension for shifts. Other shipped probes explore grids, routing, redundancy and SIMD. They remain synthetic arithmetic comparisons, not codec conformance.

The JS-number unit carries exact intermediate values in a 2^128 ring, uses routing primes and a Mersenne guard, rounds to binary64, and checks against CPU arithmetic. The separate fast C probe uses a 2^64 carrier and guards 251 and 241 for int32 operation chains. It is not the full binary64 implementation.

Field sources: ITU-T H.264, Advanced video coding for generic audiovisual services, section 8.5.12; WebM Project, VP9 Bitstream & Decoding Process Specification v0.6, section 8.7; ECMA-262 Number semantics; IEEE 754 binary floating-point arithmetic. These formulas were transcribed from specification prose, not reference implementation code.

## Measurements and limits
Historical media/idct-rings.json: clean 0 mismatches on 129600 blocks, 134.1x slower than NumPy (0.16 frames/s); 1% fault test has 1 mismatch and 4 no-majority values; 10% has 113 mismatches and 175 no-majority values. No claim of fault-proof decoding. JS scalar binary64 arithmetic is much slower than native. Full bitstream parsing, codec conformance, real playback and DRM are not established by these probes.

The runnable procedure and evidence mapping are in suite.json and STATUS.md. This disclosure describes the implemented workload and distinguishes it from unbuilt designs.

VP9 three-group negative (V7): historical `media/vp9-idct-rings/result-three-1M.json` has 17134 wrong outputs after 200000 injected faults, all in group C. The quick run uses 200000 vectors because the source prints a fixed 200000 fault label even for smaller inputs; this export does not patch the read-only source.


## Findings update — 2026-10-17

See [FINDINGS-2026-10-17.md](FINDINGS-2026-10-17.md) for later implemented opt-in APIs, exact engineering/scientific limits and a portable source-only smoke check. Earlier statements above remain the dated historical record; they are scoped to the earlier bundle. The new capsule is not a full research-gate reproduction or deployed integration.
