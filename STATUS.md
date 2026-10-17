> **Status: under development — pre-alpha.** Interfaces, numbers and file layout can change without notice.

## Works today
- V1: PASS (2026-10-17); 129600 synthetic blocks match the integer twin cleanly; faults can produce wrong samples. Evidence: `media/idct_rings.py`.
- V2: PASS (2026-10-17); Historical 134.1x slowdown; 1 and 113 mismatching samples under injected faults. Evidence: `media/idct-rings.json`.
- V3: PASS (2026-10-17); 200000 synthetic vectors match the plain integer twin. Evidence: `media/vp9-idct-rings/src/main.rs`.
- V4: PASS (2026-10-17); 20000 random pairs plus edge cases per operation match the native binary64 twin. Evidence: `jsnum/src/lib.rs`.
- V5: PASS (2026-10-17); 1000000 operations match native; injected faults caught. Evidence: `jsnum/fast_probe.c`.
- V6: EXTERNAL prerequisite; ARMILLA-WIRE-0 is a separately licensed renderer integration, outside this workload export. Evidence: `PROVENANCE.md`.

All local rows passed from the clean staged tree on 2026-10-17; external rows were not run.

V7: PASS (2026-10-17), historical negative checked by suite command; evidence media/vp9-idct-rings/result-three-1M.json.

Stage: about 30 %. Dated rerun output: `results/2026-10-17-quicktests.log`. Historical results are not new benchmarks.

## Not yet tested
Security audit, constant-time conversions, external engine/hardware integrations, and end-to-end application guarantees are outside this local suite.

## Known failures and honest negatives
The 2026-10-17 rerun also yielded 0 clean mismatches, then 1 and 113 fault-test mismatches; its timing is recorded on the loaded shared host (results/2026-10-17-quicktests.log). Timings are host/load dependent. Historical media/idct-rings.json: clean 0 mismatches on 129600 blocks, 134.1x slower than NumPy (0.16 frames/s); 1% fault test has 1 mismatch and 4 no-majority values; 10% has 113 mismatches and 175 no-majority values. No claim of fault-proof decoding. JS scalar binary64 arithmetic is much slower than native. Full bitstream parsing, codec conformance, real playback and DRM are not established by these probes.


VP9 three-group negative (V7): historical `media/vp9-idct-rings/result-three-1M.json` has 17134 wrong outputs after 200000 injected faults, all in group C. The quick run uses 200000 vectors because the source prints a fixed 200000 fault label even for smaller inputs; this export does not patch the read-only source.

## Corrections
Exactness is against the stated twin, not a security or conformance certification. Negative rows explicitly retain the original failures.

## How to re-check
Use README.md setup, then each runnable suite.json verify command from the repository root. Commands compile/write into a fresh temporary directory. The renderer and math-core index/specs are external prerequisites.

## Defensive publication
Published as defensive prior art: the author claims no patent on what is described here and asks that none be granted on it. DOI pending the first GitHub release.


## Findings update — 2026-10-17

See [FINDINGS-2026-10-17.md](FINDINGS-2026-10-17.md) for later implemented opt-in APIs, exact engineering/scientific limits and a portable source-only smoke check. Earlier statements above remain the dated historical record; they are scoped to the earlier bundle. The new capsule is not a full research-gate reproduction or deployed integration.
