# armilla-decoder — Exact decoder mathematics on residue rings

[![status: pre-alpha](https://img.shields.io/badge/status-pre--alpha%20(under%20development)-orange)](STATUS.md)
[![defensive publication](https://img.shields.io/badge/defensive%20publication-TECHNICAL--DISCLOSURE.md-blue)](TECHNICAL-DISCLOSURE.md)

Status: **under development — pre-alpha** (`STATUS.md`). Defensive publication: `TECHNICAL-DISCLOSURE.md` (dated 2026-10-16; the DOI is
minted by Zenodo on the first GitHub release).


H.264 and VP9 inverse-transform probes plus JavaScript number arithmetic expressed on rings. The shipped references compare synthetic inputs with their own plain integer or CPU twins. This is a reproducible arithmetic workload, not a full bitstream decoder or playback system.

## Run it
```bash
python3 -m venv .venv && . .venv/bin/activate && pip install numpy
PYTHONDONTWRITEBYTECODE=1 python3 tests/check.py h264
bash install.sh --no-gui
```
Additional rows need the compilers named in suite.json. No network fetch is needed once prerequisites are installed.

## What it does today
| Row | Check | Evidence |
|---|---|---|
| V1 | 129600 synthetic blocks match the integer twin cleanly; faults can produce wrong samples. | `media/idct_rings.py` |
| V2 | Historical 134.1x slowdown; 1 and 113 mismatching samples under injected faults. | `media/idct-rings.json` |
| V3 | 200000 synthetic vectors match the plain integer twin. | `media/vp9-idct-rings/src/main.rs` |
| V4 | 20000 random pairs plus edge cases per operation match the native binary64 twin. | `jsnum/src/lib.rs` |
| V5 | 1000000 operations match native; injected faults caught. | `jsnum/fast_probe.c` |
| V6 | ARMILLA-WIRE-0 is a separately licensed renderer integration, outside this workload export. | `PROVENANCE.md` |

## Not shown here
Historical media/idct-rings.json: clean 0 mismatches on 129600 blocks, 134.1x slower than NumPy (0.16 frames/s); 1% fault test has 1 mismatch and 4 no-majority values; 10% has 113 mismatches and 175 no-majority values. No claim of fault-proof decoding. JS scalar binary64 arithmetic is much slower than native. Full bitstream parsing, codec conformance, real playback and DRM are not established by these probes.

## Files
| Path | Contents |
|---|---|
| exec/ or media/ and jsnum/ | workload implementations |
| tests/check.py | portable verification wrapper |
| results/ | historical and dated local evidence |
| suite.json | runnable and external claim rows |

## Licence and stage
PolyForm Small Business 1.0.0 OR PolyForm Noncommercial 1.0.0 + Free Trial; CC BY 4.0 documents. See LICENSING.md. Stage: about 30 %. The common abacus engine is cited in PROVENANCE.md; it is not copied.

VP9 three-group negative (V7): historical `media/vp9-idct-rings/result-three-1M.json` has 17134 wrong outputs after 200000 injected faults, all in group C. The quick run uses 200000 vectors because the source prints a fixed 200000 fault label even for smaller inputs; this export does not patch the read-only source.


## Findings update — 2026-10-17

See [FINDINGS-2026-10-17.md](FINDINGS-2026-10-17.md) for later implemented opt-in APIs, exact engineering/scientific limits and a portable source-only smoke check. Earlier statements above remain the dated historical record; they are scoped to the earlier bundle. The new capsule is not a full research-gate reproduction or deployed integration.
