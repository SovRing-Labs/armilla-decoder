# AGENTS.md — armilla-decoder

For coding assistants (and people) setting this repository up. Plain steps; every command runs from the repo root.

**Status: pre-alpha, under development.** Interfaces, numbers and file layout can change without notice. Read `STATUS.md` first: what works today (each line names its test and evidence file), what is not yet tested, the known failures and the corrections.

## What it is
Stage: about 30 % (under development, pre-alpha).
Exact decoder mathematics on residue rings

## Install and check (3 commands)
```bash
git clone --depth 1 https://github.com/SovRing-Labs/armilla-decoder && cd armilla-decoder
python3 -m venv .venv && . .venv/bin/activate && python -m pip install -q numpy
PYTHONDONTWRITEBYTECODE=1 python3 tests/check.py h264   # expect: PASS
```
Needs: python3, numpy, rustc, cargo, cc.
Whole suite instead: `bash install.sh --no-gui` (python check → venv → every runnable `suite.json` test, headless; without `--no-gui` it then opens the local results page on 127.0.0.1).

GUI: `bash install.sh`

## What the tests mean
`suite.json` lists every claim this repository makes with its verify command, the expected output, what passing means, how the number was obtained and how to use it in your own work (9 tests, in `suite.json` at the repository root). `STATUS.md` says what works today, what is not yet tested and the known failures; the umbrella's `STATUS.md` has one row per repository.

## Where the evidence is
- `README.md`
- `TECHNICAL-DISCLOSURE.md` — the method, as dated prior art
- `CITATION.cff`

## Honest limits (from STATUS.md)
The 2026-10-17 rerun also yielded 0 clean mismatches, then 1 and 113 fault-test mismatches; its timing is recorded on the loaded shared host (results/2026-10-17-quicktests.log). Timings are host/load dependent. Historical media/idct-rings.json: clean 0 mismatches on 129600 blocks, 134.1x slower than NumPy (0.16 frames/s); 1% fault test has 1 mismatch and 4 no-majority values; 10% has 113 mismatches and 175 no-majority values. No claim of fault-proof decoding. JS scalar binary64 arithmetic is much slower than native. Full bitstream parsing, codec conformance, real playback and DRM are not …

## Licence
PolyForm Small Business 1.0.0 OR PolyForm Noncommercial 1.0.0 + Free Trial; CC BY 4.0 documents. Details: `LICENSING.md`.

## Guided setup through MCP
The SovRing umbrella ships an MCP server that explains, checks readiness, clones and runs these tests: https://github.com/SovRing-Labs/SovRing (README → "Paste this to your assistant"). Zero-install docs endpoint: https://gitmcp.io/SovRing-Labs/armilla-decoder
