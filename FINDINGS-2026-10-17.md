# armilla-decoder — findings addendum, 2026-10-17

This dated addition preserves the earlier disclosure and its historical limits. Code remains pre-alpha. No publication, browser/renderer/playback wiring or production service is established.

## audio

Implementation: 117feature phase-aware frontend; actual modern Hopfield M.T softmax(8Mx),3 iterations; bounded integer IFS; ARS0/RPK0 reconstructed causal history.

Engineering: 62396 worker engineering checks; not parent rerun.

Science: Hopfield top1 30.54%, fallback44.38%; combined Hopfield+IFS slower than exhaustive; model-once16.2396 bits/sample vs raw15.3561 bits/sample, no compression win; quietq16 fidelity failure.

Deployment: research prototype only; recorded fixtures/model not shipped.

Limits: 94–95% historical sieve agreement is not Hopfield accuracy; no90%/4K/codec-conformance claim; q1exact/q16maxerror8 does not imply perceptual fidelity.

## wave

Implementation: Replay/prepare, explicit Pilot/match, complementary read/admit/checked_readback and Affine/diagnose.

Engineering: 2018 parent checks; known collisions/duplicated mirrors/ineffective actuators/stale inputs.

Science: phase/direct complex tie32.26%top1,40.65%recall4 at extra cost; checksum/readback cheaper; integer IFS uncertified.

Deployment: opt-in Python, standalone synthetic video/JS fixtures, actual audio read-only adapter external.

Limits: pilot serialized/costed; only supplied common rotation invariance; no pilot-free/absolute phase or strict integer contraction.

## Run the portable subset

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 findings/smoke.py

This source-only capsule uses tiny synthetic public inputs. See findings/LINEAGE.json for exact own-source hashes and extraction/relocation changes. It does not run the full original5482/2018/1364/62396 gates. Raw recorded audio, trained model, private fixture paths and embedded secret keys are excluded. Full research gates are EXTERNAL to this release. No dependency is installed by the check. NumPy is required for crypto/decoder; defense uses Linux/POSIX flock/fsync. The capsule signs only ephemeral test material in memory and stores no private signing key. Trusted local state is an explicit defense assumption.

Original APIs are described above; some fleet-only convenience entrypoints are deliberately omitted. The smoke exercises the exported subset. Hybrid cost rows are not transferred: unmatched measurements default to positional. Positive historical numbers are descriptive local Python comparisons, not newcomer performance guarantees.
## Minimal caller example

From the released repository root, with NumPy already installed:

```python
import sys
sys.path.insert(0, "findings")
from phase import Pilot, match
query = [1+2j, 3-1j, -2+4j, 5+3j]
pilot = Pilot("synthetic", "block")
rotated_pilot = Pilot("synthetic", "block", real=0., imag=1.)
result = match(query, [[1j*x for x in query]], pilot, [rotated_pilot], common_phase=True)
assert result["scores"][0] < 1e-20
assert result["reference_bytes"] > 0
```

Pilot.serialize/deserialize, match, direct_complex and magnitude are exported own APIs. Complex payload layout is channels-by-samples; mono vectors are accepted. Explicit pilots bind reference/timebase/sample rate/origin/channel count and a supplied unit complex reference. The search considers only zero-padded offsets -2 through2. Common-phase invariance removes one declared supplied rotation across all channels; it preserves relative stereo/temporal differences. Missing/mismatched pilots and silence abstain. No pilot-free/absolute phase benefit is claimed, and retrieval proposals never replace exact payload bytes.

The original audio research uses shared RECONSTRUCTED causal receiver history, actual serialized prediction/residual bytes, q1 exactness and q16 maximum sample error8. The quiet q16 fidelity failure remains; a sample-error bound does not establish perceptual fidelity. This capsule ships no audio recording, frozen model, codec transport or browser/renderer change. Full acoustic comparisons remain EXTERNAL, and the earlier94–95% sieve/exhaustive agreement must not be read as Hopfield accuracy.
