# Provenance

Origin: Armilla workload research, authored by Juan Carlos De Santiago with the Sovereign Raptor fleet (AI-assisted). The fleet assisted implementation, test design and documentation. Publication extracts code and small evidence; manifest rewrites only relocate local paths, with no source-project edits. Third-party implementation code: none knowingly copied; mathematical field sources are credited in TECHNICAL-DISCLOSURE.md. NumPy and system compilers are prerequisites and not bundled.

Common abacus engine: https://github.com/SovRing-Labs/armilla-astrolabe at `120685686f76d0f99e6f9dd9f4c9e038ed1b956e`. Crypto Paillier already ships there as `experiments/paillier/paillier_merge.py` and its result JSON; do not duplicate it here. After obtaining that pinned repo separately, run its own documented experiment and tests. No engine download is performed by this suite.

Renderer wiring is separate: ARMILLA-WIRE-0 reached renderer revision `4b015c8`; its Apache-2.0 armilla.rs engine code is cited, not included. The recorded integration gate compared 3192191 samples, 0 different; this workload export does not rerun that gate. Video colour conversion, audio gain and JS sampling integration are distinct from complete media decode.

Clean-room decoder provenance: H.264 section 8.5.12 and VP9 section 8.7 specification prose, not libvpx/ffmpeg implementation code. Local clean-room scan evidence: results/2026-10-17-cleanroom.log (decoder export). Scanner success checks identifiers; it does not establish legal provenance by itself.


## Findings update — 2026-10-17

See [FINDINGS-2026-10-17.md](FINDINGS-2026-10-17.md) for later implemented opt-in APIs, exact engineering/scientific limits and a portable source-only smoke check. Earlier statements above remain the dated historical record; they are scoped to the earlier bundle. The new capsule is not a full research-gate reproduction or deployed integration.
