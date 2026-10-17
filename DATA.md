# Data and exclusions

Only curated workload code and small results ship. Compiled binaries, target directories, caches, shared libraries and object files do not ship; rebuild with the suite commands. No models, private datasets or unrelated source trees are included. Files above 1048576 bytes are excluded by the manifest. External math-core index/specs and renderer integration data are not bundled; obtain those projects separately. The VP9 JIT probe depends on a separate renderer emitter and is excluded; use the renderer project for that experiment.

## Selected excluded inputs

| Source-relative path | Bytes | SHA-256 prefix | Reproduction |
|---|---|---|---|
| `jsnum/fast_probe` | 16256 | `1e00c43d26f50591` | Compile the matching C source via suite command. |
| `media/vp9-idct-rings/src/bin/armilla_jit.rs` | 30124 | `8973e1d51d0363af` | Separate renderer emitter prerequisite; use renderer project. |
