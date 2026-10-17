#!/usr/bin/env python3
"""ARMILLA-MEDIA-0 probe: is a real decoder stage exact on the rings? H.264 4x4 inverse core transform (ITU-T H.264 8.5.12: row pass, column pass, (x+32)>>6) on a 1080p luma frame's worth of blocks (129,600 blocks x 16 coefficients).
Reference = plain int64 NumPy of the published formula. Ring version: values held as residues on 12 rings (3 groups x 4 of the largest primes < 256); add/sub are per-ring; the ">>1" and the final ">>6" are the only magnitude steps: read the group by Garner,
vote 2-of-3, shift, re-encode (the same step the reservoir leak does). Fault test: one wrong residue in a random value of a fraction of blocks after every stage; output must be identical. -> idct-rings.json
(Reference formula is from the published standard text as the seat knows it; it is NOT checked against JM or a real bitstream here.)"""
import json, math, sys, time
import numpy as np
P = [p for p in range(3, 256) if all(p % d for d in range(2, int(p ** .5) + 1))][::-1][:12]; G = [P[0:4], P[4:8], P[8:12]]; PV = np.array(P, np.int64)
GM = [math.prod(g) for g in G]; assert all(m > 2 ** 24 for m in GM)
def enc(x): return x[..., None] % PV                                    # int64 [...] -> residues [...,12]
def dec_group(R, g):                                                    # Garner, int64, centred
    ps = G[g]; v = R[..., 4 * g]; pre = ps[0]
    for k in range(1, 4):
        p = ps[k]; d = ((R[..., 4 * g + k] - v % p) * pow(pre % p, -1, p)) % p; v = v + d * pre; pre *= p
    return np.where(v > pre // 2, v - pre, v)
STAT = {"voted": 0, "unrecoverable": 0}
def value(R):                                                           # 2-of-3 vote; no majority = detected (counted); falls back to group 0
    a, b, c = dec_group(R, 0), dec_group(R, 1), dec_group(R, 2); ok = (a == b) | (a == c) | (b == c)
    STAT["voted"] += int(((a != b) | (a != c) | (b != c)).sum()); STAT["unrecoverable"] += int((~ok).sum()); return np.where(a == b, a, np.where(a == c, a, b))
def shr(R, k, off=0):                                                   # magnitude step: read, (x+off)>>k, re-encode
    return enc((value(R) + off) >> k)
def inject(R, p, rng):
    if p <= 0: return R
    R = R.copy(); n = R.shape[0]; sel = np.nonzero(rng.random(n) < p)[0]
    for b in sel:
        i = rng.integers(R.shape[1]); r = rng.integers(12); R[b, i, r] = (R[b, i, r] + rng.integers(1, P[r])) % P[r]
    return R
def idct_ref(d):                                                        # d [B,4,4] int64, rows then columns
    def pass1(x):                                                       # x [...,4] along last axis
        e0 = x[..., 0] + x[..., 2]; e1 = x[..., 0] - x[..., 2]; e2 = (x[..., 1] >> 1) - x[..., 3]; e3 = x[..., 1] + (x[..., 3] >> 1)
        return np.stack([e0 + e3, e1 + e2, e1 - e2, e0 - e3], -1)
    f = pass1(d); g = pass1(f.swapaxes(1, 2)).swapaxes(1, 2); return (g + 32) >> 6
def idct_ring(d, fault=0.0, seed=0):
    rng = np.random.default_rng(seed); R = enc(d.reshape(-1, 16))                 # [B,16,12]
    def add(*terms):
        out = 0
        for s, t in terms: out = out + s * t
        return out % PV
    def pass_(R, rows):                                                 # rows: True = along row index order i*4+j (j varies), False = columns
        ix = (lambda i, j: i * 4 + j) if rows else (lambda i, j: j * 4 + i)
        Rn = np.empty_like(R)
        for i in range(4):
            d0, d1, d2, d3 = (R[:, ix(i, j)] for j in range(4))
            h1, h3 = shr(d1, 1), shr(d3, 1)                             # the two >>1 (magnitude steps)
            e0 = add((1, d0), (1, d2)); e1 = add((1, d0), (-1, d2)); e2 = add((1, h1), (-1, d3)); e3 = add((1, d1), (1, h3))
            Rn[:, ix(i, 0)] = add((1, e0), (1, e3)); Rn[:, ix(i, 1)] = add((1, e1), (1, e2)); Rn[:, ix(i, 2)] = add((1, e1), (-1, e2)); Rn[:, ix(i, 3)] = add((1, e0), (-1, e3))
        return inject(Rn, fault, rng)
    R = pass_(R, True); R = pass_(R, False); out = value(shr(R, 6, 32))
    return out.reshape(-1, 4, 4)
if __name__ == "__main__":
    B = 129600; rng = np.random.default_rng(1); d = rng.integers(-3000, 3001, (B, 4, 4)).astype(np.int64)
    d[rng.random((B, 4, 4)) < 0.7] = 0                                  # typical sparse residuals
    t = time.perf_counter(); ref = idct_ref(d); t_ref = time.perf_counter() - t
    res = {"blocks": B, "frame": "1080p luma 4x4 blocks", "t_plain_numpy_s": round(t_ref, 4)}
    for name, f in (("clean", 0.0), ("fault_1pct_blocks_per_stage", 0.01), ("fault_10pct", 0.10)):
        STAT.update(voted=0, unrecoverable=0); t = time.perf_counter(); got = idct_ring(d, f); dt = time.perf_counter() - t
        res[name] = {"mismatching_samples": int((got != ref).sum()), "t_s": round(dt, 2), "values_repaired_by_vote": STAT["voted"], "no_majority_detected": STAT["unrecoverable"]}
    res["ring_vs_plain_numpy_slowdown"] = round(res["clean"]["t_s"] / t_ref, 1); res["frames_per_s_luma_ring_numpy"] = round(1 / res["clean"]["t_s"], 2)
    json.dump(res, open("idct-rings.json", "w"), indent=1); print(json.dumps(res, indent=1))
