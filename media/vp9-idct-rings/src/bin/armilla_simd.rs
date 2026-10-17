//! SIMD (J 10-15 "implement all 3", step 2): the ring arithmetic batched 16 blocks at a time in lane arrays (structure of arrays per
//! ring); x mod m by a float-reciprocal estimate + exact integer correction (no divide), so the butterflies vectorise (AVX2,
//! target-cpu=native). Exactness checked against plain i32 on the whole frame; faults injected between the passes as before.
//! Base: SIMPLEX CORE + CYCLONE COUNTING (J 10-15, readbacks confirmed): the centre JIT math core is a regular simplex at the grid's highest
//! dimension (117 systems = vertices → a 116-simplex); a superblock goes to its nearest vertex. The 20 outer T50 moduli 131…233 count
//! traffic on every link (router→group→diamond→system, 130 link ids) → live bottleneck map per level; the router steers by it.
//! Router variants: direct (cost table + counters), simplex sift (ternary query), random ternary keys (control).
//! Previous: ROUTER (J 10-15 "yes add the router"): tile column → group; block class (DC-only / sparse rows / dense) → path; within the group the
//! diamond with the lowest measured cost + current load. Compared against the fixed band split. Block mix is SYNTHETIC.
//! Base (from armilla_cores.rs): J 10-15 (readbacks confirmed): choose which coprimes sit on which rings, with a ratio per diamond so each diamond computes its share by
//! itself; math cores below the final one — layout A: one per DIAMOND (9) + final = 10 cores; layout B: one per GROUP (3) + final = 4.
//! Grid as before: 3 groups × 3 diamonds × 13 systems (layers 1-3-5-3-1). Workload: VP9 8×8 2D inverse DCT of a 1080p luma frame
//! (spec v0.6 §8.7.1–8.7.2).
//!
//! Ring placement (the per-diamond ratio): every system of a diamond carries that diamond's moduli (data rings + check ring) and its own
//! 2^64 low-bit ring. Diamond 0 = the routing primes (data 127·113·43, check 29: 2^14 ≡ ±1, Round2 is a pure route); diamonds 1–8 = 2 data
//! primes + 1 check prime from 1031…4093 (disjoint), product of data > 2^20 > the 17-bit output range. Inside a system every Round2(x,14) =
//! low bits from the 2^64 ring + exact division on the residues (no bridge). Sub-core (diamond or group math core): two independent
//! copies of each output — residue decode (CRT of the data rings) and the 2^64 ring — plus the check ring; agree → final Round2(x,5);
//! disagree → the copy the check ring confirms wins (in-place correction), none → recompute the block. It sends the block + a checksum
//! (Σ outputs mod 65521) to the FINAL core, which verifies the checksum and assembles the frame.
use std::sync::mpsc;
use std::time::Instant;

const COS64: [i64; 33] = [16384, 16364, 16305, 16207, 16069, 15893, 15679, 15426, 15137, 14811, 14449, 14053, 13623, 13160, 12665, 12140,
    11585, 11003, 10394, 9760, 9102, 8423, 7723, 7005, 6270, 5520, 4756, 3981, 3196, 2404, 1606, 804, 0];
fn cos64(a: i64) -> i64 { let a2 = a & 127; match a2 { 0..=32 => COS64[a2 as usize], 33..=64 => -COS64[(64 - a2) as usize], 65..=96 => -COS64[(a2 - 64) as usize], _ => COS64[(128 - a2) as usize] } }
fn sin64(a: i64) -> i64 { cos64(a - 32) }
fn brev(bits: u32, x: usize) -> usize { let mut t = 0; for i in 0..bits { t |= ((x >> i) & 1) << (bits - 1 - i); } t }
#[derive(Clone, Copy)] enum Op { B(usize, usize, i64, bool), H(usize, usize) }
fn idct_ops(n: usize, ops: &mut Vec<Op>) {
    let (n0, n1, n2) = (1usize << n, 1usize << (n - 1), 1usize << (n - 2)); let n3 = if n >= 3 { 1usize << (n - 3) } else { 0 };
    if n == 2 { ops.push(Op::B(0, 1, 16, true)); } else { idct_ops(n - 1, ops); }
    for i in 0..n2 { ops.push(Op::B(n1 + i, n0 - 1 - i, 32 - brev(5, n1 + i) as i64, false)); }
    if n >= 3 { for i in 0..n3 { for j in 0..2 { let (a, b) = (n1 + 4 * i + 2 * j, n1 + 1 + 4 * i + 2 * j); ops.push(if j == 1 { Op::H(b, a) } else { Op::H(a, b) }); } } }
    if n >= 3 { for i in 0..n3 { ops.push(Op::B(n0 - n3 - 1 - i, n1 + n3 + i, 16, true)); } }
    for i in 0..n1 { ops.push(Op::H(i, n0 - 1 - i)); }
}
fn round2(x: i64, n: u32) -> i64 { (x + (1 << (n - 1))) >> n }
fn inv(a: u64, m: u64) -> u64 { let (mut t, mut nt, mut r, mut nr) = (0i64, 1i64, m as i64, (a % m) as i64); while nr != 0 { let q = r / nr; (t, nt) = (nt, t - q * nt); (r, nr) = (nr, r - q * nr); } ((t % m as i64 + m as i64) % m as i64) as u64 }
fn plain_block(blk: &[i32; 64], ops: &[Op]) -> [i32; 64] {
    let mut d = *blk;
    let one = |t: &mut [i64; 8]| { let c = *t; for i in 0..8 { t[i] = c[brev(3, i)]; }
        for &op in ops { match op { Op::B(a, b, ang, f) => { let (x, y) = (t[a] * cos64(ang) - t[b] * sin64(ang), t[a] * sin64(ang) + t[b] * cos64(ang)); t[a] = round2(x, 14); t[b] = round2(y, 14); if f { t.swap(a, b); } }
                                     Op::H(a, b) => { let (x, y) = (t[a], t[b]); t[a] = x + y; t[b] = x - y; } } } };
    for i in 0..8 { let mut t = [0i64; 8]; for j in 0..8 { t[j] = d[i * 8 + j] as i64; } one(&mut t); for j in 0..8 { d[i * 8 + j] = t[j] as i32; } }
    for j in 0..8 { let mut t = [0i64; 8]; for i in 0..8 { t[i] = d[i * 8 + j] as i64; } one(&mut t); for i in 0..8 { d[i * 8 + j] = round2(t[i], 5) as i32; } }
    d
}

const KM: usize = 4; // max rings per diamond (data + check)
#[derive(Clone)]
struct Diamond { p: Vec<u64>, data: usize, inv14: Vec<u64>, crt_w: Vec<u128>, data_prod: u128 }
impl Diamond {
    fn new(p: Vec<u64>, data: usize) -> Diamond {
        let inv14 = p.iter().map(|&m| inv(16384 % m, m)).collect();
        let data_prod: u128 = p[..data].iter().map(|&m| m as u128).product();
        let crt_w = p[..data].iter().map(|&m| { let mi = data_prod / m as u128; mi * inv((mi % m as u128) as u64, m) as u128 % data_prod }).collect();
        Diamond { p, data, inv14, crt_w, data_prod }
    }
    fn n(&self) -> usize { self.p.len() }
    fn decode(&self, r: &[u64; KM]) -> i64 {
        let mut x: u128 = 0; for k in 0..self.data { x = (x + r[k] as u128 * self.crt_w[k]) % self.data_prod; }
        if x > self.data_prod / 2 { x as i128 as i64 - self.data_prod as i64 } else { x as i64 }
    }
    /// one 1D inverse DCT on this diamond's rings + the system's 2^64 ring
    fn idct(&self, r: &mut [[u64; KM]; 8], w: &mut [u64; 8], ops: &[Op]) {
        let n = self.n();
        let (cr, cw) = (*r, *w); for i in 0..8 { r[i] = cr[brev(3, i)]; w[i] = cw[brev(3, i)]; }
        for &op in ops { match op {
            Op::B(a, b, ang, flip) => {
                let (c, s) = (cos64(ang), sin64(ang));
                let xw = w[a].wrapping_mul(c as u64).wrapping_add(w[b].wrapping_mul((-s) as u64)).wrapping_add(8192);
                let yw = w[a].wrapping_mul(s as u64).wrapping_add(w[b].wrapping_mul(c as u64)).wrapping_add(8192);
                let (lx, ly) = (xw & 0x3FFF, yw & 0x3FFF);
                for k in 0..n {
                    let m = self.p[k];
                    let (cm, sm, nsm) = ((c.rem_euclid(m as i64)) as u64, (s.rem_euclid(m as i64)) as u64, ((-s).rem_euclid(m as i64)) as u64);
                    let x = (r[a][k] * cm + r[b][k] * nsm + 8192 + m * 16384 - lx) % m;
                    let y = (r[a][k] * sm + r[b][k] * cm + 8192 + m * 16384 - ly) % m;
                    r[a][k] = x * self.inv14[k] % m; r[b][k] = y * self.inv14[k] % m;
                }
                w[a] = ((xw as i64) >> 14) as u64; w[b] = ((yw as i64) >> 14) as u64;
                if flip { r.swap(a, b); w.swap(a, b); }
            }
            Op::H(a, b) => {
                for k in 0..n { let m = self.p[k]; let (x, y) = (r[a][k], r[b][k]); r[a][k] = (x + y) % m; r[b][k] = (x + m - y) % m; }
                let (x, y) = (w[a], w[b]); w[a] = x.wrapping_add(y); w[b] = x.wrapping_sub(y);
            }
        } }
    }
    /// a system computes one block (the share): output residues + 2^64 values before the final Round2. fault: (pos, which: ring k or KM = carrier bit)
    fn block(&self, blk: &[i32; 64], ops: &[Op], fault: Option<(usize, usize)>, kind: u8) -> (Vec<[u64; KM]>, [i64; 64]) {
        let mut rr = vec![[0u64; KM]; 64]; let mut ww = [0u64; 64];
        for i in 0..64 { let x = blk[i] as i64; for k in 0..self.n() { rr[i][k] = x.rem_euclid(self.p[k] as i64) as u64; } ww[i] = x as u64; }
        // row pass: kind 1/2 skip all-zero rows (the inverse DCT of a zero vector is zero: B and H of zeros stay zero, Round2(0) = 0)
        for i in 0..8 { if kind > 0 && (0..8).all(|j| blk[i * 8 + j] == 0) { continue; }
            let mut r = [[0u64; KM]; 8]; let mut w = [0u64; 8]; for j in 0..8 { r[j] = rr[i * 8 + j]; w[j] = ww[i * 8 + j]; }
            self.idct(&mut r, &mut w, ops); for j in 0..8 { rr[i * 8 + j] = r[j]; ww[i * 8 + j] = w[j]; } }
        if let Some((pos, which)) = fault { if which < self.n() { rr[pos][which] = (rr[pos][which] + 1) % self.p[which]; } else { ww[pos] ^= 1 << 7; } }
        // DC-only: every column holds the same input when row 0 came out constant (and no fault touched it) → one column, replicated
        let same_cols = kind == 2 && (1..8).all(|j| ww[j] == ww[0] && rr[j] == rr[0]) && (8..64).all(|i| ww[i] == 0 && rr[i][..self.n()].iter().all(|&x| x == 0));
        for j in 0..8 {
            if same_cols && j > 0 { for i in 0..8 { rr[i * 8 + j] = rr[i * 8]; ww[i * 8 + j] = ww[i * 8]; } continue; }
            let mut r = [[0u64; KM]; 8]; let mut w = [0u64; 8]; for i in 0..8 { r[i] = rr[i * 8 + j]; w[i] = ww[i * 8 + j]; }
            self.idct(&mut r, &mut w, ops); for i in 0..8 { rr[i * 8 + j] = r[i]; ww[i * 8 + j] = w[i]; } }
        let mut wi = [0i64; 64]; for i in 0..64 { wi[i] = ww[i] as i64; }
        (rr, wi)
    }
}

/// a sub-core (diamond or group math core): check, correct in place, final Round2(x,5); returns (block, how: 0 clean 1 corrected 2 recomputed)
fn sub_core(d: &Diamond, rr: &[[u64; KM]], w: &[i64; 64], blk: &[i32; 64], ops: &[Op]) -> ([i32; 64], u8) {
    let mut out = [0i32; 64]; let mut how = 0u8;
    for i in 0..64 {
        // fast path: the 2^64 copy agrees with every ring → it is the value (no CRT); decode only on a mismatch
        let vw = w[i];
        if (0..d.n()).all(|k| vw.rem_euclid(d.p[k] as i64) as u64 == rr[i][k]) { out[i] = round2(vw, 5) as i32; continue; }
        let vr = d.decode(&rr[i]);
        let ok_r = (d.data..d.n()).all(|k| vr.rem_euclid(d.p[k] as i64) as u64 == rr[i][k]);   // residue copy consistent with its check ring
        // carrier confirmed when it agrees with all rings but at most one (that ring is the faulty one)
        let ok_w = (0..d.n()).filter(|&k| vw.rem_euclid(d.p[k] as i64) as u64 != rr[i][k]).count() <= 1;
        let v = if vr == vw && ok_r { vr } else if ok_w { how = how.max(1); vw } else if ok_r { how = how.max(1); vr } else { return (plain_block(blk, ops), 2) };
        out[i] = round2(v, 5) as i32;
    }
    (out, how)
}
fn checksum(o: &[i32; 64]) -> u64 { o.iter().fold(0u64, |a, &x| (a + x.rem_euclid(65521) as u64) % 65521) }


const L: usize = 16;
#[inline(always)] fn modm(x: u32, m: u32, mi: f32) -> u32 {
    // x < 2^27 fits i32: signed converts vectorise on AVX2 (vcvtdq2ps / vcvttps2dq); there is no unsigned form
    let q = ((x as i32) as f32 * mi) as i32 as u32; let mut r = x.wrapping_sub(q.wrapping_mul(m)) as i32;
    if r < 0 { r += m as i32; } if r < 0 { r += m as i32; } if r >= m as i32 { r -= m as i32; } if r >= m as i32 { r -= m as i32; }
    r as u32
}
impl Diamond {
    /// one 1D inverse DCT on L lanes: t[pos][ring][lane] residues, w[pos][lane] the 2^64 rings
    fn idct_lanes(&self, t: &mut [[[u32; L]; KM]; 8], w: &mut [[u64; L]; 8], ops: &[Op]) {
        let n = self.n();
        let (ct, cw) = (*t, *w); for i in 0..8 { t[i] = ct[brev(3, i)]; w[i] = cw[brev(3, i)]; }
        let mi: Vec<f32> = self.p.iter().map(|&m| 1.0 / m as f32).collect();
        for &op in ops { match op {
            Op::B(a, b, ang, flip) => {
                let (c, s) = (cos64(ang), sin64(ang));
                let mut lx = [0u32; L]; let mut ly = [0u32; L];
                for l in 0..L {
                    let xw = w[a][l].wrapping_mul(c as u64).wrapping_add(w[b][l].wrapping_mul((-s) as u64)).wrapping_add(8192);
                    let yw = w[a][l].wrapping_mul(s as u64).wrapping_add(w[b][l].wrapping_mul(c as u64)).wrapping_add(8192);
                    lx[l] = (xw & 0x3FFF) as u32; ly[l] = (yw & 0x3FFF) as u32;
                    w[a][l] = ((xw as i64) >> 14) as u64; w[b][l] = ((yw as i64) >> 14) as u64;
                }
                for k in 0..n {
                    let m = self.p[k] as u32; let f = mi[k];
                    let (cm, sm, nsm) = (c.rem_euclid(m as i64) as u32, s.rem_euclid(m as i64) as u32, (-s).rem_euclid(m as i64) as u32);
                    let iv = self.inv14[k] as u32; let bias = 8192 + m * 16384;
                    let (ra, rb) = (t[a][k], t[b][k]);
                    let (mut ox, mut oy) = ([0u32; L], [0u32; L]);
                    for l in 0..L {
                        let x = modm(ra[l] * cm + rb[l] * nsm + bias - lx[l], m, f);
                        let y = modm(ra[l] * sm + rb[l] * cm + bias - ly[l], m, f);
                        ox[l] = modm(x * iv, m, f); oy[l] = modm(y * iv, m, f);
                    }
                    t[a][k] = ox; t[b][k] = oy;
                }
                if flip { t.swap(a, b); w.swap(a, b); }
            }
            Op::H(a, b) => {
                for k in 0..n { let m = self.p[k] as u32; let (ra, rb) = (t[a][k], t[b][k]);
                    for l in 0..L { let x = ra[l] + rb[l]; t[a][k][l] = if x >= m { x - m } else { x }; let y = ra[l] + m - rb[l]; t[b][k][l] = if y >= m { y - m } else { y }; } }
                for l in 0..L { let (x, y) = (w[a][l], w[b][l]); w[a][l] = x.wrapping_add(y); w[b][l] = x.wrapping_sub(y); }
            }
        } }
    }
    /// L blocks at once (full 2D transform), the same results as `block(.., kind 0)` per block
    fn block_lanes(&self, blks: &[&[i32; 64]], ops: &[Op], faults: &[Option<(usize, usize)>]) -> Vec<(Vec<[u64; KM]>, [i64; 64])> {
        let nb = blks.len(); let n = self.n();
        let mut rr = vec![[[0u32; L]; KM]; 64]; let mut ww = vec![[0u64; L]; 64];
        for (l, b) in blks.iter().enumerate() { for i in 0..64 { let x = b[i] as i64; for k in 0..n { rr[i][k][l] = x.rem_euclid(self.p[k] as i64) as u32; } ww[i][l] = x as u64; } }
        for i in 0..8 { let mut t = [[[0u32; L]; KM]; 8]; let mut w = [[0u64; L]; 8]; for j in 0..8 { t[j] = rr[i * 8 + j]; w[j] = ww[i * 8 + j]; }
            self.idct_lanes(&mut t, &mut w, ops); for j in 0..8 { rr[i * 8 + j] = t[j]; ww[i * 8 + j] = w[j]; } }
        for (l, f) in faults.iter().enumerate().take(nb) { if let Some((pos, which)) = f { if *which < n { rr[*pos][*which][l] = (rr[*pos][*which][l] + 1) % self.p[*which] as u32; } else { ww[*pos][l] ^= 1 << 7; } } }
        for j in 0..8 { let mut t = [[[0u32; L]; KM]; 8]; let mut w = [[0u64; L]; 8]; for i in 0..8 { t[i] = rr[i * 8 + j]; w[i] = ww[i * 8 + j]; }
            self.idct_lanes(&mut t, &mut w, ops); for i in 0..8 { rr[i * 8 + j] = t[i]; ww[i * 8 + j] = w[i]; } }
        (0..nb).map(|l| { let mut o = vec![[0u64; KM]; 64]; let mut w = [0i64; 64]; for i in 0..64 { for k in 0..n { o[i][k] = rr[i][k][l] as u64; } w[i] = ww[i][l] as i64; } (o, w) }).collect()
    }
}

fn class_of(b: &[i32; 64]) -> u8 { if b[1..].iter().all(|&x| x == 0) { 2 } else if (0..8).filter(|&i| (0..8).any(|j| b[i * 8 + j] != 0)).count() <= 3 { 1 } else { 0 } }

const SYS: usize = 117;
const OUTER: [u64; 20] = [131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233];
/// cyclone counting: one counter row per outer modulus, indexed by link id mod m (count-min: the estimate is the min over rows)
struct Outer { c: Vec<Vec<f64>> }
impl Outer {
    fn new() -> Outer { Outer { c: OUTER.iter().map(|&m| vec![0.0; m as usize]).collect() } }
    fn add(&mut self, link: usize, v: f64) { for (k, &m) in OUTER.iter().enumerate() { self.c[k][link % m as usize] += v; } }
    fn est(&self, link: usize) -> f64 { OUTER.iter().enumerate().map(|(k, &m)| self.c[k][link % m as usize]).fold(f64::MAX, f64::min) }
}
// link ids: 0..3 router→group g, 3..12 group→diamond d, 12..129 diamond→system s
fn link_g(g: usize) -> usize { g } fn link_d(d: usize) -> usize { 3 + d } fn link_s(s: usize) -> usize { 12 + s }

/// regular simplex keys: vertex i = e_i − centroid in R^SYS (rank SYS−1), unit length; <v_i, v_j> = −1/(SYS−1) for i ≠ j
fn simplex_keys() -> Vec<Vec<f64>> {
    let c = 1.0 / SYS as f64; let nrm = ((1.0 - c) * (1.0 - c) + (SYS as f64 - 1.0) * c * c).sqrt();
    (0..SYS).map(|i| (0..SYS).map(|j| (if i == j { 1.0 - c } else { -c }) / nrm).collect()).collect()
}
fn random_ternary_keys(seed: u64) -> Vec<Vec<f64>> {
    let mut s = seed; let mut nx = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    (0..SYS).map(|_| (0..SYS).map(|_| match nx() % 3 { 0 => -1.0, 1 => 0.0, _ => 1.0 }).collect()).collect()
}
/// ternary query: the superposition Σ desirability_i · key_i, each coordinate quantised to {−1, 0, +1} by a threshold at ±⅓ of the max
fn ternary(q: &[f64]) -> Vec<f64> { let mx = q.iter().fold(0f64, |a, &x| a.max(x.abs())); q.iter().map(|&x| if x > mx / 3.0 { 1.0 } else if x < -mx / 3.0 { -1.0 } else { 0.0 }).collect() }

fn run(mode: &str, diamonds: &[Diamond], blocks: &[[i32; 64]], faults: &[Option<(usize, usize)>], plain: &[[i32; 64]], ops: &[Op], bw: usize, bh: usize, cost: &[[f64; 3]], intra: &[bool], keys: &[Vec<f64>]) -> (String, Vec<usize>) {
    let n = blocks.len();
    let (sbw, sbh) = ((bw + 7) / 8, (bh + 7) / 8);
    let mut plan: Vec<Vec<(usize, u8)>> = vec![vec![]; 9];
    let mut outer = Outer::new();
    let mut choice = vec![0usize; sbw * sbh];
    for sb in 0..sbw * sbh {
        let (sx, sy) = (sb % sbw, sb / sbw);
        let members: Vec<usize> = (0..8).flat_map(|dy| (0..8).map(move |dx| (sy * 8 + dy, sx * 8 + dx))).filter(|&(y, x)| y < bh && x < bw).map(|(y, x)| y * bw + x).collect();
        let tile = sx * 3 / sbw;
        // this superblock's cost on each diamond (sum of its blocks' class costs)
        let ucost: Vec<f64> = (0..9).map(|d| members.iter().map(|&bi| cost[d][class_of(&blocks[bi]) as usize]).sum()).collect();
        // desirability of system s: −(unit cost on its diamond + traffic estimate on its links, from the outer rings); intra → own tile's group only
        let des: Vec<f64> = (0..SYS).map(|s| { let d = s / 13; let g = d / 3;
            if intra[sb] && g != tile { return f64::NEG_INFINITY; }
            -(ucost[d] + outer.est(link_s(s)) * 13.0 + outer.est(link_d(d)) + outer.est(link_g(g)) / 3.0) }).collect();
        let s_pick = match mode {
            "fixed" => { let g = tile; let d = g * 3 + (sy * 3 / sbh).min(2); d * 13 + sb % 13 }
            "direct" => (0..SYS).max_by(|&a, &b| des[a].partial_cmp(&des[b]).unwrap()).unwrap(),
            _ => {
                // simplex / random keys: query = Σ (des_i − min) · key_i over eligible, quantised ternary; sift = argmax <q, key_j> over eligible
                let fin: Vec<f64> = des.iter().cloned().filter(|x| x.is_finite()).collect();
                let lo = fin.iter().cloned().fold(f64::MAX, f64::min); let hi = fin.iter().cloned().fold(f64::MIN, f64::max);
                let w: Vec<f64> = des.iter().map(|&x| if x.is_finite() && hi > lo { ((x - lo) / (hi - lo)).powi(8) } else { 0.0 }).collect();
                let mut q = vec![0f64; SYS]; for i in 0..SYS { if w[i] > 0.0 { for j in 0..SYS { q[j] += w[i] * keys[i][j]; } } }
                let q = ternary(&q);
                (0..SYS).filter(|&s| des[s].is_finite()).max_by(|&a, &b| { let sa: f64 = (0..SYS).map(|j| q[j] * keys[a][j]).sum(); let sb2: f64 = (0..SYS).map(|j| q[j] * keys[b][j]).sum(); sa.partial_cmp(&sb2).unwrap() }).unwrap()
            }
        };
        choice[sb] = s_pick;
        let d = s_pick / 13; let g = d / 3;
        outer.add(link_g(g), ucost[d]); outer.add(link_d(d), ucost[d]); outer.add(link_s(s_pick), ucost[d]);
        for &bi in &members { let kind = if mode == "fixed" { class_of(&blocks[bi]) } else { class_of(&blocks[bi]) }; plan[d].push((bi, kind)); }
    }
    // bottleneck map from the outer rings: utilisation per level = hottest link / mean link of that level
    let lvl = |r: std::ops::Range<usize>| { let v: Vec<f64> = r.map(|l| outer.est(l)).collect(); let mean = v.iter().sum::<f64>() / v.len() as f64; v.iter().cloned().fold(0f64, f64::max) / mean.max(1e-12) };
    let (ug, ud) = (lvl(0..3), lvl(3..12));
    let inline = std::env::var_os("CORE_INLINE").is_some(); // layout A: the math core runs inside each diamond (no hop to a group core)
    let t1 = Instant::now();
    let (final_tx, final_rx) = mpsc::channel::<(usize, [i32; 64], u64, u8)>();
    let (frame, corrected, recomputed, bad_sum, busy) = std::thread::scope(|sc| {
        let mut core_tx = vec![]; let mut core_hs = vec![];
        for _c in 0..3 { let (tx, rx) = mpsc::channel::<(usize, usize, Vec<[u64; KM]>, [i64; 64])>(); core_tx.push(tx); let ftx = final_tx.clone();
            core_hs.push(sc.spawn(move || { for (d, bi, rr, w) in rx { let (o, how) = sub_core(&diamonds[d], &rr, &w, &blocks[bi], ops); ftx.send((bi, o, checksum(&o), how)).unwrap(); } })); }
        let dh: Vec<_> = (0..9).map(|d| { let (plan_d, txs, ftx_d) = (plan[d].clone(), core_tx.clone(), final_tx.clone());
            sc.spawn(move || { let t = Instant::now();
                let mut batch: Vec<usize> = Vec::with_capacity(L);
                let flush = |batch: &mut Vec<usize>| { if batch.is_empty() { return; }
                    let bl: Vec<&[i32; 64]> = batch.iter().map(|&bi| &blocks[bi]).collect(); let fl: Vec<Option<(usize, usize)>> = batch.iter().map(|&bi| faults[bi]).collect();
                    for (k, (rr, w)) in diamonds[d].block_lanes(&bl, ops, &fl).into_iter().enumerate() {
                        if inline { let bi = batch[k]; let (o, how) = sub_core(&diamonds[d], &rr, &w, &blocks[bi], ops); ftx_d.send((bi, o, checksum(&o), how)).unwrap(); }
                        else { txs[d / 3].send((d, batch[k], rr, w)).unwrap(); } }
                    batch.clear(); };
                for (bi, kind) in plan_d {
                    if kind == 2 { let (rr, w) = diamonds[d].block(&blocks[bi], ops, faults[bi], kind);
                        if inline { let (o, how) = sub_core(&diamonds[d], &rr, &w, &blocks[bi], ops); ftx_d.send((bi, o, checksum(&o), how)).unwrap(); } else { txs[d / 3].send((d, bi, rr, w)).unwrap(); } continue; }
                    batch.push(bi); if batch.len() == L { flush(&mut batch); }
                }
                flush(&mut batch);
                t.elapsed().as_secs_f64() }) }).collect();
        drop(core_tx); drop(final_tx);
        let mut frame = vec![[0i32; 64]; n]; let (mut c1, mut c2, mut bad) = (0, 0, 0);
        for (bi, o, cs, how) in final_rx { if checksum(&o) != cs { bad += 1; } frame[bi] = o; if how == 1 { c1 += 1; } if how == 2 { c2 += 1; } }
        let busy: Vec<f64> = dh.into_iter().map(|h| h.join().unwrap()).collect();
        for h in core_hs { h.join().unwrap(); }
        (frame, c1, c2, bad, busy)
    });
    let ta = t1.elapsed().as_secs_f64();
    let wrong = frame.iter().zip(plain).filter(|(a, b)| a != b).count();
    let gbusy: Vec<String> = (0..3).map(|g| format!("{:.4}", busy[g * 3..g * 3 + 3].iter().cloned().fold(0., f64::max))).collect();
    (format!("{{\"mode\":\"{mode}\",\"frames_per_s\":{:.2},\"wall_s\":{ta:.4},\"blocks_wrong\":{wrong},\"corrected_in_place\":{corrected},\"recomputed\":{recomputed},\"final_checksum_failures\":{bad_sum},\"group_busy_s\":[{}],\"outer_ring_bottleneck_group_level\":{ug:.2},\"outer_ring_bottleneck_diamond_level\":{ud:.2}}}", 1.0 / ta, gbusy.join(",")), choice)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let fault_rate: f64 = a.get(1).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let mut ops = vec![]; idct_ops(3, &mut ops);
    let isp = |n: u64| n > 1 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0);
    let pool: Vec<u64> = (1031u64..4096).filter(|&n| isp(n)).collect();
    let mut diamonds = vec![Diamond::new(vec![127, 113, 43, 29], 3)];
    for k in 1..9 { let j = (k - 1) * 3; diamonds.push(Diamond::new(vec![pool[j * 7], pool[j * 7 + 50], pool[j * 7 + 100]], 2)); }
    let (bw, bh) = (240usize, 135usize);
    let mut s: u64 = 0x9E3779B97F4A7C15;
    let mut rnd = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    let blocks: Vec<[i32; 64]> = (0..bw * bh).map(|bi| { let x = bi % bw; let busy = x * 3 / bw == 2; let r = rnd() % 100;
        let class = if busy { if r < 15 { 2 } else if r < 40 { 1 } else { 0 } } else if r < 50 { 2 } else if r < 88 { 1 } else { 0 };
        let mut b = [0i32; 64]; b[0] = (rnd() % 1024) as i32 - 512;
        if class == 1 { for _ in 0..(1 + rnd() % 5) { let i = (rnd() % 3) as usize; let j = (rnd() % 4) as usize; b[i * 8 + j] = (rnd() % 256) as i32 - 128; } }
        if class == 0 { for (k, c) in b.iter_mut().enumerate() { if k < 20 || rnd() % 4 == 0 { *c = (rnd() % 1024) as i32 - 512; } } }
        b }).collect();
    let plain: Vec<[i32; 64]> = blocks.iter().map(|b| plain_block(b, &ops)).collect();
    let (sbw, sbh) = ((bw + 7) / 8, (bh + 7) / 8);
    let intra: Vec<bool> = (0..sbw * sbh).map(|_| rnd() % 100 < 30).collect(); // SYNTHETIC: 30 % intra superblocks pinned to their tile's group
    let mut faults: Vec<Option<(usize, usize)>> = vec![None; bw * bh];
    if fault_rate > 0.0 { for f in faults.iter_mut() { if (rnd() % 1_000_000) as f64 / 1e6 < fault_rate { *f = Some(((rnd() % 64) as usize, (rnd() % 4) as usize)); } } }
    let mut cost = vec![[0f64; 3]; 9];
    for d in 0..9 { for kind in 0..3u8 { let sample: Vec<&[i32; 64]> = blocks.iter().filter(|b| class_of(b) == kind).take(200).collect();
        let t = Instant::now(); for b in &sample { let _ = diamonds[d].block(b, &ops, None, kind); } cost[d][kind as usize] = t.elapsed().as_secs_f64() / sample.len().max(1) as f64; } }
    let sk = simplex_keys(); let rk = random_ternary_keys(7);
    // geometry check: simplex keys pairwise inner products (should all be −1/116)
    let ip = (0..SYS).flat_map(|i| (0..SYS).filter(move |&j| j != i).map(move |j| (i, j))).map(|(i, j)| (0..SYS).map(|k| sk[i][k] * sk[j][k]).sum::<f64>()).fold((f64::MAX, f64::MIN), |a, x| (a.0.min(x), a.1.max(x)));
    println!("{{\"systems_vertices\":{SYS},\"simplex_dimension\":{},\"simplex_inner_product_min_max\":[{:.5},{:.5}],\"expected\":{:.5},\"faults_injected\":{}}}", SYS - 1, ip.0, ip.1, -1.0 / (SYS as f64 - 1.0), faults.iter().filter(|f| f.is_some()).count());
    let mut choices = vec![];
    for (mode, keys) in [("fixed", &sk), ("direct", &sk), ("simplex", &sk), ("random-keys", &rk)] {
        let (line, ch) = run(mode, &diamonds, &blocks, &faults, &plain, &ops, bw, bh, &cost, &intra, keys); println!("{line}"); choices.push((mode, ch));
    }
    let agree = |a: &Vec<usize>, b: &Vec<usize>| a.iter().zip(b).filter(|(x, y)| x == y).count() as f64 / a.len() as f64;
    let agree_d = |a: &Vec<usize>, b: &Vec<usize>| a.iter().zip(b).filter(|(x, y)| *x / 13 == *y / 13).count() as f64 / a.len() as f64;
    println!("{{\"agreement_with_direct_system\":{{\"simplex\":{:.3},\"random_keys\":{:.3}}},\"agreement_with_direct_diamond\":{{\"simplex\":{:.3},\"random_keys\":{:.3}}}}}",
        agree(&choices[2].1, &choices[1].1), agree(&choices[3].1, &choices[1].1), agree_d(&choices[2].1, &choices[1].1), agree_d(&choices[3].1, &choices[1].1));
}
