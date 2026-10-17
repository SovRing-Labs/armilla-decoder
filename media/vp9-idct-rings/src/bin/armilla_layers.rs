//! J 10-15 19:40 (readback confirmed, "yes, explore"): each diamond becomes 1-3-7-3-1 ring systems (15), primes placed by SIZE per layer —
//! the 7-system middle layer mainly small primes, the two 3-system layers medium primes, the single top and bottom systems large primes.
//! Grid as before: 3 groups × 3 diamonds → group math cores → final math core (layout B of armilla_cores.rs). Workload: VP9 8×8 2D inverse
//! DCT of a 1080p luma frame (spec v0.6 §8.7.1–8.7.2), checked block by block against plain integers.
//!
//! Reading V (vertical), mode `v137`: one value spans the diamond, each system holding one ring of it.
//!   middle 7  = information rings 3·5·7·11·13·17·19 (≈ 2^22.2: the value range; all ≤ 19, i.e. byte-shuffle sized)
//!   upper 3   = redundant routing primes 29·43·127 (2^14 ≡ ∓1, ±1, +1: Round2 is a pure route on them)
//!   lower 3   = redundant medium primes 113·241·251
//!   top 1     = redundant large prime 2^31 − 1 (2^31 ≡ 1: ×2^k is a rotation)
//!   bottom 1  = the 2^64 carrier (low bits for every Round2; also a full copy of the value)
//!   Redundant residue number system (RRNS): 7 redundant odd moduli, every one larger than every information modulus (the condition under
//!   which r redundant residues correct ⌊r/2⌋ wrong ones and detect r) → up to 3 wrong rings per value corrected in place.
//! Mode `base`: the 1-3-5-3-1 per-diamond placement of armilla_cores.rs (diamond 0 routing 127·113·43 + check 29; diamonds 1–8 two data +
//!   one check prime 1031…4093) with its sub-core (one check ring), in this same binary for a same-conditions comparison.
//! Faults: FAULT_RATE of blocks get K wrong rings (K distinct rings of one output value, the carrier counted as a ring) between the row
//! and column passes. Usage: armilla_layers base|v137 FAULT_RATE K [WORKERS]
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
fn inv(a: u64, m: u64) -> u64 { let (mut t, mut nt, mut r, mut nr) = (0i128, 1i128, m as i128, (a % m) as i128); while nr != 0 { let q = r / nr; (t, nt) = (nt, t - q * nt); (r, nr) = (nr, r - q * nr); } ((t % m as i128 + m as i128) % m as i128) as u64 }
fn plain_block(blk: &[i32; 64], ops: &[Op]) -> [i32; 64] {
    let mut d = *blk;
    let one = |t: &mut [i64; 8]| { let c = *t; for i in 0..8 { t[i] = c[brev(3, i)]; }
        for &op in ops { match op { Op::B(a, b, ang, f) => { let (x, y) = (t[a] * cos64(ang) - t[b] * sin64(ang), t[a] * sin64(ang) + t[b] * cos64(ang)); t[a] = round2(x, 14); t[b] = round2(y, 14); if f { t.swap(a, b); } }
                                     Op::H(a, b) => { let (x, y) = (t[a], t[b]); t[a] = x + y; t[b] = x - y; } } } };
    for i in 0..8 { let mut t = [0i64; 8]; for j in 0..8 { t[j] = d[i * 8 + j] as i64; } one(&mut t); for j in 0..8 { d[i * 8 + j] = t[j] as i32; } }
    for j in 0..8 { let mut t = [0i64; 8]; for i in 0..8 { t[i] = d[i * 8 + j] as i64; } one(&mut t); for i in 0..8 { d[i * 8 + j] = round2(t[i], 5) as i32; } }
    d
}

const KM: usize = 14; // max odd rings per diamond
type Res = [u64; KM];
/// CRT of the rings `idx` → signed value (None when the product overflows u128 headroom — not used with our sets)
fn crt(p: &[u64], r: &Res, idx: &[usize]) -> i128 {
    let prod: u128 = idx.iter().map(|&k| p[k] as u128).product();
    let mut x: u128 = 0;
    for &k in idx { let m = p[k] as u128; let mi = prod / m; let w = mi * inv((mi % m) as u64, p[k]) as u128 % prod; x = (x + (r[k] as u128 % m) * w % prod) % prod; }
    if x > prod / 2 { x as i128 - prod as i128 } else { x as i128 }
}
#[derive(Clone)]
struct Diamond { p: Vec<u64>, data: usize, inv14: Vec<u64>, range: i128, layer: Vec<&'static str>, w: Vec<u128>, prod: u128, helpers: usize }
impl Diamond {
    fn new(p: Vec<u64>, data: usize, layer: Vec<&'static str>) -> Diamond {
        let inv14 = p.iter().map(|&m| inv(16384 % m, m)).collect();
        let prod: u128 = p[..data].iter().map(|&m| m as u128).product();
        let w = p[..data].iter().map(|&m| { let mi = prod / m as u128; mi * inv((mi % m as u128) as u64, m) as u128 % prod }).collect();
        Diamond { p, data, inv14, range: (prod / 2) as i128, layer, w, prod, helpers: 1 }
    }
    fn with_helpers(mut self, h: usize) -> Diamond { self.helpers = h; self }
    /// CRT of the information rings with precomputed weights (the common path)
    fn decode(&self, r: &Res) -> i128 {
        let mut x: u128 = 0; for k in 0..self.data { x = (x + r[k] as u128 * self.w[k]) % self.prod; }
        if x > self.prod / 2 { x as i128 - self.prod as i128 } else { x as i128 }
    }
    fn n(&self) -> usize { self.p.len() }
    /// one 1D inverse DCT on the rings + rounding helper A (w) + helper B (w2, when helpers >= 2). Ring k takes its Round2 low bits from
    /// helper A when k < data (the main rings) and from helper B when k >= data (the reserved rings) — each helper feeds its own rings.
    fn idct(&self, r: &mut [Res; 8], w: &mut [u64; 8], w2: &mut [u64; 8], ops: &[Op]) {
        let n = self.n(); let two = self.helpers > 1;
        let (cr, cw, cw2) = (*r, *w, *w2); for i in 0..8 { r[i] = cr[brev(3, i)]; w[i] = cw[brev(3, i)]; w2[i] = cw2[brev(3, i)]; }
        for &op in ops { match op {
            Op::B(a, b, ang, flip) => {
                let (c, s) = (cos64(ang), sin64(ang));
                let xw = w[a].wrapping_mul(c as u64).wrapping_add(w[b].wrapping_mul((-s) as u64)).wrapping_add(8192);
                let yw = w[a].wrapping_mul(s as u64).wrapping_add(w[b].wrapping_mul(c as u64)).wrapping_add(8192);
                let (lx, ly) = (xw & 0x3FFF, yw & 0x3FFF);
                let xb = w2[a].wrapping_mul(c as u64).wrapping_add(w2[b].wrapping_mul((-s) as u64)).wrapping_add(8192);
                let yb = w2[a].wrapping_mul(s as u64).wrapping_add(w2[b].wrapping_mul(c as u64)).wrapping_add(8192);
                for k in 0..n {
                    let m = self.p[k];
                    let (lx, ly) = if two && k >= self.data { (xb & 0x3FFF, yb & 0x3FFF) } else { (lx, ly) };
                    let (cm, sm, nsm) = ((c.rem_euclid(m as i64)) as u64, (s.rem_euclid(m as i64)) as u64, ((-s).rem_euclid(m as i64)) as u64);
                    let x = ((r[a][k] * cm) % m + (r[b][k] * nsm) % m + 8192 % m + m - lx % m) % m;
                    let y = ((r[a][k] * sm) % m + (r[b][k] * cm) % m + 8192 % m + m - ly % m) % m;
                    r[a][k] = x * self.inv14[k] % m; r[b][k] = y * self.inv14[k] % m;
                }
                w[a] = ((xw as i64) >> 14) as u64; w[b] = ((yw as i64) >> 14) as u64;
                w2[a] = ((xb as i64) >> 14) as u64; w2[b] = ((yb as i64) >> 14) as u64;
                if flip { r.swap(a, b); w.swap(a, b); w2.swap(a, b); }
            }
            Op::H(a, b) => {
                for k in 0..n { let m = self.p[k]; let (x, y) = (r[a][k], r[b][k]); r[a][k] = (x + y) % m; r[b][k] = (x + m - y) % m; }
                let (x, y) = (w[a], w[b]); w[a] = x.wrapping_add(y); w[b] = x.wrapping_sub(y);
                let (x, y) = (w2[a], w2[b]); w2[a] = x.wrapping_add(y); w2[b] = x.wrapping_sub(y);
            }
        } }
    }
    /// fault: (output position, rings to corrupt — index n = the carrier)
    /// fault: (output position, slots to corrupt — ring k < n, n = rounding helper A (feeds every ring's rounding), n + 1 = helper B)
    fn block(&self, blk: &[i32; 64], ops: &[Op], fault: &Option<(usize, Vec<usize>)>) -> (Vec<Res>, [i64; 64], [i64; 64]) {
        let mut rr = vec![[0u64; KM]; 64]; let mut ww = [0u64; 64];
        for i in 0..64 { let x = blk[i] as i64; for k in 0..self.n() { rr[i][k] = x.rem_euclid(self.p[k] as i64) as u64; } ww[i] = x as u64; }
        let mut wb = ww;
        let pass = |rr: &mut Vec<Res>, ww: &mut [u64; 64], wb: &mut [u64; 64], rows: bool| for o in 0..8 {
            let at = |q: usize| if rows { o * 8 + q } else { q * 8 + o };
            let (mut r, mut w, mut w2) = ([[0u64; KM]; 8], [0u64; 8], [0u64; 8]);
            for q in 0..8 { r[q] = rr[at(q)]; w[q] = ww[at(q)]; w2[q] = wb[at(q)]; }
            self.idct(&mut r, &mut w, &mut w2, ops);
            for q in 0..8 { rr[at(q)] = r[q]; ww[at(q)] = w[q]; wb[at(q)] = w2[q]; } };
        pass(&mut rr, &mut ww, &mut wb, true);
        if let Some((pos, which)) = fault { for &k in which { if k < self.n() { rr[*pos][k] = (rr[*pos][k] + 1) % self.p[k]; } else if k == self.n() { ww[*pos] ^= 1 << 7; } else { wb[*pos] ^= 1 << 9; } } }
        pass(&mut rr, &mut ww, &mut wb, false);
        let (mut wi, mut wbi) = ([0i64; 64], [0i64; 64]); for i in 0..64 { wi[i] = ww[i] as i64; wbi[i] = wb[i] as i64; }
        (rr, wi, wbi)
    }
    fn disagree(&self, v: i128, r: &Res) -> usize {
        if let Ok(v) = i64::try_from(v) { (0..self.n()).filter(|&k| v.rem_euclid(self.p[k] as i64) as u64 != r[k]).count() }
        else { (0..self.n()).filter(|&k| v.rem_euclid(self.p[k] as i128) as u64 != r[k]).count() }
    }
}

/// base sub-core (armilla_cores.rs): data CRT vs carrier, one check ring decides
fn sub_core_base(d: &Diamond, rr: &[Res], w: &[i64; 64], blk: &[i32; 64], ops: &[Op]) -> ([i32; 64], u8) {
    let mut out = [0i32; 64]; let mut how = 0u8;
    for i in 0..64 {
        let vr = d.decode(&rr[i]) as i64; let vw = w[i];
        let ok_r = (d.data..d.n()).all(|k| vr.rem_euclid(d.p[k] as i64) as u64 == rr[i][k]);
        let ok_w = d.disagree(vw as i128, &rr[i]) <= 1;
        let v = if vr == vw && ok_r { vr } else if ok_w { how = how.max(1); vw } else if ok_r { how = how.max(1); vr } else { return (plain_block(blk, ops), 2) };
        out[i] = round2(v, 5) as i32;
    }
    (out, how)
}
/// v137 sub-core: RRNS with 7 redundant rings + the carrier. Accept a candidate value that lies in the information range and disagrees
/// with at most 3 rings (the wrong rings are corrected in place); candidates: the carrier, then CRT decodes of the information set with
/// 0, 1 or 2 of its rings swapped for redundant rings (covers 3 wrong rings with the carrier among them, or up to 2 wrong information rings).
fn sub_core_v(d: &Diamond, rr: &[Res], w: &[i64; 64], blk: &[i32; 64], ops: &[Op]) -> ([i32; 64], u8) {
    let mut out = [0i32; 64]; let mut how = 0u8; let t = (d.n() - d.data) / 2;
    for i in 0..64 {
        let vw = w[i] as i128;
        let ok = |v: i128| v.abs() <= d.range && d.disagree(v, &rr[i]) <= t;
        let v = if d.disagree(vw, &rr[i]) == 0 { vw } else if ok(vw) { how = how.max(1); vw } else {
            let mut found = None;
            'search: for drop in 0..=2usize {
                let info: Vec<usize> = (0..d.data).collect();
                let combos: Vec<Vec<usize>> = match drop { 0 => vec![vec![]], 1 => (0..d.data).map(|a| vec![a]).collect(),
                    _ => (0..d.data).flat_map(|a| (a + 1..d.data).map(move |b| vec![a, b])).collect() };
                for out_set in combos {
                    let mut idx: Vec<usize> = info.iter().copied().filter(|k| !out_set.contains(k)).collect();
                    idx.extend(d.data..d.data + out_set.len());
                    let v = crt(&d.p, &rr[i], &idx);
                    if ok(v) { found = Some(v); break 'search; }
                }
            }
            match found { Some(v) => { how = how.max(1); v } None => return (plain_block(blk, ops), 2) }
        };
        out[i] = round2(v as i64, 5) as i32;
    }
    (out, how)
}
/// mix sub-core (J 10-15 20:05). Per system: main rings = information, 2 reserved rings (the other size classes) = checks, 1 or 2 rounding
/// helpers. With 2 helpers each feeds its own rings: helper A the main rings' rounding, helper B the reserved rings' rounding — so a helper
/// fault bends only its own rings and the other side stays an independent witness. LEGAL RANGE: |coefficient| ≤ 512 and each 8-point
/// pass grows a value at most 8× → |v| ≤ 2^15 (checked with margin 2^16).
///   helpers agree and are legal → that value (disagreeing rings corrected in place)
///   2 helpers, they differ: A is confirmed when it is legal and a majority of the main rings agree (all of them when there is one); B is confirmed when it is legal and
///     agrees with both reserved rings; exactly one confirmed → it; otherwise the main rings' decode if every reserved ring confirms it;
///     else recompute (e.g. both helpers hit — detected, never guessed)
///   1 helper: agrees with every ring → it; main-ring decode confirmed by both reserved rings → it; helper agrees with all but ≤ 1 ring →
///     it; else recompute
fn sub_core_mix(d: &Diamond, rr: &[Res], w: &[i64; 64], w2: &[i64; 64], blk: &[i32; 64], ops: &[Op]) -> ([i32; 64], u8) {
    const LIM: i128 = 1 << 16;
    let mut out = [0i32; 64]; let mut how = 0u8;
    let legal = |v: i128| v.abs() <= LIM;
    for i in 0..64 {
        let r = &rr[i]; let (a, b) = (w[i] as i128, w2[i] as i128);
        let dis_main = |v: i128| (0..d.data).filter(|&k| v.rem_euclid(d.p[k] as i128) as u64 != r[k]).count();
        let dis_res = |v: i128| (d.data..d.n()).filter(|&k| v.rem_euclid(d.p[k] as i128) as u64 != r[k]).count();
        let vr = d.decode(r); let rings_ok = legal(vr) && dis_res(vr) == 0;
        let v = if d.helpers > 1 {
            if a == b && legal(a) { if dis_main(a) + dis_res(a) > 0 { how = how.max(1); } Some(a) }
            else {
                let ca = legal(a) && dis_main(a) <= (d.data - 1) / 2; let cb = legal(b) && dis_res(b) == 0;   // a minority of main rings may be wrong
                how = how.max(1);
                if ca && !cb { Some(a) } else if cb && !ca { Some(b) } else if rings_ok { Some(vr) } else { None }
            }
        } else {
            let da = dis_main(a) + dis_res(a);
            if da == 0 && legal(a) { Some(a) } else if rings_ok { how = how.max(1); Some(vr) } else if da <= 1 && legal(a) { how = how.max(1); Some(a) } else { None }
        };
        match v { Some(v) => out[i] = round2(v as i64, 5) as i32, None => return (plain_block(blk, ops), 2) }
    }
    (out, how)
}
fn checksum(o: &[i32; 64]) -> u64 { o.iter().fold(0u64, |a, &x| (a + x.rem_euclid(65521) as u64) % 65521) }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mode = a.get(1).map(|s| s.as_str()).unwrap_or("v137").to_string();
    let fault_rate: f64 = a.get(2).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let kf: usize = a.get(3).and_then(|v| v.parse().ok()).unwrap_or(1);
    let workers: usize = a.get(4).and_then(|v| v.parse().ok()).unwrap_or(8);
    let mut ops = vec![]; idct_ops(3, &mut ops);
    let isp = |n: u64| n > 1 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0);
    // systems per diamond: base / v137 = one ring set per diamond (every lane of the diamond carries it); mix = 15 systems, each its own set
    let sys: Vec<Vec<Diamond>> = if mode == "base" {
        let pool: Vec<u64> = (1031u64..4096).filter(|&n| isp(n)).collect();
        let mut v = vec![vec![Diamond::new(vec![127, 113, 43, 29], 3, vec!["data", "data", "data", "check"])]];
        for k in 1..9 { let j = (k - 1) * 3; v.push(vec![Diamond::new(vec![pool[j * 7], pool[j * 7 + 50], pool[j * 7 + 100]], 2, vec!["data", "data", "check"])]); }
        v
    } else if mode == "v137" {
        let p: Vec<u64> = vec![3, 5, 7, 11, 13, 17, 19, 29, 43, 127, 113, 241, 251, 2147483647];
        let layer: Vec<&'static str> = vec!["middle"; 7].into_iter().chain(vec!["upper"; 3]).chain(vec!["lower"; 3]).chain(vec!["top"]).collect();
        (0..9).map(|_| vec![Diamond::new(p.clone(), 7, layer.clone())]).collect()
    } else {
        // mix (J 10-15 20:05): every system = its main size class (the information rings: the value range) + one reserved ring of each other
        // size class (checks) + 1 or 2 rounding helpers. middle 7 = small main + medium + large; upper / lower 3 = medium main + small + large;
        // top / bottom = large main + medium + small.
        let h: usize = if mode == "mix2" { 2 } else { 1 };
        let large: Vec<u64> = (0..).map(|k| 2147483647u64 - k).filter(|&n| isp(n)).take(15).collect();
        let sp: Vec<u64> = vec![3, 5, 7, 11, 13, 17, 19, 23, 29, 31];
        let medr = [251u64, 241, 239, 233, 229, 227, 223];
        let medmain: [[u64; 3]; 6] = [[211, 199, 197], [193, 191, 181], [179, 173, 167], [163, 157, 151], [149, 139, 137], [131, 127, 113]];
        let smallr = [17u64, 13, 11, 19, 7, 5];
        let mut d = vec![];
        let top = Diamond::new(vec![large[13], 109, 3], 1, vec!["top"; 3]).with_helpers(h);
        d.push(top);
        for u in 0..3 { let mut p = medmain[u].to_vec(); p.extend([smallr[u], large[7 + u]]); d.push(Diamond::new(p, 3, vec!["upper"; 5]).with_helpers(h)); }
        for m in 0..7 {
            let mut main = vec![]; let mut prod = 1u64;
            for k in 0..sp.len() { if prod >= 1 << 22 { break; } let q = sp[(m + k) % sp.len()]; main.push(q); prod *= q; }
            main.sort(); let nm = main.len(); main.extend([medr[m], large[m]]);
            d.push(Diamond::new(main, nm, vec!["middle"; nm + 2]).with_helpers(h));
        }
        for u in 3..6 { let mut p = medmain[u].to_vec(); p.extend([smallr[u], large[7 + u]]); d.push(Diamond::new(p, 3, vec!["lower"; 5]).with_helpers(h)); }
        d.push(Diamond::new(vec![large[14], 107, 5], 1, vec!["bottom"; 3]).with_helpers(h));
        (0..9).map(|_| d.clone()).collect()
    };
    let (bw, bh) = (240usize, 135usize);
    let mut s: u64 = 0x9E3779B97F4A7C15;
    let mut rnd = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    let blocks: Vec<[i32; 64]> = (0..bw * bh).map(|_| { let mut b = [0i32; 64]; for (k, c) in b.iter_mut().enumerate() { *c = if k < 20 || rnd() % 4 == 0 { (rnd() % 1024) as i32 - 512 } else { 0 }; } b }).collect();
    let t0 = Instant::now(); let plain: Vec<[i32; 64]> = blocks.iter().map(|b| plain_block(b, &ops)).collect(); let tp = t0.elapsed().as_secs_f64();
    let diamond_of: Vec<usize> = (0..bw * bh).map(|bi| { let row = bi / bw; let g = row * 3 / bh; let d = (row * 3 % bh) * 3 / bh % 3; g * 3 + d }).collect();
    // block → system inside its diamond: smooth weighted round robin, weight = 1 / (rings + helpers) (the scalar cost)
    let mut sys_of = vec![0usize; bw * bh];
    { let mut credit: Vec<Vec<f64>> = sys.iter().map(|d| vec![0.0; d.len()]).collect();
      for bi in 0..bw * bh { let dd = diamond_of[bi]; let n = sys[dd].len(); if n == 1 { continue; }
        let wts: Vec<f64> = sys[dd].iter().map(|x| 1.0 / (x.n() + x.helpers) as f64).collect(); let tot: f64 = wts.iter().sum();
        for k in 0..n { credit[dd][k] += wts[k]; }
        let best = (0..n).max_by(|&x, &y| credit[dd][x].partial_cmp(&credit[dd][y]).unwrap()).unwrap(); credit[dd][best] -= tot; sys_of[bi] = best; } }
    let mut faults: Vec<Option<(usize, Vec<usize>)>> = vec![None; bw * bh];
    if fault_rate > 0.0 { for (bi, f) in faults.iter_mut().enumerate() { if (rnd() % 1_000_000) as f64 / 1e6 < fault_rate {
        let x = &sys[diamond_of[bi]][sys_of[bi]]; let slots = x.n() + x.helpers; let mut which: Vec<usize> = vec![];
        while which.len() < kf.min(slots) { let k = (rnd() % slots as u64) as usize; if !which.contains(&k) { which.push(k); } }
        *f = Some(((rnd() % 64) as usize, which)); } } }
    let md = mode.clone();
    let t1 = Instant::now();
    let (final_tx, final_rx) = mpsc::channel::<(usize, [i32; 64], u64, u8)>();
    let report = std::thread::scope(|sc| {
        let mut core_tx = vec![]; let mut core_hs = vec![];
        for _c in 0..3 {
            let (tx, rx) = mpsc::channel::<(usize, Vec<Res>, [i64; 64], [i64; 64])>();
            core_tx.push(tx);
            let (sys, blocks, diamond_of, sys_of, ops, ftx, md) = (&sys, &blocks, &diamond_of, &sys_of, &ops, final_tx.clone(), md.clone());
            core_hs.push(sc.spawn(move || { for (bi, rr, w, w2) in rx { let d = &sys[diamond_of[bi]][sys_of[bi]];
                let (o, how) = match md.as_str() { "base" => sub_core_base(d, &rr, &w, &blocks[bi], ops), "v137" => sub_core_v(d, &rr, &w, &blocks[bi], ops),
                                                   _ => sub_core_mix(d, &rr, &w, &w2, &blocks[bi], ops) };
                ftx.send((bi, o, checksum(&o), how)).unwrap(); } }));
        }
        drop(final_tx);
        let chunk = (bw * bh + workers - 1) / workers;
        let whs: Vec<_> = (0..workers).map(|t| { let (sys, blocks, diamond_of, sys_of, faults, ops, txs) = (&sys, &blocks, &diamond_of, &sys_of, &faults, &ops, core_tx.clone());
            sc.spawn(move || for bi in t * chunk..((t + 1) * chunk).min(bw * bh) { let d = diamond_of[bi]; let (rr, w, w2) = sys[d][sys_of[bi]].block(&blocks[bi], ops, &faults[bi]); txs[d / 3].send((bi, rr, w, w2)).unwrap(); }) }).collect();
        drop(core_tx);
        let mut frame = vec![[0i32; 64]; bw * bh]; let mut how_of = vec![0u8; bw * bh]; let mut bad_sum = 0usize;
        for (bi, o, cs, how) in final_rx { if checksum(&o) != cs { bad_sum += 1; } frame[bi] = o; how_of[bi] = how; }
        for h in whs { h.join().unwrap(); } for h in core_hs { h.join().unwrap(); }
        (frame, how_of, bad_sum)
    });
    let ta = t1.elapsed().as_secs_f64();
    let (frame, how_of, bad_sum) = report;
    let wrong_v: Vec<bool> = frame.iter().zip(&plain).map(|(a, b)| a != b).collect();
    let injected = faults.iter().filter(|f| f.is_some()).count();
    let helper_hit = (0..bw * bh).filter(|&bi| { let x = &sys[diamond_of[bi]][sys_of[bi]]; faults[bi].as_ref().map_or(false, |(_, w)| w.iter().any(|&k| k >= x.n())) }).count();
    // per layer: blocks, faulted, corrected, recomputed, wrong
    let mut per: std::collections::BTreeMap<&str, [usize; 5]> = Default::default();
    for bi in 0..bw * bh { let x = &sys[diamond_of[bi]][sys_of[bi]]; let e = per.entry(x.layer[0]).or_default();
        e[0] += 1; if faults[bi].is_some() { e[1] += 1; } if how_of[bi] == 1 { e[2] += 1; } if how_of[bi] == 2 { e[3] += 1; } if wrong_v[bi] { e[4] += 1; } }
    let corrected = how_of.iter().filter(|&&h| h == 1).count(); let recomputed = how_of.iter().filter(|&&h| h == 2).count(); let wrong = wrong_v.iter().filter(|&&w| w).count();
    let systems: Vec<String> = sys[0].iter().map(|x| format!("{{\"layer\":\"{}\",\"main\":{:?},\"reserved\":{:?},\"helpers\":{}}}", x.layer[0], &x.p[..x.data], &x.p[x.data..], x.helpers)).collect();
    let per_s: Vec<String> = per.iter().map(|(k, v)| format!("\"{k}\":{{\"blocks\":{},\"faulted\":{},\"corrected\":{},\"recomputed\":{},\"wrong\":{}}}", v[0], v[1], v[2], v[3], v[4])).collect();
    println!("{{\"probe\":\"armilla-layers-vp9-idct8x8\",\"mode\":\"{mode}\",\"workers\":{workers},\"blocks\":{},\"fault_rate\":{fault_rate},\"rings_per_fault\":{kf},\"blocks_wrong\":{wrong},\"faults_injected\":{injected},\"faults_hitting_a_helper\":{helper_hit},\"corrected_in_place\":{corrected},\"recomputed\":{recomputed},\"final_checksum_failures\":{bad_sum},\"armilla_s\":{ta:.4},\"plain_s\":{tp:.4},\"ratio\":{:.1},\"frames_per_s\":{:.2},\"per_layer\":{{{}}},\"diamond0_systems\":[{}]}}",
        bw * bh, ta / tp, 1.0 / ta, per_s.join(","), systems.join(","));
}
