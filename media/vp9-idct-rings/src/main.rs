//! MEDIA-PLAYER-0 Armilla lead (J 10-15, readback confirmed: "the decode as exact integer arithmetic on Armilla's rings, in the
//! renderer's core"). Cheapest probe: the VP9 8-point inverse DCT (VP9 Bitstream & Decoding Process Specification v0.6 §8.7.1:
//! butterflies B / H, cos64 table, Round2(x, 14), the §8.7.1.2 bit-reversal permutation, the §8.7.1.3 recursive process) computed
//! (a) on plain i32 and (b) with every value held as residues on 5 odd-prime rings + one 2^16 ring. Round2 on the rings: the 2^16 ring
//! gives the low 14 bits, the odd rings divide exactly by 2^14, and the 2^16 residue of the quotient is rebuilt from the odd rings by
//! mixed-radix base extension. Question: bit-exact? how fast vs plain?
use std::time::Instant;

const COS64: [i64; 33] = [16384, 16364, 16305, 16207, 16069, 15893, 15679, 15426, 15137, 14811, 14449, 14053, 13623, 13160, 12665, 12140,
    11585, 11003, 10394, 9760, 9102, 8423, 7723, 7005, 6270, 5520, 4756, 3981, 3196, 2404, 1606, 804, 0];
fn cos64(a: i64) -> i64 { let a2 = a & 127; match a2 { 0..=32 => COS64[a2 as usize], 33..=64 => -COS64[(64 - a2) as usize], 65..=96 => -COS64[(a2 - 64) as usize], _ => COS64[(128 - a2) as usize] } }
fn sin64(a: i64) -> i64 { cos64(a - 32) }
fn brev(bits: u32, x: usize) -> usize { let mut t = 0; for i in 0..bits { t |= ((x >> i) & 1) << (bits - 1 - i); } t }

#[derive(Clone, Copy, Debug)]
enum Op { B(usize, usize, i64, bool), H(usize, usize) }

/// §8.7.1.3 as an op list (the spec procedure, recursion included), n = log2 length
fn idct_ops(n: usize, ops: &mut Vec<Op>) {
    let (n0, n1, n2) = (1usize << n, 1usize << (n - 1), 1usize << (n - 2));
    let n3 = if n >= 3 { 1usize << (n - 3) } else { 0 };
    if n == 2 { ops.push(Op::B(0, 1, 16, true)); } else { idct_ops(n - 1, ops); }
    for i in 0..n2 { ops.push(Op::B(n1 + i, n0 - 1 - i, 32 - brev(5, n1 + i) as i64, false)); }
    if n >= 3 { for i in 0..n3 { for j in 0..2 { let (a, b) = (n1 + 4 * i + 2 * j, n1 + 1 + 4 * i + 2 * j); ops.push(if j == 1 { Op::H(b, a) } else { Op::H(a, b) }); } } }
    assert!(n <= 3, "probe covers n = 2, 3 (steps 4 and 5 are for 16 / 32 points)");
    if n >= 3 { for i in 0..n3 { ops.push(Op::B(n0 - n3 - 1 - i, n1 + n3 + i, 16, true)); } }
    for i in 0..n1 { ops.push(Op::H(i, n0 - 1 - i)); }
}
fn round2(x: i64, n: u32) -> i64 { (x + (1 << (n - 1))) >> n }

fn plain(t: &mut [i32; 8], ops: &[Op]) {
    let c = *t; for i in 0..8 { t[i] = c[brev(3, i)]; }
    for &op in ops { match op {
        Op::B(a, b, ang, flip) => { let (ta, tb) = (t[a] as i64, t[b] as i64);
            let x = ta * cos64(ang) - tb * sin64(ang); let y = ta * sin64(ang) + tb * cos64(ang);
            t[a] = round2(x, 14) as i32; t[b] = round2(y, 14) as i32; if flip { t.swap(a, b); } }
        Op::H(a, b) => { let (x, y) = (t[a], t[b]); t[a] = x + y; t[b] = x - y; }
    } }
}

// ── rings: 5 odd primes (product ≈ 2^39.5 > the ±2^31 range of x) + 2^16 ──
const M: [u64; 5] = [251, 241, 239, 233, 229];
const P2: u64 = 1 << 16;
struct Rings { inv14: [u64; 5], mrc_inv: [[u64; 5]; 5], prod: u128 }
fn inv(a: u64, m: u64) -> u64 { let (mut t, mut nt, mut r, mut nr) = (0i64, 1i64, m as i64, (a % m) as i64); while nr != 0 { let q = r / nr; (t, nt) = (nt, t - q * nt); (r, nr) = (nr, r - q * nr); } ((t % m as i64 + m as i64) % m as i64) as u64 }
impl Rings {
    fn new() -> Rings {
        let mut mrc_inv = [[0u64; 5]; 5];
        for i in 0..5 { for j in 0..i { mrc_inv[i][j] = inv(M[j] % M[i], M[i]); } }
        Rings { inv14: M.map(|m| inv(16384 % m, m)), mrc_inv, prod: M.iter().map(|&m| m as u128).product() }
    }
    fn enc(&self, x: i64) -> ([u64; 5], u64) { (M.map(|m| x.rem_euclid(m as i64) as u64), (x as u64) & (P2 - 1)) }
    /// mixed-radix reconstruction of the signed value from the odd rings (the base-extension bridge)
    fn dec(&self, r: &[u64; 5]) -> i64 {
        let mut a = [0u64; 5];
        for i in 0..5 { let mut v = r[i]; for j in 0..i { v = (v + M[i] - a[j] % M[i]) % M[i] * self.mrc_inv[i][j] % M[i]; } a[i] = v; }
        let mut x: u128 = 0; let mut w: u128 = 1;
        for i in 0..5 { x += a[i] as u128 * w; w *= M[i] as u128; }
        if x > self.prod / 2 { x as i128 as i64 - self.prod as i64 } else { x as i64 }
    }
}

/// one Round2(x, 14) on residues: (odd residues, 2^16 residue) of x → those of Round2(x, 14)
fn ring_round2(rg: &Rings, r: &mut [u64; 5], p: &mut u64) {
    for i in 0..5 { r[i] = (r[i] + 8192) % M[i]; }
    *p = (*p + 8192) & (P2 - 1);
    let low = *p & 0x3FFF;
    for i in 0..5 { r[i] = (r[i] + M[i] - low % M[i]) % M[i] * rg.inv14[i] % M[i]; }
    *p = (rg.dec(r) as u64) & (P2 - 1);
}

fn rings(rg: &Rings, t: &[i32; 8], ops: &[Op]) -> [i32; 8] {
    let mut r = [[0u64; 5]; 8]; let mut p = [0u64; 8];
    for i in 0..8 { let (a, b) = rg.enc(t[brev(3, i)] as i64); r[i] = a; p[i] = b; }
    let cm = |v: i64, m: u64| v.rem_euclid(m as i64) as u64;
    for &op in ops { match op {
        Op::B(a, b, ang, flip) => {
            let (c, s) = (cos64(ang), sin64(ang));
            let (mut xr, mut yr) = ([0u64; 5], [0u64; 5]);
            for i in 0..5 { let m = M[i]; xr[i] = (r[a][i] * cm(c, m) + r[b][i] * cm(-s, m)) % m; yr[i] = (r[a][i] * cm(s, m) + r[b][i] * cm(c, m)) % m; }
            let mut xp = (p[a].wrapping_mul(c as u64).wrapping_add(p[b].wrapping_mul((-s) as u64))) & (P2 - 1);
            let mut yp = (p[a].wrapping_mul(s as u64).wrapping_add(p[b].wrapping_mul(c as u64))) & (P2 - 1);
            ring_round2(rg, &mut xr, &mut xp); ring_round2(rg, &mut yr, &mut yp);
            r[a] = xr; r[b] = yr; p[a] = xp; p[b] = yp;
            if flip { r.swap(a, b); p.swap(a, b); }
        }
        Op::H(a, b) => {
            let (ra, rb) = (r[a], r[b]);
            for i in 0..5 { r[a][i] = (ra[i] + rb[i]) % M[i]; r[b][i] = (ra[i] + M[i] - rb[i]) % M[i]; }
            let (pa, pb) = (p[a], p[b]); p[a] = (pa + pb) & (P2 - 1); p[b] = pa.wrapping_sub(pb) & (P2 - 1);
        }
    } }
    let mut out = [0i32; 8]; for i in 0..8 { out[i] = rg.dec(&r[i]) as i32; } out
}


// ── mode 2 (J 10-15: "choose a sorting / organization array in the rings, routed so the rounding is removed"): rings that divide
// 2^28 − 1 = 3·5·29·43·113·127. There 2^14 ≡ +1 (3, 43, 127) or −1 (5, 29, 113), so "× inverse of 2^14" is the identity or a negation:
// Round2 on these rings = add 8192, subtract the low 14 bits, route straight or negated. The low bits come from a 2^64 ring (wrapping
// u64): every T value is < 2^16 (spec bound), so the 2^64 ring always holds the exact value and never needs a base extension.
const R: [u64; 6] = [3, 43, 127, 5, 29, 113];
const NEG: [bool; 6] = [false, false, false, true, true, true];
fn routed(t: &[i32; 8], ops: &[Op]) -> ([i32; 8], bool) { routed_f(t, ops, None) }
/// fault: Some((op index, slot, bit)) flips that bit of the 2^64-ring value after that op (a corrupted register / memory word)
fn routed_f(t: &[i32; 8], ops: &[Op], fault: Option<(usize, usize, u32)>) -> ([i32; 8], bool) {
    let mut r = [[0u64; 6]; 8]; let mut w = [0u64; 8];
    for i in 0..8 { let x = t[brev(3, i)] as i64; for k in 0..6 { r[i][k] = x.rem_euclid(R[k] as i64) as u64; } w[i] = x as u64; }
    let cm = |v: i64, m: u64| v.rem_euclid(m as i64) as u64;
    let rnd = |rr: &mut [u64; 6], ww: &mut u64| {
        *ww = ww.wrapping_add(8192); let low = *ww & 0x3FFF;
        for k in 0..6 { let m = R[k]; let v = (rr[k] + 8192 % m + m - low % m) % m; rr[k] = if NEG[k] { (m - v) % m } else { v }; }
        *ww = ((*ww as i64) >> 14) as u64; // the 2^64 ring: exact (|x| < 2^31), so the shift IS Round2
    };
    for (oi, &op) in ops.iter().enumerate() { if let Some((fo, fs, fb)) = fault { if fo == oi { w[fs] ^= 1u64 << fb; } } match op {
        Op::B(a, b, ang, flip) => {
            let (c, s) = (cos64(ang), sin64(ang));
            let (mut xr, mut yr) = ([0u64; 6], [0u64; 6]);
            for k in 0..6 { let m = R[k]; xr[k] = (r[a][k] * cm(c, m) + r[b][k] * cm(-s, m)) % m; yr[k] = (r[a][k] * cm(s, m) + r[b][k] * cm(c, m)) % m; }
            let mut xw = w[a].wrapping_mul(c as u64).wrapping_add(w[b].wrapping_mul((-s) as u64));
            let mut yw = w[a].wrapping_mul(s as u64).wrapping_add(w[b].wrapping_mul(c as u64));
            rnd(&mut xr, &mut xw); rnd(&mut yr, &mut yw);
            r[a] = xr; r[b] = yr; w[a] = xw; w[b] = yw;
            if flip { r.swap(a, b); w.swap(a, b); }
        }
        Op::H(a, b) => {
            let (ra, rb) = (r[a], r[b]);
            for k in 0..6 { r[a][k] = (ra[k] + rb[k]) % R[k]; r[b][k] = (ra[k] + R[k] - rb[k]) % R[k]; }
            let (wa, wb) = (w[a], w[b]); w[a] = wa.wrapping_add(wb); w[b] = wa.wrapping_sub(wb);
        }
    } }
    // out = the 2^64 ring's value; the six routed rings must agree with it (a free fault check: any corrupted value disagrees)
    let mut out = [0i32; 8]; let mut agree = true;
    for i in 0..8 { let x = w[i] as i64; out[i] = x as i32; for k in 0..6 { if x.rem_euclid(R[k] as i64) as u64 != r[i][k] { agree = false; } } }
    (out, agree)
}

// ── mode 3 (J 10-15: "the 3x3x3 with math core … the output all routed"): the C5 three-group format of ARMILLA-CALC-0 on T50 primes.
// A = routing rings 3·43·127·5·29·113 (2^14 ≡ ±1: Round2 = route), B = T50 primes 7·11·13·17·19·23 (Round2 = × inverse of 2^14, one
// constant multiply), C = guard 2^64 ring (low bits; the value itself). Centre: each group decodes its own output (A range 2^28−1,
// B range 7,436,429, both > the 17-bit output range) and a 2-of-3 vote CORRECTS any single faulty group.
const GB: [u64; 6] = [7, 11, 13, 17, 19, 23];
fn crt(r: &[u64], m: &[u64]) -> i64 {
    let big: u128 = m.iter().map(|&x| x as u128).product();
    let mut x: u128 = 0;
    for i in 0..m.len() { let mi = big / m[i] as u128; let inv_ = inv((mi % m[i] as u128) as u64, m[i]) as u128; x = (x + r[i] as u128 * mi % big * inv_) % big; }
    if x > big / 2 { x as i128 as i64 - big as i64 } else { x as i64 }
}
/// fault: Some((group 0=A 1=B 2=C, op, slot, which residue / bit))
fn three(t: &[i32; 8], ops: &[Op], fault: Option<(usize, usize, usize, u32)>) -> ([i32; 8], u32) {
    let invb: [u64; 6] = GB.map(|m| inv(16384 % m, m));
    let mut a = [[0u64; 6]; 8]; let mut b = [[0u64; 6]; 8]; let mut w = [0u64; 8];
    for i in 0..8 { let x = t[brev(3, i)] as i64; for k in 0..6 { a[i][k] = x.rem_euclid(R[k] as i64) as u64; b[i][k] = x.rem_euclid(GB[k] as i64) as u64; } w[i] = x as u64; }
    let cm = |v: i64, m: u64| v.rem_euclid(m as i64) as u64;
    for (oi, &op) in ops.iter().enumerate() {
        if let Some((g, fo, fs, fb)) = fault { if fo == oi { match g { 0 => { let k = fb as usize % 6; a[fs][k] = (a[fs][k] + 1) % R[k]; } 1 => { let k = fb as usize % 6; b[fs][k] = (b[fs][k] + 1) % GB[k]; } _ => w[fs] ^= 1u64 << (fb % 32) } } }
        match op {
        Op::B(p, q, ang, flip) => {
            let (c, s) = (cos64(ang), sin64(ang));
            let mut wx = w[p].wrapping_mul(c as u64).wrapping_add(w[q].wrapping_mul((-s) as u64)).wrapping_add(8192);
            let mut wy = w[p].wrapping_mul(s as u64).wrapping_add(w[q].wrapping_mul(c as u64)).wrapping_add(8192);
            let (lx, ly) = (wx & 0x3FFF, wy & 0x3FFF);
            wx = ((wx as i64) >> 14) as u64; wy = ((wy as i64) >> 14) as u64;
            let (mut ax, mut ay, mut bx, mut by) = ([0u64; 6], [0u64; 6], [0u64; 6], [0u64; 6]);
            for k in 0..6 {
                let m = R[k];
                let vx = (a[p][k] * cm(c, m) + a[q][k] * cm(-s, m) + 8192 % m + m - lx % m) % m;
                let vy = (a[p][k] * cm(s, m) + a[q][k] * cm(c, m) + 8192 % m + m - ly % m) % m;
                ax[k] = if NEG[k] { (m - vx) % m } else { vx }; ay[k] = if NEG[k] { (m - vy) % m } else { vy };
                let m = GB[k];
                bx[k] = (b[p][k] * cm(c, m) + b[q][k] * cm(-s, m) + 8192 % m + m - lx % m) % m * invb[k] % m;
                by[k] = (b[p][k] * cm(s, m) + b[q][k] * cm(c, m) + 8192 % m + m - ly % m) % m * invb[k] % m;
            }
            a[p] = ax; a[q] = ay; b[p] = bx; b[q] = by; w[p] = wx; w[q] = wy;
            if flip { a.swap(p, q); b.swap(p, q); w.swap(p, q); }
        }
        Op::H(p, q) => {
            for k in 0..6 { let (x, y) = (a[p][k], a[q][k]); a[p][k] = (x + y) % R[k]; a[q][k] = (x + R[k] - y) % R[k];
                            let (x, y) = (b[p][k], b[q][k]); b[p][k] = (x + y) % GB[k]; b[q][k] = (x + GB[k] - y) % GB[k]; }
            let (x, y) = (w[p], w[q]); w[p] = x.wrapping_add(y); w[q] = x.wrapping_sub(y);
        }
    } }
    // centre: three decodes, 2-of-3 vote; flags: bit0 = a group disagreed (detected), bit1 = no majority
    let mut out = [0i32; 8]; let mut flags = 0u32;
    for i in 0..8 {
        let (va, vb, vc) = (crt(&a[i], &R), crt(&b[i], &GB), w[i] as i64);
        if va != vb || vb != vc { flags |= 1; }
        out[i] = if va == vb || va == vc { va as i32 } else if vb == vc { vb as i32 } else { flags |= 2; vc as i32 };
    }
    (out, flags)
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(1_000_000);
    let mut ops = vec![]; idct_ops(3, &mut ops);
    let rg = Rings::new();
    // coefficients: 8-bit depth → the spec bounds T to 16-bit signed; inputs drawn so every stage stays in range (|c| ≤ 2^12)
    let mut s: u64 = 0x9E3779B97F4A7C15;
    let mut rnd = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; ((s >> 11) % 8192) as i32 - 4096 };
    let input: Vec<[i32; 8]> = (0..n).map(|_| { let mut v = [0i32; 8]; for x in v.iter_mut() { *x = rnd(); } v }).collect();
    let t0 = Instant::now();
    let plain_out: Vec<[i32; 8]> = input.iter().map(|v| { let mut t = *v; plain(&mut t, &ops); t }).collect();
    let tp = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let ring_out: Vec<[i32; 8]> = input.iter().map(|v| rings(&rg, v, &ops)).collect();
    let tr = t1.elapsed().as_secs_f64();
    let bad = plain_out.iter().zip(&ring_out).filter(|(a, b)| a != b).count();
    let t2 = Instant::now();
    let routed_out: Vec<([i32; 8], bool)> = input.iter().map(|v| routed(v, &ops)).collect();
    let tq = t2.elapsed().as_secs_f64();
    let bad2 = plain_out.iter().zip(&routed_out).filter(|(a, b)| **a != b.0).count();
    let disagree = routed_out.iter().filter(|x| !x.1).count();
    let t3 = Instant::now();
    let three_out: Vec<([i32; 8], u32)> = input.iter().map(|v| three(v, &ops, None)).collect();
    let tt = t3.elapsed().as_secs_f64();
    let bad3 = plain_out.iter().zip(&three_out).filter(|(a, b)| **a != b.0).count();
    let alarms = three_out.iter().filter(|x| x.1 != 0).count();
    let (mut f_changed_raw, mut f_out_ok, mut f_detected) = (0usize, 0usize, 0usize);
    let mut s3: u64 = 0x5DEECE66D; let mut nx3 = || { s3 ^= s3 << 13; s3 ^= s3 >> 7; s3 ^= s3 << 17; s3 };
    let (mut wrong_by_group, mut inj_by_group) = ([0usize; 3], [0usize; 3]);
    for (k, v) in input.iter().enumerate().take(200_000) {
        let f = ((nx3() % 3) as usize, nx3() as usize % ops.len(), nx3() as usize % 8, (nx3() % 64) as u32);
        let (o, fl) = three(v, &ops, Some(f));
        if fl & 1 != 0 { f_detected += 1; }
        if o == plain_out[k] { f_out_ok += 1; } else { f_changed_raw += 1; wrong_by_group[f.0] += 1; }
        inj_by_group[f.0] += 1;
    }
    println!("{{\"probe\":\"vp9-idct8-three-groups\",\"vectors\":{n},\"mismatches\":{bad3},\"clean_false_alarms\":{alarms},\"three_s\":{tt:.4},\"ratio_vs_plain\":{:.2},\"ns_per_vec\":{:.1},\"faults_injected\":200000,\"faults_detected\":{f_detected},\"outputs_correct_after_vote\":{f_out_ok},\"outputs_wrong\":{f_changed_raw},\"wrong_by_group_ABC\":{:?},\"injected_by_group_ABC\":{:?}}}", tt / tp, tt * 1e9 / n as f64, wrong_by_group, inj_by_group);
    // fault check: one bit flipped in one 2^64-ring value, at a random op / slot / bit (0..31), per vector
    let mut s2: u64 = 0x2545F4914F6CDD1D; let mut nx = || { s2 ^= s2 << 13; s2 ^= s2 >> 7; s2 ^= s2 << 17; s2 };
    let (mut changed, mut caught) = (0usize, 0usize);
    for (k, v) in input.iter().enumerate().take(200_000) {
        let f = (nx() as usize % ops.len(), nx() as usize % 8, (nx() % 32) as u32);
        let (o, agree) = routed_f(v, &ops, Some(f));
        if o != plain_out[k] { changed += 1; if !agree { caught += 1; } }
    }
    println!("{{\"fault_injection\":200000,\"outputs_changed\":{changed},\"caught_by_routed_rings\":{caught}}}");
    // fault check: flip one bit of one 2^64-ring value mid-way is modelled by corrupting the output and re-checking the rings
    println!("{{\"probe\":\"vp9-idct8-routed\",\"vectors\":{n},\"mismatches\":{bad2},\"ring_disagreements\":{disagree},\"routed_s\":{tq:.4},\"ratio_vs_plain\":{:.2},\"routed_ns_per_vec\":{:.1}}}", tq / tp, tq * 1e9 / n as f64);
    let ck: i64 = plain_out.iter().flat_map(|v| v.iter()).map(|&x| x as i64).sum();
    println!("{{\"probe\":\"vp9-idct8-rings\",\"vectors\":{n},\"ops\":{},\"mismatches\":{bad},\"plain_s\":{tp:.4},\"rings_s\":{tr:.4},\"ratio\":{:.1},\"plain_ns_per_vec\":{:.1},\"rings_ns_per_vec\":{:.1},\"checksum\":{ck}}}",
        ops.len(), tr / tp, tp * 1e9 / n as f64, tr * 1e9 / n as f64);
}
