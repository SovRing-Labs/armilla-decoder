//! J 10-15 (readbacks confirmed): choose which coprimes sit on which rings, with a ratio per diamond so each diamond computes its share by
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
    fn block(&self, blk: &[i32; 64], ops: &[Op], fault: Option<(usize, usize)>) -> (Vec<[u64; KM]>, [i64; 64]) {
        let mut rr = vec![[0u64; KM]; 64]; let mut ww = [0u64; 64];
        for i in 0..64 { let x = blk[i] as i64; for k in 0..self.n() { rr[i][k] = x.rem_euclid(self.p[k] as i64) as u64; } ww[i] = x as u64; }
        for i in 0..8 { let mut r = [[0u64; KM]; 8]; let mut w = [0u64; 8]; for j in 0..8 { r[j] = rr[i * 8 + j]; w[j] = ww[i * 8 + j]; }
            self.idct(&mut r, &mut w, ops); for j in 0..8 { rr[i * 8 + j] = r[j]; ww[i * 8 + j] = w[j]; } }
        if let Some((pos, which)) = fault { if which < self.n() { rr[pos][which] = (rr[pos][which] + 1) % self.p[which]; } else { ww[pos] ^= 1 << 7; } }
        for j in 0..8 { let mut r = [[0u64; KM]; 8]; let mut w = [0u64; 8]; for i in 0..8 { r[i] = rr[i * 8 + j]; w[i] = ww[i * 8 + j]; }
            self.idct(&mut r, &mut w, ops); for i in 0..8 { rr[i * 8 + j] = r[i]; ww[i * 8 + j] = w[i]; } }
        let mut wi = [0i64; 64]; for i in 0..64 { wi[i] = ww[i] as i64; }
        (rr, wi)
    }
}

/// a sub-core (diamond or group math core): check, correct in place, final Round2(x,5); returns (block, how: 0 clean 1 corrected 2 recomputed)
fn sub_core(d: &Diamond, rr: &[[u64; KM]], w: &[i64; 64], blk: &[i32; 64], ops: &[Op]) -> ([i32; 64], u8) {
    let mut out = [0i32; 64]; let mut how = 0u8;
    for i in 0..64 {
        let vr = d.decode(&rr[i]); let vw = w[i];
        let ok_r = (d.data..d.n()).all(|k| vr.rem_euclid(d.p[k] as i64) as u64 == rr[i][k]);   // residue copy consistent with its check ring
        // carrier confirmed when it agrees with all rings but at most one (that ring is the faulty one)
        let ok_w = (0..d.n()).filter(|&k| vw.rem_euclid(d.p[k] as i64) as u64 != rr[i][k]).count() <= 1;
        let v = if vr == vw && ok_r { vr } else if ok_w { how = how.max(1); vw } else if ok_r { how = how.max(1); vr } else { return (plain_block(blk, ops), 2) };
        out[i] = round2(v, 5) as i32;
    }
    (out, how)
}
fn checksum(o: &[i32; 64]) -> u64 { o.iter().fold(0u64, |a, &x| (a + x.rem_euclid(65521) as u64) % 65521) }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let layout = a.get(1).map(|s| s.as_str()).unwrap_or("A"); // A: 9 diamond cores + final; B: 3 group cores + final
    let fault_rate: f64 = a.get(2).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let workers: usize = a.get(3).and_then(|v| v.parse().ok()).unwrap_or(8);
    let mut ops = vec![]; idct_ops(3, &mut ops);
    // ring placement: diamond 0 = routing primes; diamonds 1..8 = 2 data + 1 check from the primes 1031..4093, disjoint
    let isp = |n: u64| n > 1 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0);
    let pool: Vec<u64> = (1031u64..4096).filter(|&n| isp(n)).collect();
    let mut diamonds = vec![Diamond::new(vec![127, 113, 43, 29], 3)];
    for k in 1..9 { let j = (k - 1) * 3; diamonds.push(Diamond::new(vec![pool[j * 7], pool[j * 7 + 50], pool[j * 7 + 100]], 2)); }
    let (bw, bh) = (240usize, 135usize);
    let mut s: u64 = 0x9E3779B97F4A7C15;
    let mut rnd = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    let blocks: Vec<[i32; 64]> = (0..bw * bh).map(|_| { let mut b = [0i32; 64]; for (k, c) in b.iter_mut().enumerate() { *c = if k < 20 || rnd() % 4 == 0 { (rnd() % 1024) as i32 - 512 } else { 0 }; } b }).collect();
    let t0 = Instant::now(); let plain: Vec<[i32; 64]> = blocks.iter().map(|b| plain_block(b, &ops)).collect(); let tp = t0.elapsed().as_secs_f64();
    // shares: band of block rows → group; sub-band → diamond (systems = lanes inside the diamond)
    let diamond_of: Vec<usize> = (0..bw * bh).map(|bi| { let row = bi / bw; let g = row * 3 / bh; let d = (row * 3 % bh) * 3 / bh % 3; g * 3 + d }).collect();
    let mut faults: Vec<Option<(usize, usize)>> = vec![None; bw * bh];
    if fault_rate > 0.0 { for (bi, f) in faults.iter_mut().enumerate() { if (rnd() % 1_000_000) as f64 / 1e6 < fault_rate { let n = diamonds[diamond_of[bi]].n(); *f = Some(((rnd() % 64) as usize, (rnd() % (n as u64 + 1)) as usize)); } } }
    let ncore = if layout == "A" { 9 } else { 3 };
    let core_of = |d: usize| if layout == "A" { d } else { d / 3 };
    let t1 = Instant::now();
    let (final_tx, final_rx) = mpsc::channel::<(usize, [i32; 64], u64, u8)>();
    let report = std::thread::scope(|sc| {
        let mut core_tx = vec![]; let mut core_hs = vec![];
        for _c in 0..ncore {
            let (tx, rx) = mpsc::channel::<(usize, Vec<[u64; KM]>, [i64; 64])>();
            core_tx.push(tx);
            let (diamonds, blocks, diamond_of, ops, ftx) = (&diamonds, &blocks, &diamond_of, &ops, final_tx.clone());
            core_hs.push(sc.spawn(move || { for (bi, rr, w) in rx { let d = &diamonds[diamond_of[bi]]; let (o, how) = sub_core(d, &rr, &w, &blocks[bi], ops); ftx.send((bi, o, checksum(&o), how)).unwrap(); } }));
        }
        drop(final_tx);
        // the systems (lanes) compute their shares on the worker threads and send to their sub-core
        let chunk = (bw * bh + workers - 1) / workers;
        let whs: Vec<_> = (0..workers).map(|t| { let (diamonds, blocks, diamond_of, faults, ops, txs) = (&diamonds, &blocks, &diamond_of, &faults, &ops, core_tx.clone());
            sc.spawn(move || for bi in t * chunk..((t + 1) * chunk).min(bw * bh) { let d = diamond_of[bi]; let (rr, w) = diamonds[d].block(&blocks[bi], ops, faults[bi]); txs[core_of(d)].send((bi, rr, w)).unwrap(); }) }).collect();
        drop(core_tx);
        // FINAL core: verify each sub-core's checksum, assemble the frame
        let mut frame = vec![[0i32; 64]; bw * bh]; let (mut corrected, mut recomputed, mut bad_sum) = (0usize, 0usize, 0usize);
        for (bi, o, cs, how) in final_rx { if checksum(&o) != cs { bad_sum += 1; } frame[bi] = o; if how == 1 { corrected += 1; } if how == 2 { recomputed += 1; } }
        for h in whs { h.join().unwrap(); } for h in core_hs { h.join().unwrap(); }
        (frame, corrected, recomputed, bad_sum)
    });
    let ta = t1.elapsed().as_secs_f64();
    let (frame, corrected, recomputed, bad_sum) = report;
    let wrong = frame.iter().zip(&plain).filter(|(a, b)| a != b).count();
    let injected = faults.iter().filter(|f| f.is_some()).count();
    let rings: Vec<String> = diamonds.iter().map(|d| format!("{:?}/{}", d.p, d.data)).collect();
    println!("{{\"probe\":\"armilla-cores-vp9-idct8x8\",\"layout\":\"{layout}\",\"cores\":{},\"workers\":{workers},\"blocks\":{},\"blocks_wrong\":{wrong},\"faults_injected\":{injected},\"corrected_in_place\":{corrected},\"recomputed\":{recomputed},\"final_checksum_failures\":{bad_sum},\"armilla_s\":{ta:.4},\"plain_s\":{tp:.4},\"ratio\":{:.1},\"frames_per_s\":{:.2},\"rings_per_diamond\":{:?}}}",
        ncore + 1, bw * bh, ta / tp, 1.0 / ta, rings);
}
