//! J's Armilla layout (10-15, readback confirmed): one DIAMOND = 5 layers of T50 ring systems 1-3-5-3-1 (13 systems); a GROUP = 3 diamonds;
//! 3 groups (117 systems); every group computes its own share and outputs to the CENTRE core = the math core.
//! Workload: the VP9 8×8 2D inverse DCT of a 1080p luma frame (VP9 spec v0.6 §8.7.1.2–8.7.1.3 + §8.7.2: rows, columns, final
//! Round2(x, Min(6, n+2)) = Round2(x, 5)).
//!
//! Mapping (seat's, stated to J): system s holds every value as residues on its own 50 odd primes (system 0 = T50 3…233, then the next
//! primes, disjoint) plus its own 2^64 low-bit ring. Round2(x,14) inside a system = low 14 bits from its 2^64 ring + exact division on its
//! 50 primes (× inverse of 2^14; on the T50 routing primes 3·43·127 / 5·29·113 that inverse is +1 / −1 = a pure route). The 2D transform
//! never leaves the residues. Shares: frame → 3 bands (groups) → 3 sub-bands (diamonds) → blocks round-robin over the 13 systems in
//! layer order 1,3,5,3,1. Centre: per block, check all 50 residues of every output against the system's 2^64 value (the guard), apply
//! the final Round2(x,5), and recompute a block that fails (repair).
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

const K: usize = 50;
struct System { p: [u64; K], inv14: [u64; K], cos_r: Vec<[u64; K]>, sin_r: Vec<[u64; K]>, nsin_r: Vec<[u64; K]> }
impl System {
    fn new(p: [u64; K], ops: &[Op]) -> System {
        let cm = |v: i64, m: u64| v.rem_euclid(m as i64) as u64;
        let angs: Vec<i64> = ops.iter().map(|o| if let Op::B(_, _, a, _) = o { *a } else { 0 }).collect();
        System { p, inv14: p.map(|m| inv(16384 % m, m)),
            cos_r: angs.iter().map(|&a| p.map(|m| cm(cos64(a), m))).collect(), sin_r: angs.iter().map(|&a| p.map(|m| cm(sin64(a), m))).collect(),
            nsin_r: angs.iter().map(|&a| p.map(|m| cm(-sin64(a), m))).collect() }
    }
    /// one 1D inverse DCT on residues (r: 8 values × K rings) + the 2^64 carrier (w)
    fn idct(&self, r: &mut [[u64; K]; 8], w: &mut [u64; 8], ops: &[Op]) {
        let (cr, cw) = (*r, *w); for i in 0..8 { r[i] = cr[brev(3, i)]; w[i] = cw[brev(3, i)]; }
        for (oi, &op) in ops.iter().enumerate() { match op {
            Op::B(a, b, ang, flip) => {
                let (c, s) = (cos64(ang), sin64(ang));
                let xw = w[a].wrapping_mul(c as u64).wrapping_add(w[b].wrapping_mul((-s) as u64)).wrapping_add(8192);
                let yw = w[a].wrapping_mul(s as u64).wrapping_add(w[b].wrapping_mul(c as u64)).wrapping_add(8192);
                let (lx, ly) = (xw & 0x3FFF, yw & 0x3FFF);
                let (cs, sn, ns) = (&self.cos_r[oi], &self.sin_r[oi], &self.nsin_r[oi]);
                for k in 0..K {
                    let m = self.p[k];
                    let x = (r[a][k] * cs[k] + r[b][k] * ns[k] + 8192 + m * 16384 - lx) % m;
                    let y = (r[a][k] * sn[k] + r[b][k] * cs[k] + 8192 + m * 16384 - ly) % m;
                    r[a][k] = x * self.inv14[k] % m; r[b][k] = y * self.inv14[k] % m;
                }
                w[a] = ((xw as i64) >> 14) as u64; w[b] = ((yw as i64) >> 14) as u64;
                if flip { r.swap(a, b); w.swap(a, b); }
            }
            Op::H(a, b) => {
                for k in 0..K { let m = self.p[k]; let (x, y) = (r[a][k], r[b][k]); r[a][k] = (x + y) % m; r[b][k] = (x + m - y) % m; }
                let (x, y) = (w[a], w[b]); w[a] = x.wrapping_add(y); w[b] = x.wrapping_sub(y);
            }
        } }
    }
    /// the system's share: a whole 8×8 block in residues → (64 output residue rows, 64 carrier values) before the final Round2
    fn block(&self, blk: &[i32; 64], ops: &[Op], fault: Option<(usize, usize)>) -> (Vec<[u64; K]>, [i64; 64]) {
        let mut rr = vec![[0u64; K]; 64]; let mut ww = [0u64; 64];
        for i in 0..64 { let x = blk[i] as i64; for k in 0..K { rr[i][k] = x.rem_euclid(self.p[k] as i64) as u64; } ww[i] = x as u64; }
        for i in 0..8 { let mut r = [[0u64; K]; 8]; let mut w = [0u64; 8]; for j in 0..8 { r[j] = rr[i * 8 + j]; w[j] = ww[i * 8 + j]; }
            self.idct(&mut r, &mut w, ops); for j in 0..8 { rr[i * 8 + j] = r[j]; ww[i * 8 + j] = w[j]; } }
        if let Some((pos, k)) = fault { rr[pos][k] = (rr[pos][k] + 1) % self.p[k]; } // a corrupted residue between the passes
        for j in 0..8 { let mut r = [[0u64; K]; 8]; let mut w = [0u64; 8]; for i in 0..8 { r[i] = rr[i * 8 + j]; w[i] = ww[i * 8 + j]; }
            self.idct(&mut r, &mut w, ops); for i in 0..8 { rr[i * 8 + j] = r[i]; ww[i * 8 + j] = w[i]; } }
        let mut wi = [0i64; 64]; for i in 0..64 { wi[i] = ww[i] as i64; }
        (rr, wi)
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let fault_rate: f64 = a.get(1).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let threads: usize = a.get(2).and_then(|v| v.parse().ok()).unwrap_or(8);
    let mut ops = vec![]; idct_ops(3, &mut ops);
    // 117 systems × 50 disjoint odd primes; system 0 = T50 (3…233)
    let mut primes = vec![]; let mut n = 3u64;
    while primes.len() < 117 * K { if (2..).take_while(|d| d * d <= n).all(|d| n % d != 0) { primes.push(n); } n += 2; }
    let systems: Vec<System> = (0..117).map(|s| { let mut p = [0u64; K]; p.copy_from_slice(&primes[s * K..(s + 1) * K]); System::new(p, &ops) }).collect();
    // layout: group g (0..3), diamond d (0..3), layer sizes 1,3,5,3,1 → system index
    let layers = [1usize, 3, 5, 3, 1];
    // 1080p luma frame of 8×8 blocks with coefficients in the 8-bit-depth range (|c| ≤ 2^12)
    let (bw, bh) = (240usize, 135usize);
    let mut s: u64 = 0x9E3779B97F4A7C15;
    let mut rnd = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    let blocks: Vec<[i32; 64]> = (0..bw * bh).map(|_| { let mut b = [0i32; 64]; for (k, c) in b.iter_mut().enumerate() { *c = if k < 20 || rnd() % 4 == 0 { (rnd() % 1024) as i32 - 512 } else { 0 }; } b }).collect();
    let t0 = Instant::now();
    let plain: Vec<[i32; 64]> = blocks.iter().map(|b| plain_block(b, &ops)).collect();
    let tp = t0.elapsed().as_secs_f64();
    // shares: block index → (group, diamond, system) — bands by block row
    let assign: Vec<usize> = (0..bw * bh).map(|bi| {
        let row = bi / bw; let g = row * 3 / bh; let d = (row * 3 % bh) * 3 / bh % 3; let within = bi % 13;
        let mut acc = 0; let mut sys = 0; for (l, &w) in layers.iter().enumerate() { if within < acc + w { sys = acc + (within - acc); let _ = l; break; } acc += w; }
        g * 39 + d * 13 + sys }).collect();
    let mut faults: Vec<Option<(usize, usize)>> = vec![None; bw * bh];
    if fault_rate > 0.0 { for f in faults.iter_mut() { if (rnd() % 1_000_000) as f64 / 1e6 < fault_rate { *f = Some(((rnd() % 64) as usize, (rnd() % K as u64) as usize)); } } }
    let t1 = Instant::now();
    // the groups / diamonds / systems compute their shares in parallel (threads), the centre checks + rounds + repairs
    let chunk = (bw * bh + threads - 1) / threads;
    let results: Vec<([i32; 64], u8)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads).map(|t| { let (systems, blocks, assign, faults, ops) = (&systems, &blocks, &assign, &faults, &ops);
            sc.spawn(move || (t * chunk..((t + 1) * chunk).min(bw * bh)).map(|bi| {
                let sys = &systems[assign[bi]];
                let (rr, w) = sys.block(&blocks[bi], ops, faults[bi]);
                // CENTRE (math core): guard = every residue of every output equals the carrier value mod its prime
                let ok = (0..64).all(|i| (0..K).all(|k| w[i].rem_euclid(sys.p[k] as i64) as u64 == rr[i][k]));
                if ok { let mut o = [0i32; 64]; for i in 0..64 { o[i] = round2(w[i], 5) as i32; } (o, 0u8) }
                else { (plain_block(&blocks[bi], ops), 1u8) } // repair: recompute the block exactly at the centre
            }).collect::<Vec<_>>()) }).collect();
        hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    let ta = t1.elapsed().as_secs_f64();
    let mismatches = results.iter().zip(&plain).filter(|(r, p)| r.0 != **p).count();
    let repaired = results.iter().filter(|r| r.1 == 1).count();
    let injected = faults.iter().filter(|f| f.is_some()).count();
    let mut per_sys = vec![0usize; 117]; for &s in &assign { per_sys[s] += 1; }
    println!("{{\"probe\":\"armilla-grid-vp9-idct8x8\",\"layout\":\"3 groups x 3 diamonds x (1-3-5-3-1) = 117 T50-class systems x 50 primes + own 2^64 ring; centre = check + final Round2 + repair\",\"frame\":\"1080p luma, {} 8x8 blocks\",\"threads\":{threads},\"blocks_mismatching_plain\":{mismatches},\"faults_injected\":{injected},\"blocks_repaired_by_centre\":{repaired},\"plain_s\":{tp:.4},\"armilla_s\":{ta:.4},\"ratio\":{:.1},\"armilla_frames_per_s\":{:.2},\"blocks_per_system_min_max\":[{},{}],\"largest_prime\":{}}}",
        bw * bh, ta / tp, 1.0 / ta, per_sys.iter().min().unwrap(), per_sys.iter().max().unwrap(), primes.last().unwrap());
}
