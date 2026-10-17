//! Armilla JS-number unit (J 10-15, readbacks confirmed: "the Armilla should do any kind of arithmetic … it has the capability, it's now just a
//! matter of organizing the rings, coprimes, how they interact within each unit and group, and the routing"; "the math core should also have
//! algorithms it can use for any kind of computation").
//!
//! Ring organisation of one unit:
//! - the POWER-OF-TWO ring Z/2^128 carries the exact value: rounding, normalisation (leading bit), bitwise ops and shifts read its low bits
//!   directly — no bridge;
//! - the ROUTING rings 3·43·127 (2^14 ≡ +1) and 5·29·113 (2^14 ≡ −1) — the factors of 2^28 − 1 — carry the same integers as residues; on them
//!   × or ÷ by any power of two is a rotation of a 14-cycle (identity / negation for 2^14), so every alignment, shift and rounding split is routing;
//! - the routing rings are the GUARD: each operation's exact integer core (aligned sum, product, division with remainder, the rounding split
//!   M = q·2^s + rem) is checked residue by residue; a corrupted 2^128 value disagrees (miss probability ≤ 1/(2^28 − 1) per check);
//! - division is a MATH-CORE ALGORITHM on the 2^128 ring (long division), verified on the routing rings by q·d + r = n.
//!
//! Semantics: ECMA-262 §6.1.6.1 Number (IEEE 754-2019 binary64, roundTiesToEven, −0, ±∞, NaN) and §6.1.6.1.x int32 operations
//! (ToInt32 wrap, &, |, ^, <<, >>, >>>). Clean room: from IEEE 754 / ECMA-262 definitions; the reference for testing is the CPU's own f64.

pub const R: [u64; 6] = [3, 43, 127, 5, 29, 113];
/// 2^k mod p for k = 0..13 on each routing ring (the 14-cycle; 2^14 ≡ ±1)
fn pow2_mod(k: u64, p: u64) -> u64 { let mut v = 1u64; let s = if [3u64, 43, 127].contains(&p) { 1 } else { p - 1 }; for _ in 0..(k % 14) { v = v * 2 % p; } let t = (k / 14) % 2; if t == 1 && s != 1 { (p - v) % p } else { v } }

/// residues of an unsigned 128-bit integer on the routing rings
#[inline] pub fn res(x: u128) -> [u64; 6] { R.map(|p| (x % p as u128) as u64) }
#[inline] fn shl_r(r: [u64; 6], k: u64) -> [u64; 6] { let mut o = r; for i in 0..6 { o[i] = r[i] * pow2_mod(k, R[i]) % R[i]; } o }
#[inline] fn add_r(a: [u64; 6], b: [u64; 6]) -> [u64; 6] { let mut o = a; for i in 0..6 { o[i] = (a[i] + b[i]) % R[i]; } o }
#[inline] fn sub_r(a: [u64; 6], b: [u64; 6]) -> [u64; 6] { let mut o = a; for i in 0..6 { o[i] = (a[i] + R[i] - b[i]) % R[i]; } o }
#[inline] fn mul_r(a: [u64; 6], b: [u64; 6]) -> [u64; 6] { let mut o = a; for i in 0..6 { o[i] = a[i] * b[i] % R[i]; } o }

/// Mersenne guard ring M61 = 2^61 − 1 (J 10-15 ring organisation): routing rings alone are blind when an operand is a multiple of 2^28 − 1
/// (repeating binary fractions such as 4/17 have such mantissas). A mantissa (< 2^53, ≠ 0) is never ≡ 0 mod M61 and 2 is invertible, so any
/// single-bit fault 2^b·x is visible; ×2^k on M61 is a 61-bit rotation (2^61 ≡ 1) — routing again.
pub const M61: u128 = (1u128 << 61) - 1;
#[inline] fn m61(x: u128) -> u128 { let mut v = (x & M61) + ((x >> 61) & M61) + (x >> 122); while v >= M61 { v -= M61; } v }
#[inline] fn m61_shl(x: u128, k: u64) -> u128 { let k = (k % 61) as u32; if k == 0 { x } else { m61(((x << k) & M61) | (x >> (61 - k))) } }
#[inline] fn m61_mul(a: u128, b: u128) -> u128 { m61(a * b) }
/// a guard failure: the 2^128 ring and the routing rings disagree
#[derive(Debug, Default, Clone, Copy)]
pub struct Guard { pub checks: u64, pub alarms: u64 }
impl Guard {
    fn check(&mut self, ring: u128, want: [u64; 6]) { self.checks += 1; if res(ring) != want { self.alarms += 1; } }
    fn check61(&mut self, ring: u128, want: u128) { self.checks += 1; if m61(ring) != want { self.alarms += 1; } }
}

/// fault model for tests: flip bit `bit` of the 2^128 value right after the integer core of the next operation
pub static FAULT: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(u32::MAX);
fn maybe_fault(m: &mut u128) { let b = FAULT.swap(u32::MAX, std::sync::atomic::Ordering::Relaxed); if b < 128 { *m ^= 1u128 << b; } }

/// finite double → (sign, integer mantissa, exponent): value = (−1)^s · m · 2^e
fn decode(x: f64) -> (bool, u128, i32) {
    let b = x.to_bits(); let s = b >> 63 == 1; let ex = ((b >> 52) & 0x7FF) as i32; let fr = (b & ((1u64 << 52) - 1)) as u128;
    if ex == 0 { (s, fr, -1074) } else { (s, fr | (1u128 << 52), ex - 1075) }
}

/// round (−1)^s · M · 2^e (+ a sticky "something below") to binary64, roundTiesToEven, with gradual underflow and overflow to ±∞.
/// The split M = q·2^shift + rem is guarded on the routing rings.
fn round(s: bool, mut m: u128, e: i32, mut sticky: bool, g: &mut Guard) -> f64 {
    let sgn = if s { -0.0f64 } else { 0.0f64 };
    if m == 0 { return if sticky { f64::from_bits(sgn.to_bits()) } else { sgn }; }
    let bl = 128 - m.leading_zeros() as i32;
    let mut shift = (bl - 53).max(-1074 - e);
    // deeper than the value's own length: everything is below half the smallest subnormal → ±0 (shift = bl rounds via the split below)
    if shift > bl { return sgn; }
    let mut q: u128;
    if shift > 0 {
        let sh = shift as u32;
        q = m >> sh; let rem = m & ((1u128 << sh) - 1);
        // guard: M = q·2^shift + rem on the routing rings (the shift is a 14-cycle rotation there)
        g.check(m, add_r(shl_r(res(q), shift as u64), res(rem)));
        let half = 1u128 << (sh - 1);
        if rem > half || (rem == half && (sticky || q & 1 == 1)) { q += 1; }
        sticky = false;
        if q == 1u128 << 53 { q >>= 1; shift += 1; }
    } else { q = m << (-shift) as u32; m = q; let _ = m; }
    let _ = sticky;
    let fe = e + shift;
    if fe > 971 { return if s { f64::NEG_INFINITY } else { f64::INFINITY }; }
    let bits = if q < (1u128 << 52) { q as u64 } else { (((fe + 1075) as u64) << 52) | (q as u64 & ((1u64 << 52) - 1)) };
    f64::from_bits(bits | ((s as u64) << 63))
}

/// ECMA-262 Number::add (IEEE addition, roundTiesToEven)
pub fn add(x: f64, y: f64, g: &mut Guard) -> f64 {
    if x.is_nan() || y.is_nan() { return f64::NAN; }
    if x.is_infinite() || y.is_infinite() { return if x.is_infinite() && y.is_infinite() && x.signum() != y.signum() { f64::NAN } else if x.is_infinite() { x } else { y }; }
    if x == 0.0 && y == 0.0 { return if x.is_sign_negative() && y.is_sign_negative() { -0.0 } else { 0.0 }; }
    if x == 0.0 { return y; } if y == 0.0 { return x; }
    let (s1, m1, e1) = decode(x); let (s2, m2, e2) = decode(y);
    // big = the operand with the larger exponent
    let ((sb, mb, eb), (ss, ms, es)) = if e1 >= e2 { ((s1, m1, e1), (s2, m2, e2)) } else { ((s2, m2, e2), (s1, m1, e1)) };
    let d = (eb - es) as u32;
    if d > 60 {
        // the small operand is below 1/16 of the big one's extended last place: 3 guard bits + ±1 represent it exactly for rounding
        let mut m = mb << 3; let mr = shl_r(res(mb), 3);
        let (sum, want) = if sb == ss { (m + 1, add_r(mr, res(1))) } else { (m - 1, sub_r(mr, res(1))) };
        m = sum; maybe_fault(&mut m); g.check(m, want);
        g.check61(m, m61(if sb == ss { m61_shl(m61(mb), 3) + 1 } else { m61_shl(m61(mb), 3) + M61 - 1 }));
        return round(sb, m, eb - 3, false, g);
    }
    let a = mb << d; let ar = shl_r(res(mb), d as u64); let br = res(ms);
    let (s, mut m, want) = if sb == ss { (sb, a + ms, add_r(ar, br)) }
        else if a >= ms { (sb, a - ms, sub_r(ar, br)) } else { (ss, ms - a, sub_r(br, ar)) };
    maybe_fault(&mut m); g.check(m, want);
    let (a61, b61) = (m61_shl(m61(mb), d as u64), m61(ms));
    g.check61(m, if sb == ss { m61(a61 + b61) } else if a >= ms { m61(a61 + M61 - b61) } else { m61(b61 + M61 - a61) });
    if m == 0 { return 0.0; } // x + (−x) = +0 under roundTiesToEven
    round(s, m, es, false, g)
}
pub fn sub(x: f64, y: f64, g: &mut Guard) -> f64 { add(x, -y, g) }

/// Number::multiply
pub fn mul(x: f64, y: f64, g: &mut Guard) -> f64 {
    if x.is_nan() || y.is_nan() { return f64::NAN; }
    let s = x.is_sign_negative() ^ y.is_sign_negative();
    if x.is_infinite() || y.is_infinite() { return if x == 0.0 || y == 0.0 { f64::NAN } else if s { f64::NEG_INFINITY } else { f64::INFINITY }; }
    if x == 0.0 || y == 0.0 { return if s { -0.0 } else { 0.0 }; }
    let (_, m1, e1) = decode(x); let (_, m2, e2) = decode(y);
    let mut m = m1 * m2; maybe_fault(&mut m); g.check(m, mul_r(res(m1), res(m2))); g.check61(m, m61_mul(m61(m1), m61(m2)));
    round(s, m, e1 + e2, false, g)
}

/// Number::divide — math-core algorithm: long division on the 2^128 ring, verified q·d + r = n on the routing rings
pub fn div(x: f64, y: f64, g: &mut Guard) -> f64 {
    if x.is_nan() || y.is_nan() { return f64::NAN; }
    let s = x.is_sign_negative() ^ y.is_sign_negative();
    if x.is_infinite() { return if y.is_infinite() { f64::NAN } else if s { f64::NEG_INFINITY } else { f64::INFINITY }; }
    if y.is_infinite() { return if s { -0.0 } else { 0.0 }; }
    if y == 0.0 { return if x == 0.0 { f64::NAN } else if s { f64::NEG_INFINITY } else { f64::INFINITY }; }
    if x == 0.0 { return if s { -0.0 } else { 0.0 }; }
    let (_, m1, e1) = decode(x); let (_, m2, e2) = decode(y);
    // normalise both to 53 bits (subnormals) so the quotient always has ≥ 55 significant bits
    let (n1, n2) = (m1.leading_zeros() as i32 - 75, m2.leading_zeros() as i32 - 75);
    let (m1, e1) = (m1 << n1 as u32, e1 - n1); let (m2, e2) = (m2 << n2 as u32, e2 - n2);
    let k = 74u32; let n = m1 << k;
    let mut q = n / m2; let r = n % m2;
    maybe_fault(&mut q);
    g.check(n, add_r(mul_r(res(q), res(m2)), res(r))); g.check61(n, m61(m61_mul(m61(q), m61(m2)) + m61(r)));
    round(s, q, e1 - e2 - k as i32, r != 0, g)
}

/// ECMA-262 ToInt32 (§7.1.6): the 2^32 ring of the value's integer part
pub fn to_int32(x: f64) -> i32 {
    if !x.is_finite() || x == 0.0 { return 0; }
    let t = x.trunc(); let (s, m, e) = decode(t);
    let v: u32 = if e >= 32 { 0 } else if e >= 0 { (m << e as u32) as u32 } else { (m >> (-e) as u32) as u32 };
    (if s { v.wrapping_neg() } else { v }) as i32
}
/// int32 ops on the power-of-two ring (bitwise lives there); shifts guarded on the routing rings where they are ring maps
pub fn band(a: f64, b: f64) -> f64 { (to_int32(a) & to_int32(b)) as f64 }
pub fn bor(a: f64, b: f64) -> f64 { (to_int32(a) | to_int32(b)) as f64 }
pub fn bxor(a: f64, b: f64) -> f64 { (to_int32(a) ^ to_int32(b)) as f64 }
pub fn shl(a: f64, b: f64) -> f64 { (to_int32(a).wrapping_shl(to_int32(b) as u32 & 31)) as f64 }
pub fn sar(a: f64, b: f64, g: &mut Guard) -> f64 {
    let (v, k) = (to_int32(a), to_int32(b) as u32 & 31); let q = v >> k;
    // guard: v = q·2^k + low (floor division) — on the routing rings, with v and q taken in a non-negative window (+2^31 offset scaled)
    let low = (v as i64 - ((q as i64) << k)) as u128; let off = 1u128 << 40;
    g.check(((v as i64 as i128) + (off << k) as i128) as u128, add_r(shl_r(res(((q as i64 as i128) + off as i128) as u128), k as u64), res(low)));
    q as f64
}
pub fn shr(a: f64, b: f64) -> f64 { ((to_int32(a) as u32) >> (to_int32(b) as u32 & 31)) as f64 }
/// Number::lessThan on the power-of-two ring (sign of the difference)
pub fn lt(a: f64, b: f64) -> bool { if a.is_nan() || b.is_nan() { false } else { let (x, y) = (a.to_bits() as i64, b.to_bits() as i64); let k = |v: i64| if v < 0 { i64::MIN - v } else { v }; if a == 0.0 && b == 0.0 { false } else { k(x) < k(y) } } }
