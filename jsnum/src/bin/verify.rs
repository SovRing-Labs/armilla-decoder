//! Bit-exact check of the Armilla JS-number unit against the CPU's IEEE 754 f64 (the same semantics as ECMA-262 Number), on random operands
//! across the whole exponent range plus edge cases; guard alarms on clean runs; fault injection; ns/op vs native.
use armilla_jsnum::*;
use std::time::Instant;

fn same(a: f64, b: f64) -> bool { (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits() }

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(2_000_000);
    let mut s: u64 = 0x9E3779B97F4A7C15;
    let mut nx = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    let edge = [0.0, -0.0, 1.0, -1.0, 0.5, 3.0, f64::MIN_POSITIVE, -f64::MIN_POSITIVE, 5e-324, -5e-324, f64::MAX, -f64::MAX, f64::INFINITY, f64::NEG_INFINITY,
        f64::NAN, 1e308, 1e-308, 2.0f64.powi(53), 2.0f64.powi(53) + 2.0, 0.1, 0.2, 0.3, 1.0 / 3.0, 2.2250738585072009e-308, 4.9406564584124654e-324 * 3.0];
    // random doubles: random bit patterns (whole range incl. subnormals), plus random "near" pairs (cancellation, ties)
    let mut xs = vec![]; for &a in &edge { for &b in &edge { xs.push((a, b)); } }
    while xs.len() < n {
        let a = f64::from_bits(nx()); let b = match nx() % 4 { 0 => f64::from_bits(nx()), 1 => f64::from_bits(a.to_bits() ^ (nx() & 0xFFFF)), 2 => -a * (1.0 + (nx() % 7) as f64 * f64::EPSILON), _ => ((nx() % 2001) as f64 - 1000.0) / ((nx() % 999) + 1) as f64 };
        let a = if nx() % 5 == 0 { ((nx() % 100001) as f64 - 50000.0) / 7.0 } else { a };
        xs.push((a, b));
    }
    let mut g = Guard::default();
    type F = fn(f64, f64, &mut Guard) -> f64;
    let ops: [(&str, F, fn(f64, f64) -> f64); 4] = [("add", add, |a, b| a + b), ("sub", sub, |a, b| a - b), ("mul", mul, |a, b| a * b), ("div", div, |a, b| a / b)];
    for (name, f, nat) in ops {
        let mut bad = 0usize; let mut ex = None;
        let t = Instant::now(); let mut acc = 0u64;
        for &(a, b) in &xs { let r = f(a, b, &mut g); acc ^= r.to_bits(); if !same(r, nat(a, b)) { bad += 1; if ex.is_none() { ex = Some((a, b, r, nat(a, b))); } } }
        let ta = t.elapsed().as_secs_f64();
        let t = Instant::now(); for &(a, b) in &xs { acc ^= std::hint::black_box(nat(std::hint::black_box(a), b)).to_bits(); } let tn = t.elapsed().as_secs_f64();
        println!("{{\"op\":\"{name}\",\"cases\":{},\"mismatches\":{bad},\"first_mismatch\":\"{:?}\",\"armilla_ns\":{:.1},\"native_ns\":{:.2},\"acc\":{acc}}}", xs.len(), ex, ta * 1e9 / xs.len() as f64, tn * 1e9 / xs.len() as f64);
    }
    // int32 + compare
    let i32ops: [(&str, fn(f64, f64) -> f64, fn(i32, i32) -> i32); 4] = [("and", band, |a, b| a & b), ("or", bor, |a, b| a | b), ("xor", bxor, |a, b| a ^ b), ("shl", shl, |a, b| a.wrapping_shl(b as u32 & 31))];
    let ints: Vec<(f64, f64)> = xs.iter().map(|&(a, b)| (if a.is_finite() && a.abs() < 1e12 { a } else { (nx() as i64 >> 20) as f64 }, (nx() % 40) as f64)).collect();
    let toi = |x: f64| -> i32 { if !x.is_finite() { 0 } else { (x.trunc().rem_euclid(4294967296.0) as u64 as u32) as i32 } };
    for (name, f, nat) in i32ops { let bad = ints.iter().filter(|&&(a, b)| f(a, b) != nat(toi(a), toi(b)) as f64).count(); println!("{{\"op\":\"{name}\",\"cases\":{},\"mismatches\":{bad}}}", ints.len()); }
    let bad = ints.iter().filter(|&&(a, b)| sar(a, b, &mut g) != (toi(a) >> (toi(b) as u32 & 31)) as f64).count(); println!("{{\"op\":\"sar\",\"cases\":{},\"mismatches\":{bad}}}", ints.len());
    let bad = ints.iter().filter(|&&(a, b)| shr(a, b) != ((toi(a) as u32) >> (toi(b) as u32 & 31)) as f64).count(); println!("{{\"op\":\"shr\",\"cases\":{},\"mismatches\":{bad}}}", ints.len());
    let bad = xs.iter().filter(|&&(a, b)| lt(a, b) != (a < b)).count(); println!("{{\"op\":\"lt\",\"cases\":{},\"mismatches\":{bad}}}", xs.len());
    println!("{{\"guard_checks_clean\":{},\"guard_alarms_clean\":{}}}", g.checks, g.alarms);
    // fault injection: one bit of the 2^128 value flipped after the integer core of a random op
    let mut gf = Guard::default(); let (mut changed, mut caught) = (0usize, 0usize);
    for (k, &(a, b)) in xs.iter().enumerate().take(500_000) {
        if !(a.is_finite() && b.is_finite() && a != 0.0 && b != 0.0) { continue; }
        let op = k % 3; let before = gf.alarms;
        FAULT.store((nx() % 110) as u32, std::sync::atomic::Ordering::Relaxed);
        let r = match op { 0 => add(a, b, &mut gf), 1 => mul(a, b, &mut gf), _ => div(a, b, &mut gf) };
        FAULT.store(u32::MAX, std::sync::atomic::Ordering::Relaxed);
        let want = match op { 0 => a + b, 1 => a * b, _ => a / b };
        let alarm = gf.alarms > before;
        if !same(r, want) { changed += 1; if alarm { caught += 1; } else if changed - caught <= 5 { eprintln!("UNCAUGHT op {op} a {a:e} b {b:e} got {r:e} want {want:e}"); } }
    }
    println!("{{\"faults_injected_runs\":\"≤500000 finite\",\"outputs_changed\":{changed},\"caught_by_routing_rings\":{caught}}}");
}
