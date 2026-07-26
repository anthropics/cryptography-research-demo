// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
mod aes;
mod bridge;
mod exp_basic;
mod exp_entropy;
mod exp_honest;
mod exp_timing;
mod exp_wrongkey;
mod gf;
mod rng;
mod sr744;

use std::env;
use std::fs;
use std::time::Instant;

fn stamp() -> String {
    let out = std::process::Command::new("date").args(["-u", "+%Y-%m-%d %H:%M:%S UTC"]).output();
    match out { Ok(o) => String::from_utf8_lossy(&o.stdout).trim().to_string(), Err(_) => "?".into() }
}

fn save(name: &str, body: &str, secs: f64) {
    let text = format!("[{}] wall-clock {:.1}s\n{}\n", stamp(), secs, body);
    print!("{}", text);
    fs::create_dir_all("results").ok();
    fs::write(format!("results/{}.txt", name), text).ok();
}

fn run<F: FnOnce() -> String>(name: &str, f: F) {
    let t0 = Instant::now();
    let out = f();
    save(name, &out, t0.elapsed().as_secs_f64());
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let which = args.get(1).map(|s| s.as_str()).unwrap_or("all");
    let n_arg: u64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
    let gf = gf::Gf::new();
    let all = which == "all";
    if all || which == "e1" { run("e1_bridge_identity", || exp_basic::e1(&gf)); }
    if all || which == "e2" { run("e2_sigma", || exp_basic::e2(&gf)); }
    if all || which == "e3" { run("e3_exponents", || exp_basic::e3(&gf)); }
    if all || which == "e8" { run("e8_ddt", || exp_basic::e8(&gf)); }
    if all || which == "e9" { run("e9_aes_sanity", || exp_basic::e9(&gf)); }
    if all || which == "e4" {
        let n = if n_arg > 0 { n_arg } else { 200_000 };
        run("e4_honest_bridge", || exp_honest::e4(&gf, n));
    }
    if all || which == "e7" {
        let n = if n_arg > 0 { n_arg } else { 5_000_000 };
        run("e7_agl_invariance", || exp_entropy::e7(&gf, n));
    }
    if which == "e5" {
        let lg = if n_arg > 0 { n_arg as u32 } else { 32 };
        run(&format!("e5_wrongkey_2^{}", lg), || exp_wrongkey::e5(&gf, lg));
    }
    if which == "e6" {
        let n = if n_arg > 0 { n_arg } else { 100_000_000 };
        let pt = bridge::PowTab::new(&gf);
        let sub = args.get(3).map(|s| s.as_str()).unwrap_or("all");
        let off: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
        let order: Vec<usize> = args.get(5).map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_default();
        let ord_tag = if order.is_empty() { String::new() } else { format!("_ord{}", args.get(5).unwrap()) };
        let l255 = 255f64.log2();
        if sub == "all" || sub == "imn_odd" {
            run(&format!("e6_imn_odd_{}_off{}{}", n, off, ord_tag), || {
                let r = exp_entropy::measure(&gf, &pt, exp_entropy::Kind::Imn, false, n, 1, off, &order);
                exp_entropy::e6_report(&format!("E6a I_(m,n) 12-byte fingerprint, |D|=255 (odd), bytes {}..", off+1), &r, l255, "log2 255, bytes are nonzero")
            });
        }
        if sub == "all" || sub == "imn_even" {
            run(&format!("e6_imn_even_{}_off{}{}", n, off, ord_tag), || {
                let r = exp_entropy::measure(&gf, &pt, exp_entropy::Kind::Imn, true, n, 2, off, &order);
                exp_entropy::e6_report(&format!("E6b I_(m,n) 12-byte fingerprint, |D|=256 (appended zero, even), bytes {}..", off+1), &r, l255, "log2 255")
            });
        }
        let p_odd = 0.5 * (1.0 - (127.0f64 / 128.0).powi(255));
        let cbits = -(p_odd * p_odd + (1.0 - p_odd) * (1.0 - p_odd)).log2() * 8.0;
        if sub == "all" || sub == "chi_odd" {
            run(&format!("e6_chi_odd_{}_off{}{}", n, off, ord_tag), || {
                let r = exp_entropy::measure(&gf, &pt, exp_entropy::Kind::ChiStar, false, n, 3, off, &order);
                exp_entropy::e6_report(&format!("E6c chi-star canonical vector, |D|=255 (odd), bytes {}..", off+1), &r, cbits, "8 x -log2(p^2+(1-p)^2), p = 0.4323 -> 7.79, the paper's independence assumption")
            });
        }
        if sub == "all" || sub == "chi_even" {
            run(&format!("e6_chi_even_{}_off{}{}", n, off, ord_tag), || {
                let r = exp_entropy::measure(&gf, &pt, exp_entropy::Kind::ChiStar, true, n, 4, off, &order);
                exp_entropy::e6_report(&format!("E6d chi-star canonical vector, |D|=256 (appended zero, even), bytes {}..", off+1), &r, cbits, "7.79 per byte if independent")
            });
        }
    }
    if all || which == "e10" { run("e10_timing", || exp_timing::e10(&gf)); }
    if which == "e4d" { run("e4_debug", || exp_honest::e4_debug(&gf)); }
    if which == "e12" { let lg = if n_arg > 0 { n_arg as u32 } else { 28 }; run(&format!("e12_fp_direct_2^{}", lg), || exp_entropy::e12(&gf, lg)); }
    if which.starts_with("e13") {
        let g16 = sr744::Gf16::new();
        if which == "e13a" { run("e13a_sr744_sanity", || sr744::e13a(&g16)); }
        if which == "e13s" { let n = if n_arg > 0 { n_arg } else { 8 }; run("e13s_sizing", || sr744::e13s(&g16, n)); }
        if which == "e13b" { let n = if n_arg > 0 { n_arg } else { 200_000 }; run("e13b_sr744_bridge", || sr744::e13b(&g16, n)); }
        if which == "e13d" { run("e13d_stages", || sr744::e13d(&g16)); }
        if which == "e13e" { let n = if n_arg > 0 { n_arg } else { 3 }; let k: usize = args.get(3).and_then(|x| x.parse().ok()).unwrap_or(8); run(&format!("e13e_honest_{}keys_K{}", n, k), || sr744::e13e_sweep(&g16, n, k, 0xE13E)); }
        if which == "e13c" { run("e13c_sr744_keyrecovery", || sr744::e13c(&g16, n_arg)); }
    }
    if which == "e5c" { let lg = if n_arg > 0 { n_arg as u32 } else { 32 }; run(&format!("e5c_chi_wrongkey_2^{}", lg), || exp_wrongkey::e5_chi(&gf, lg)); }
}
