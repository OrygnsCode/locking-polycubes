// polycubes that lock with a rotated copy of themselves. this is the Rust one of the two
// programs, its result lines are labelled "B".
//
// modes:
//   selftest <nmax> [wupto]                      one line per n: counts and checksums
//   count <nmax> [wupto]                         same walk, counts in named fields
//   report <nmax> <wupto> <outdir>               both, plus JSON totals and the samples
//   shardtable <d> <outfile>                     the shard table, one REF shape per line
//   shard <nmax> <wupto> <table> <id> <outdir>   one shard of a sharded run
//   classify <file>                              reads shapes, prints a verdict for each
// wupto is the largest n where witnesses get counted (that needs all 24 matrices, so it's only
// for small n).
//
// TWINB_ABLATE=0 skips the locking test and =1 builds D_R without scanning it. both print WRONG
// counts on purpose (and say so on stderr), they're only for timing the two halves. report and
// shard mode refuse anything but the default 2, and the value goes in every JSON line.
//
// build: cargo build --release
mod dr;

use dr::{Cell, Mat, Work, NEIGH};
use std::env;
use std::io::Write;

const IDENTITY: Mat = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// dr.rs packs D_R at stride 34 per axis and D_R spans at most 2n - 1, so 17 is the limit.
const NMAX_SUPPORTED: usize = 17;

// cpu time

/// cpu seconds (user + kernel) from GetProcessTimes. cpu time because it adds up across shards.
#[cfg(windows)]
mod cputime {
    #[repr(C)]
    #[derive(Default, Clone, Copy)]
    pub struct FileTime { low: u32, high: u32 }

    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> isize;
        fn GetProcessTimes(h: isize, creation: *mut FileTime, exit: *mut FileTime,
                           kernel: *mut FileTime, user: *mut FileTime) -> i32;
    }

    pub fn cpu_seconds() -> f64 {
        let mut c = FileTime::default();
        let mut e = FileTime::default();
        let mut k = FileTime::default();
        let mut u = FileTime::default();
        let ok = unsafe {
            GetProcessTimes(GetCurrentProcess(), &mut c, &mut e, &mut k, &mut u)
        };
        if ok == 0 { return -1.0; }
        // FILETIME is in 100 ns units
        let to = |f: FileTime| (((f.high as u64) << 32) | f.low as u64) as f64 * 1e-7;
        to(k) + to(u)
    }
}

#[cfg(not(windows))]
mod cputime {
    pub fn cpu_seconds() -> f64 { -1.0 }
}

// sha256, hand-rolled so the crate has no dependencies. checked against sha256sum.

mod sha256 {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];

    pub fn hex(data: &[u8]) -> String {
        let mut h: [u32; 8] = [
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
        ];
        let bitlen = (data.len() as u64).wrapping_mul(8);
        let mut msg = data.to_vec();
        msg.push(0x80);
        while msg.len() % 64 != 56 { msg.push(0); }
        msg.extend_from_slice(&bitlen.to_be_bytes());
        for chunk in msg.chunks(64) {
            let mut w = [0u32; 64];
            for i in 0..16 {
                w[i] = u32::from_be_bytes([chunk[4*i], chunk[4*i+1], chunk[4*i+2], chunk[4*i+3]]);
            }
            for i in 16..64 {
                let a = w[i-15];
                let b = w[i-2];
                let s0 = a.rotate_right(7) ^ a.rotate_right(18) ^ (a >> 3);
                let s1 = b.rotate_right(17) ^ b.rotate_right(19) ^ (b >> 10);
                w[i] = w[i-16].wrapping_add(s0).wrapping_add(w[i-7]).wrapping_add(s1);
            }
            let mut v = h;
            for i in 0..64 {
                let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
                let ch = (v[4] & v[5]) ^ ((!v[4]) & v[6]);
                let t1 = v[7].wrapping_add(s1).wrapping_add(ch)
                             .wrapping_add(K[i]).wrapping_add(w[i]);
                let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
                let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
                let t2 = s0.wrapping_add(maj);
                v[7] = v[6]; v[6] = v[5]; v[5] = v[4];
                v[4] = v[3].wrapping_add(t1);
                v[3] = v[2]; v[2] = v[1]; v[1] = v[0];
                v[0] = t1.wrapping_add(t2);
            }
            for i in 0..8 { h[i] = h[i].wrapping_add(v[i]); }
        }
        let mut s = String::with_capacity(64);
        for x in h.iter() { s.push_str(&format!("{:08x}", x)); }
        s
    }
}

fn own_sha256() -> String {
    match env::current_exe().and_then(std::fs::read) {
        Ok(bytes) => sha256::hex(&bytes),
        Err(e) => format!("UNREADABLE({})", e),
    }
}

// symmetry

fn all_mats() -> Vec<Mat> {
    let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let (mut proper, mut improper) = (Vec::new(), Vec::new());
    for p in perms.iter() {
        let mut swaps = 0;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if p[i] > p[j] { swaps += 1; }
            }
        }
        let psign: i32 = if swaps % 2 == 0 { 1 } else { -1 };
        for sx in [1i32, -1] { for sy in [1i32, -1] { for sz in [1i32, -1] {
            let s = [sx, sy, sz];
            let mut m: Mat = [[0; 3]; 3];
            for i in 0..3 { m[i][p[i]] = s[i]; }
            if psign * sx * sy * sz == 1 { proper.push(m) } else { improper.push(m) }
        }}}
    }
    assert_eq!(proper.len(), 24);
    assert_eq!(improper.len(), 24);
    proper.extend(improper);
    proper
}

fn det(m: &Mat) -> i32 {
    m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
  - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0])
  + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0])
}

/// for a signed permutation matrix the inverse is the transpose
fn transpose(m: &Mat) -> Mat {
    let mut t: Mat = [[0; 3]; 3];
    for i in 0..3 { for j in 0..3 { t[i][j] = m[j][i]; } }
    t
}

/// one matrix from each {R, R^-1} pair, in input order so the identity stays at index 0
/// (decide() reads the identity witness from there).
///
/// fine for verdicts: (P, R, v) locks iff (P, R^-1, -R^-1 v) does. map both copies by
/// x -> R^-1(x - v) and they swap places, and overlap, contact and blocking along the six axes
/// all survive a signed permutation (blocking is symmetric too: Q + t*d meets P iff P - t*d
/// meets Q). NOT fine for witness counts, each pair stands for two (R, v), so debug mode turns
/// it off.
fn inverse_pair_reps(ms: &[Mat]) -> Vec<Mat> {
    let mut out = Vec::new();
    for (i, m) in ms.iter().enumerate() {
        let inv = transpose(m);
        let j = ms.iter().position(|x| *x == inv)
            .expect("matrix set is not closed under inverse");
        if j >= i { out.push(*m); }
    }
    out
}

#[inline]
fn apply(m: &Mat, c: &Cell) -> Cell {
    [m[0][0]*c[0] + m[0][1]*c[1] + m[0][2]*c[2],
     m[1][0]*c[0] + m[1][1]*c[1] + m[1][2]*c[2],
     m[2][0]*c[0] + m[2][1]*c[1] + m[2][2]*c[2]]
}

/// translate to the min corner and sort, also returns the shift
fn normalize_with_shift(cells: &[Cell]) -> (Vec<Cell>, Cell) {
    let mut mn = [i32::MAX; 3];
    for c in cells {
        for k in 0..3 { if c[k] < mn[k] { mn[k] = c[k]; } }
    }
    let mut out: Vec<Cell> = cells.iter()
        .map(|c| [c[0]-mn[0], c[1]-mn[1], c[2]-mn[2]]).collect();
    out.sort_unstable();
    (out, mn)
}

fn normalize(cells: &[Cell]) -> Vec<Cell> { normalize_with_shift(cells).0 }

/// REF is the lex-least of the 48 normalized sorted images, plus |Aut| (images equal to it)
/// and whether it's achiral (some improper matrix fixes it).
fn ref_aut_achiral(cells: &[Cell], mats: &[Mat]) -> (Vec<Cell>, u32, bool) {
    let (rf, aut, ach, _, _) = ref_full(cells, mats);
    (rf, aut, ach)
}

/// like ref_aut_achiral, but also returns the first index into mats whose image is REF, and
/// that image's min corner. any matrix that maps to REF would do: two of them differ by a
/// stabiliser element of REF, so the mapped cell lands in the same orbit, and acceptance is
/// orbit membership.
fn ref_full(cells: &[Cell], mats: &[Mat]) -> (Vec<Cell>, u32, bool, usize, Cell) {
    let base = normalize(cells);
    let mut imgs: Vec<Vec<Cell>> = Vec::with_capacity(mats.len());
    let mut shifts: Vec<Cell> = Vec::with_capacity(mats.len());
    let mut achiral = false;
    for (i, m) in mats.iter().enumerate() {
        let t: Vec<Cell> = cells.iter().map(|c| apply(m, c)).collect();
        let (img, mn) = normalize_with_shift(&t);
        if i >= 24 && img == base { achiral = true; }
        imgs.push(img);
        shifts.push(mn);
    }
    let mut bi = 0usize;
    for i in 1..imgs.len() {
        if imgs[i] < imgs[bi] { bi = i; }
    }
    let aut = imgs.iter().filter(|x| **x == imgs[bi]).count() as u32;
    let best = imgs[bi].clone();
    let shift = shifts[bi];
    (best, aut, achiral, bi, shift)
}

// hashing

#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

const SEED1: u64 = 0x243F6A8885A308D3;
const SEED2: u64 = 0x13198A2E03707344;

fn hash_ref(rf: &[Cell], seed: u64) -> u64 {
    let mut h = seed;
    for c in rf {
        h = splitmix64(h ^ ((c[0] as u64) | ((c[1] as u64) << 8) | ((c[2] as u64) << 16)));
    }
    splitmix64(h ^ (rf.len() as u64))
}

// enumeration

fn in_set(sorted: &[Cell], c: &Cell) -> bool { sorted.binary_search(c).is_ok() }

fn stabiliser_maps(rf: &[Cell], mats: &[Mat]) -> Vec<(Mat, Cell)> {
    let mut out = Vec::new();
    for m in mats {
        let t: Vec<Cell> = rf.iter().map(|c| apply(m, c)).collect();
        let (img, mn) = normalize_with_shift(&t);
        if img == *rf { out.push((*m, [-mn[0], -mn[1], -mn[2]])); }
    }
    out
}

fn orbit_reps(cells: &[Cell], maps: &[(Mat, Cell)]) -> Vec<Cell> {
    let mut seen: Vec<Cell> = Vec::new();
    let mut reps = Vec::new();
    for c in cells {
        if seen.binary_search(c).is_ok() { continue; }
        reps.push(*c);
        for (m, t) in maps {
            let x = apply(m, c);
            let y = [x[0]+t[0], x[1]+t[1], x[2]+t[2]];
            if let Err(pos) = seen.binary_search(&y) { seen.insert(pos, y); }
        }
        if let Err(pos) = seen.binary_search(c) { seen.insert(pos, *c); }
    }
    reps
}

fn connected_without(rf: &[Cell], skip: usize) -> bool {
    let n = rf.len();
    if n <= 2 { return true; }
    let mut keep: Vec<Cell> = Vec::with_capacity(n - 1);
    for (i, c) in rf.iter().enumerate() { if i != skip { keep.push(*c); } }
    let mut seen = vec![false; keep.len()];
    let mut stack = vec![0usize];
    seen[0] = true;
    let mut cnt = 1;
    while let Some(i) = stack.pop() {
        for d in NEIGH.iter() {
            let t = [keep[i][0]+d[0], keep[i][1]+d[1], keep[i][2]+d[2]];
            if let Ok(j) = keep.binary_search(&t) {
                if !seen[j] { seen[j] = true; cnt += 1; stack.push(j); }
            }
        }
    }
    cnt == keep.len()
}

/// the canonical deletion orbit: the orbit of the lex-greatest cell whose removal keeps the
/// shape connected. REF is canonical, so it's the same for every representative, which is
/// what makes the subtrees disjoint and sharding exact.
fn canonical_orbit(rf: &[Cell], mats: &[Mat]) -> Vec<Cell> {
    // rf is sorted, so scanning down from the top and stopping at the first removable cell
    // gives the lex-greatest one, usually after one or two connectivity checks
    let mut pick: Option<Cell> = None;
    for i in (0..rf.len()).rev() {
        if connected_without(rf, i) { pick = Some(rf[i]); break; }
    }
    let pick = pick.unwrap();
    let maps = stabiliser_maps(rf, mats);
    let mut orbit: Vec<Cell> = Vec::new();
    for (m, t) in &maps {
        let x = apply(m, &pick);
        let y = [x[0]+t[0], x[1]+t[1], x[2]+t[2]];
        if let Err(pos) = orbit.binary_search(&y) { orbit.insert(pos, y); }
    }
    if let Err(pos) = orbit.binary_search(&pick) { orbit.insert(pos, pick); }
    orbit
}

struct Acc {
    all: u64,
    all_h: (u64, u64),
    twin: u64,
    twin_h: (u64, u64),
    mirror: u64,
    mirror_h: (u64, u64),
    wcount: u64,
    wh1: u64,
    fixed_equiv: u64,
    achiral: u64,
    translate: u64,
    sampled: u64,
}

impl Acc {
    fn new() -> Acc {
        Acc { all: 0, all_h: (0,0), twin: 0, twin_h: (0,0), mirror: 0, mirror_h: (0,0),
              wcount: 0, wh1: 0, fixed_equiv: 0, achiral: 0, translate: 0, sampled: 0 }
    }
}

/// the matrix lists: `proper` and `improper` are all 24 each (witness counting, debug mode),
/// `proper_red` and `improper_red` are the 17 {R, R^-1} reps used for verdicts.
struct Mats {
    all: Vec<Mat>,
    proper: Vec<Mat>,
    improper: Vec<Mat>,
    proper_red: Vec<Mat>,
    improper_red: Vec<Mat>,
}

impl Mats {
    fn build() -> Mats {
        let all = all_mats();
        let proper: Vec<Mat> = all[..24].to_vec();
        let improper: Vec<Mat> = all[24..].to_vec();

        // decide() expects the identity at index 0
        assert_eq!(all[0], IDENTITY, "mats[0] must be the identity");
        assert_eq!(proper[0], IDENTITY, "proper[0] must be the identity");
        assert!(proper.iter().all(|m| det(m) == 1), "proper coset has a det -1 matrix");
        assert!(improper.iter().all(|m| det(m) == -1), "improper coset has a det +1 matrix");

        let proper_red = inverse_pair_reps(&proper);
        let improper_red = inverse_pair_reps(&improper);
        // 10 self-inverse plus 7 pairs in each coset. an improper M is -R for a proper R and
        // M^2 = R^2, so M is self-inverse exactly when R is.
        assert_eq!(proper_red.len(), 17, "proper inverse-pair reps should be 17");
        assert_eq!(improper_red.len(), 17, "improper inverse-pair reps should be 17");
        assert_eq!(proper_red[0], IDENTITY, "the identity must stay at index 0");

        Mats { all, proper, improper, proper_red, improper_red }
    }
}

struct Ctx {
    ablate: u8,
    nmin: usize,
    nmax: usize,
    witness_upto: usize,
    acc: Vec<Acc>,
    seen: u64,
    work: Work,
    samples: Option<std::io::BufWriter<std::fs::File>>,
    /// shardtable mode: collect the REF shapes of exactly this size
    emit_at: usize,
    emitted: Vec<String>,
}

/// everything decided about one shape
struct Decision {
    twin: bool,
    mirror: bool,
    w: u64,
    by_identity: bool,
    witness: Option<(Mat, Cell)>,
}

impl Decision {
    fn none() -> Decision {
        Decision { twin: false, mirror: false, w: 0, by_identity: false, witness: None }
    }
}

/// mirror = twin, or some improper matrix locks, so a shape that locks by rotation never tries
/// the improper ones. in debug mode (full) nothing is reduced and every witness over the 24
/// proper matrices is counted.
///
/// the witness has det +1 when twin and -1 when only a reflected copy locks.
fn decide(p: &[Cell], m: &Mats, full: bool, ablate: u8, w: &mut Work) -> Decision {
    if ablate == 0 { return Decision::none(); }
    if ablate == 1 {
        let list = if full { &m.proper } else { &m.proper_red };
        for mm in list.iter() { w.build_only(p, mm); }
        return Decision::none();
    }
    if full {
        let (twin, wc, by_identity, first) = w.witnesses(p, &m.proper);
        if twin {
            return Decision {
                twin: true, mirror: true, w: wc, by_identity,
                witness: first.map(|(mi, v)| (m.proper[mi], v)),
            };
        }
        let imp = w.verdict(p, &m.improper);
        return Decision {
            twin: false, mirror: imp.is_some(), w: wc, by_identity,
            witness: imp.map(|(mi, v)| (m.improper[mi], v)),
        };
    }
    match w.verdict(p, &m.proper_red) {
        Some((mi, v)) => Decision {
            twin: true, mirror: true, w: 0, by_identity: mi == 0,
            witness: Some((m.proper_red[mi], v)),
        },
        None => {
            let imp = w.verdict(p, &m.improper_red);
            Decision {
                twin: false, mirror: imp.is_some(), w: 0, by_identity: false,
                witness: imp.map(|(mi, v)| (m.improper_red[mi], v)),
            }
        }
    }
}

fn record(ctx: &mut Ctx, m: &Mats, rf: &[Cell], aut: u32, achiral: bool) {
    let n = rf.len();
    ctx.seen += 1;
    if ctx.seen % (1 << 18) == 0 {
        eprintln!("... {} shapes recorded, level counts: {}",
            ctx.seen,
            (ctx.nmin..=ctx.nmax).map(|k| ctx.acc[k].all.to_string())
                .collect::<Vec<_>>().join(","));
        std::io::stderr().flush().ok();
    }
    if ctx.emit_at > 0 {
        // shardtable mode: just collect the level
        if n == ctx.emit_at {
            let flat: Vec<String> = rf.iter()
                .map(|c| format!("{} {} {}", c[0], c[1], c[2])).collect();
            ctx.emitted.push(format!("{} {}", n, flat.join(" ")));
        }
        return;
    }
    assert_eq!(48 % aut, 0, "|Aut| does not divide 48");
    let h1 = hash_ref(rf, SEED1);
    let h2 = hash_ref(rf, SEED2);
    let full = n <= ctx.witness_upto;
    let ablate = ctx.ablate;
    let d = decide(rf, m, full, ablate, &mut ctx.work);

    if ctx.samples.is_some() && (h2 & 0xFFFF) == 0 {
        // one shape in 65536, picked by its own second checksum and not by walk order, so both
        // programs dump the same set and a third can recheck it. every mirror-yes sample needs a
        // witness whose det matches the verdict (+1 by rotation, -1 only reflected).
        assert_eq!(d.mirror, d.witness.is_some(), "mirror verdict without a witness");
        if let Some((r, _)) = d.witness {
            assert_eq!(det(&r), if d.twin { 1 } else { -1 }, "witness determinant");
        }
        let cells: Vec<String> = rf.iter()
            .map(|c| format!("[{},{},{}]", c[0], c[1], c[2])).collect();
        let wit = match d.witness {
            Some((r, v)) => format!(
                ",\"witness\":{{\"R\":[[{},{},{}],[{},{},{}],[{},{},{}]],\"v\":[{},{},{}]}}",
                r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2],
                r[2][0], r[2][1], r[2][2], v[0], v[1], v[2]),
            None => String::new(),
        };
        let line = format!(
            "{{\"prog\":\"B\",\"n\":{},\"cells\":[{}],\"twin\":{},\"mirror\":{},\"h1\":\"{:016x}\",\"h2\":\"{:016x}\"{}}}",
            n, cells.join(","), d.twin, d.mirror, h1, h2, wit);
        let fh = ctx.samples.as_mut().unwrap();
        writeln!(fh, "{}", line).expect("write sample");
        ctx.acc[n].sampled += 1;
    }

    let a = &mut ctx.acc[n];
    a.all += 1;
    a.all_h.0 = a.all_h.0.wrapping_add(h1);
    a.all_h.1 = a.all_h.1.wrapping_add(h2);
    a.fixed_equiv += (48 / aut) as u64;
    if achiral { a.achiral += 1; }
    if full {
        a.wcount += d.w;
        a.wh1 = a.wh1.wrapping_add(d.w.wrapping_mul(h1));
    }
    if d.twin {
        a.twin += 1;
        a.twin_h.0 = a.twin_h.0.wrapping_add(h1);
        a.twin_h.1 = a.twin_h.1.wrapping_add(h2);
    }
    if d.mirror {
        a.mirror += 1;
        a.mirror_h.0 = a.mirror_h.0.wrapping_add(h1);
        a.mirror_h.1 = a.mirror_h.1.wrapping_add(h2);
    }
    if d.by_identity { a.translate += 1; }
}

/// canonical augmentation (McKay): a shape is only generated from its canonical parent, the
/// shape left after deleting one distinguished cell. so the free classes form a tree and this
/// walk hits each one exactly once, holding only the path from the single cell to the current
/// shape.
fn dfs(ctx: &mut Ctx, m: &Mats, rf: Vec<Cell>, aut: u32, achiral: bool) {
    record(ctx, m, &rf, aut, achiral);
    if rf.len() == ctx.nmax { return; }
    let mats = &m.all;
    let maps = stabiliser_maps(&rf, mats);
    // perimeter
    let mut perim: Vec<Cell> = Vec::new();
    for c in &rf {
        for d in NEIGH.iter() {
            let t = [c[0]+d[0], c[1]+d[1], c[2]+d[2]];
            if in_set(&rf, &t) { continue; }
            if let Err(pos) = perim.binary_search(&t) { perim.insert(pos, t); }
        }
    }
    for c in orbit_reps(&perim, &maps) {
        let mut raw: Vec<Cell> = rf.clone();
        raw.push(c);
        let (child, mn0) = normalize_with_shift(&raw);
        let cs = [c[0]-mn0[0], c[1]-mn0[1], c[2]-mn0[2]];
        let (crf, caut, cach, bi, bshift) = ref_full(&child, mats);
        let orbit = canonical_orbit(&crf, mats);
        // map the added cell into the REF frame with the matrix ref_full already found
        let x = apply(&mats[bi], &cs);
        let ic = [x[0]-bshift[0], x[1]-bshift[1], x[2]-bshift[2]];
        if orbit.binary_search(&ic).is_ok() { dfs(ctx, m, crf, caut, cach); }
    }
}

// reporting

#[allow(clippy::too_many_arguments)]
fn result_line(a: &Acc, n: usize, shard: &str, shards_total: usize,
               witness_upto: usize, build: &str, cpu: f64, wall: f64,
               ablate: u8) -> String {
    // no "leaves" here, this program never builds fixed polycubes (fixed_equiv does the A001931
    // check instead). "seconds" is cpu time since that adds up across shards, wall_seconds is
    // next to it.
    format!(
        "{{\"prog\":\"B\",\"n\":{},\"shard\":\"{}\",\"shards_total\":{},\
         \"all\":{{\"count\":{},\"h1\":\"{:016x}\",\"h2\":\"{:016x}\"}},\
         \"twin\":{{\"count\":{},\"h1\":\"{:016x}\",\"h2\":\"{:016x}\"}},\
         \"mirror\":{{\"count\":{},\"h1\":\"{:016x}\",\"h2\":\"{:016x}\"}},\
         \"fixed_equiv\":{},\"achiral\":{},\"translate\":{},\"sampled\":{},\
         \"wcount\":{},\"wh1\":\"{:016x}\",\"witness_upto\":{},\"ablate\":{},\
         \"build\":\"{}\",\"seconds\":{:.3},\"wall_seconds\":{:.3}}}\n",
        n, shard, shards_total,
        a.all, a.all_h.0, a.all_h.1,
        a.twin, a.twin_h.0, a.twin_h.1,
        a.mirror, a.mirror_h.0, a.mirror_h.1,
        a.fixed_equiv, a.achiral, a.translate, a.sampled,
        a.wcount, a.wh1, witness_upto, ablate, build, cpu, wall)
}

fn selftest_line(a: &Acc, n: usize) -> String {
    format!("{} {} {:016x} {:016x} {} {:016x} {:016x} {} {:016x} {:016x}",
        n, a.all, a.all_h.0, a.all_h.1, a.twin, a.twin_h.0, a.twin_h.1,
        a.mirror, a.mirror_h.0, a.mirror_h.1)
}

fn load_shard_table(path: &str) -> Vec<Vec<Cell>> {
    let text = std::fs::read_to_string(path).expect("read shard table");
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let nums: Vec<i32> = line.split_whitespace()
            .filter_map(|s| s.parse().ok()).collect();
        if nums.is_empty() { continue; }
        let n = nums[0] as usize;
        assert_eq!(nums.len(), 1 + 3 * n, "shard table line has the wrong field count");
        let mut cells: Vec<Cell> = Vec::with_capacity(n);
        for i in 0..n { cells.push([nums[1+3*i], nums[2+3*i], nums[3+3*i]]); }
        out.push(cells);
    }
    assert!(!out.is_empty(), "shard table is empty");
    let d = out[0].len();
    assert!(out.iter().all(|s| s.len() == d), "shard table mixes sizes");
    out
}

// main

fn main() {
    let args: Vec<String> = env::args().collect();
    let mode = args.get(1).cloned().unwrap_or_else(|| "selftest".into());
    let build = own_sha256();
    eprintln!("# twin_b5 build sha256 {}", build);
    std::io::stderr().flush().ok();

    let m = Mats::build();
    let ablate: u8 = env::var("TWINB_ABLATE").ok()
        .and_then(|v| v.parse().ok()).unwrap_or(2);
    if ablate != 2 {
        eprintln!("# ABLATION LEVEL {} -- timings only, counts are NOT valid", ablate);
    }

    if mode == "classify" {
        let path = args.get(2).expect("classify needs a file");
        let text = std::fs::read_to_string(path).expect("read");
        let mut w = Work::new();
        for line in text.lines() {
            let nums: Vec<i32> = line.split_whitespace()
                .filter_map(|s| s.parse().ok()).collect();
            if nums.is_empty() { continue; }
            let n = nums[0] as usize;
            assert!(n <= NMAX_SUPPORTED, "classify shape of size {} exceeds {}", n, NMAX_SUPPORTED);
            let mut cells: Vec<Cell> = Vec::with_capacity(n);
            for i in 0..n {
                cells.push([nums[1+3*i], nums[2+3*i], nums[3+3*i]]);
            }
            let (rf, aut, ach) = ref_aut_achiral(&cells, &m.all);
            // full: all 24 proper matrices, witnesses counted, no reduction. fast: the 17 reps with early
            // exit on both cosets. the labelled columns let the fast path's mirror verdict and identity
            // flag be checked per shape.
            let full = decide(&rf, &m, true, 2, &mut w);
            let fast = decide(&rf, &m, false, 2, &mut w);
            println!("{} {} {} {} {} {} {} {:016x} {:016x} mirror_fast={} translate_fast={}",
                full.twin as u8, fast.twin as u8, full.mirror as u8, full.w,
                full.by_identity as u8, aut, ach as u8,
                hash_ref(&rf, SEED1), hash_ref(&rf, SEED2),
                fast.mirror as u8, fast.by_identity as u8);
        }
        return;
    }

    if mode == "shardtable" {
        let d: usize = args.get(2).and_then(|s| s.parse().ok()).expect("shardtable needs d");
        let out = args.get(3).expect("shardtable needs an output file");
        assert!(d >= 1 && d <= NMAX_SUPPORTED, "d out of range");
        let mut ctx = Ctx {
            ablate, nmin: 1, nmax: d, witness_upto: 0,
            acc: (0..=d).map(|_| Acc::new()).collect(),
            seen: 0, work: Work::new(), samples: None,
            emit_at: d, emitted: Vec::new(),
        };
        let start = vec![[0, 0, 0]];
        let (rf, aut, ach) = ref_aut_achiral(&start, &m.all);
        dfs(&mut ctx, &m, rf, aut, ach);
        let mut text = String::new();
        text.push_str(&format!("# shard table, d = {}, {} shards, build {}\n",
                               d, ctx.emitted.len(), build));
        for line in ctx.emitted.iter() { text.push_str(line); text.push('\n'); }
        std::fs::write(out, text.as_bytes()).expect("write shard table");
        eprintln!("# {} shards of size {} written to {}", ctx.emitted.len(), d, out);
        return;
    }

    if mode == "shard" {
        // a shard is the subtree rooted at one free class of size d. every shape with n >= d has
        // exactly one ancestor of that size, so the shards are disjoint and cover everything >= d
        // (smaller shapes come from a separate run). write order is results, then samples, then the
        // done marker, so a half-written shard never gets collected and resuming is safe.
        let nmax: usize = args.get(2).and_then(|s| s.parse().ok()).expect("shard needs nmax");
        let wupto: usize = args.get(3).and_then(|s| s.parse().ok()).expect("shard needs wupto");
        let table = args.get(4).expect("shard needs a table path");
        let id: usize = args.get(5).and_then(|s| s.parse().ok()).expect("shard needs an id");
        let dir = args.get(6).expect("shard needs an output directory");
        assert!(nmax <= NMAX_SUPPORTED, "nmax {} exceeds {}", nmax, NMAX_SUPPORTED);

        // only level 2 is a real run, refuse anything else here
        assert_eq!(ablate, 2,
            "shard mode refuses TWINB_ABLATE={}; only 2 produces a reportable run", ablate);

        let shapes = load_shard_table(table);
        assert!(id < shapes.len(), "shard id {} is outside 0..{}", id, shapes.len() - 1);
        let d = shapes[0].len();
        assert!(d <= nmax, "shard depth {} exceeds nmax {}", d, nmax);
        std::fs::create_dir_all(dir).expect("create outdir");

        // a shard with a done marker is finished, skip it and never touch its files. one without a
        // marker gets rerun from scratch, which is safe since the driver only collects marked ones.
        let done = format!("{}/done_{}", dir, id);
        if std::path::Path::new(&done).exists() {
            eprintln!("# shard {} already done, nothing to do", id);
            return;
        }

        let rpath = format!("{}/B_results_{}.jsonl", dir, id);
        let spath = format!("{}/B_samples_{}.jsonl", dir, id);
        let fh = std::fs::File::create(&spath).expect("create samples file");
        let mut ctx = Ctx {
            ablate, nmin: d, nmax, witness_upto: wupto,
            acc: (0..=nmax).map(|_| Acc::new()).collect(),
            seen: 0, work: Work::new(),
            samples: Some(std::io::BufWriter::new(fh)),
            emit_at: 0, emitted: Vec::new(),
        };
        let (rf, aut, ach) = ref_aut_achiral(&shapes[id], &m.all);
        assert_eq!(rf, shapes[id], "shard table entry {} is not its own REF", id);
        let t0 = std::time::Instant::now();
        dfs(&mut ctx, &m, rf, aut, ach);
        let wall = t0.elapsed().as_secs_f64();
        let cpu = cputime::cpu_seconds();
        if let Some(mut s) = ctx.samples.take() { s.flush().expect("flush samples"); }

        let mut text = String::new();
        for n in d..=nmax {
            text.push_str(&result_line(&ctx.acc[n], n, &id.to_string(), shapes.len(),
                                       wupto, &build, cpu, wall, ablate));
        }
        std::fs::write(&rpath, text.as_bytes()).expect("write results");
        // marker goes last, so if it's there the results file is complete
        std::fs::write(&done, format!("{} {:.3}\n", id, cpu).as_bytes()).expect("write done");
        eprintln!("# shard {} of {} done, {} shapes, {:.3} CPU-s, {:.3} wall-s",
                  id, shapes.len(), ctx.seen, cpu, wall);
        return;
    }

    let nmax: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(9);
    let witness_upto: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(10);
    // dr.rs's stride only works up to NMAX_SUPPORTED, so check it once here
    assert!(nmax >= 1 && nmax <= NMAX_SUPPORTED,
            "nmax {} outside 1..{} (dr.rs packs D_R at stride 34)", nmax, NMAX_SUPPORTED);

    let (results_path, samples) = if mode == "report" {
        let dir = args.get(4).expect("report needs an output directory");
        assert_eq!(ablate, 2,
            "report mode refuses TWINB_ABLATE={}; only 2 produces a reportable run", ablate);
        std::fs::create_dir_all(dir).expect("create outdir");
        let rp = format!("{}/B_results.jsonl", dir);
        // don't overwrite results that are already there
        assert!(!std::path::Path::new(&rp).exists(),
                "{} already exists; move it aside or pick another directory", rp);
        let sp = format!("{}/B_samples.jsonl", dir);
        let fh = std::fs::File::create(&sp).expect("create samples file");
        (Some(rp), Some(std::io::BufWriter::new(fh)))
    } else {
        (None, None)
    };

    let mut ctx = Ctx {
        ablate, nmin: 1,
        nmax, witness_upto,
        acc: (0..=nmax).map(|_| Acc::new()).collect(),
        seen: 0,
        work: Work::new(),
        samples,
        emit_at: 0, emitted: Vec::new(),
    };
    let start = vec![[0, 0, 0]];
    let t0 = std::time::Instant::now();
    let (rf, aut, ach) = ref_aut_achiral(&start, &m.all);
    dfs(&mut ctx, &m, rf, aut, ach);
    let wall = t0.elapsed().as_secs_f64();
    let cpu = cputime::cpu_seconds();

    if let Some(mut fh) = ctx.samples.take() { fh.flush().expect("flush samples"); }

    let mut results = String::new();
    for n in 1..=nmax {
        let a = &ctx.acc[n];
        // the selftest line goes to stdout in both selftest and report mode, so one run gives the
        // diff against data/selftest_lines.txt and the input for a tally
        if mode == "count" {
            println!("n={} all={} twin={} mirror={} wcount={} wh1={:016x} fixed_equiv={} achiral={} translate={}",
                n, a.all, a.twin, a.mirror, a.wcount, a.wh1, a.fixed_equiv, a.achiral, a.translate);
        } else {
            println!("{}", selftest_line(a, n));
        }
        std::io::stdout().flush().ok();
        if results_path.is_some() {
            results.push_str(&result_line(a, n, "0", 1, witness_upto, &build, cpu, wall, ablate));
        }
    }

    if let Some(p) = results_path {
        std::fs::write(&p, results.as_bytes()).expect("write results");
        eprintln!("# wrote {} and the sample file, {:.3} CPU-s, {:.3} wall-s", p, cpu, wall);
    } else {
        eprintln!("# {:.3} CPU-s, {:.3} wall-s", cpu, wall);
    }
}
