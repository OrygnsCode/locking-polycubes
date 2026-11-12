// the locking test, done in offset space on the difference set D_R.
//
// D_R = { a - s : a in P, s in R(P) }, n^2 elements. for an offset v:
//   overlap(v)            <=>  v in D_R
//   face contact(v)       <=>  v + e in D_R for some signed unit e
//   blocked along d(v)    <=>  v + t*d in D_R for some t >= 1
// so each one is a bit test or a single comparison, no search over cells. why they hold:
//   - overlap: a cell c of R(P) lands on some a in P iff c + v = a iff v = a - c is in D_R.
//   - contact: a in P is next to s + v iff a - s - v = -e for a signed unit e, iff
//     v = (a - s) - e. NEIGH has both signs, so d - e over d in D_R and e in NEIGH gives
//     exactly the contact offsets.
//   - blocking: the copy slid by t along d hits P iff a = s + v + t*d for some a, s, iff
//     v + t*d is in D_R. along +x that's just: the max x in v's (y, z) lane of D_R is > v_x.
//
// the candidates are generated, not filtered: every offset with contact is d - e, so they're the
// 6n^2 neighbours of D_R, deduped with a bitset. tables get reset through lists of touched
// indices instead of being cleared (at most n^2 dirty entries against 39304 bits).

pub type Cell = [i32; 3];
pub type Mat = [[i32; 3]; 3];

pub const NEIGH: [Cell; 6] = [
    [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
];

/// stride per axis. D_R spans at most (box of P) + (box of R(P)) - 1 <= 33 for n <= 17, so 34
/// leaves one spare.
const STR: usize = 34;
const BITWORDS: usize = (STR * STR * STR + 63) / 64;
const LANE: usize = STR * STR;

#[inline]
fn widx(q0: usize, q1: usize, q2: usize) -> usize { (q0 * STR + q1) * STR + q2 }

pub struct Work {
    dr: Vec<Cell>,              // the n^2 differences, in q-space
    span: [usize; 3],
    lo: Cell,
    bits: Vec<u64>,             // membership of D_R
    bits_touched: Vec<usize>,
    vis: Vec<u64>,              // candidate dedupe
    vis_touched: Vec<usize>,
    lmin: Vec<i32>,             // 3 * LANE, lane minima in q-space
    lmax: Vec<i32>,
    lane_touched: Vec<usize>,   // indices into lmin/lmax
    s: Vec<Cell>,
}

impl Work {
    pub fn new() -> Work {
        Work {
            dr: Vec::with_capacity(512),
            span: [0; 3],
            lo: [0; 3],
            bits: vec![0u64; BITWORDS],
            bits_touched: Vec::with_capacity(512),
            vis: vec![0u64; BITWORDS],
            vis_touched: Vec::with_capacity(2048),
            lmin: vec![i32::MAX; 3 * LANE],
            lmax: vec![i32::MIN; 3 * LANE],
            lane_touched: Vec::with_capacity(2048),
            s: Vec::with_capacity(32),
        }
    }

    fn reset(&mut self) {
        for &w in self.bits_touched.iter() { self.bits[w] = 0; }
        self.bits_touched.clear();
        for &w in self.vis_touched.iter() { self.vis[w] = 0; }
        self.vis_touched.clear();
        for &i in self.lane_touched.iter() {
            self.lmin[i] = i32::MAX;
            self.lmax[i] = i32::MIN;
        }
        self.lane_touched.clear();
        self.dr.clear();
    }

    #[inline]
    fn bit_set(&mut self, idx: usize) {
        let w = idx >> 6;
        let b = 1u64 << (idx & 63);
        if self.bits[w] == 0 { self.bits_touched.push(w); }
        self.bits[w] |= b;
    }

    #[inline]
    fn bit_get(&self, idx: usize) -> bool { self.bits[idx >> 6] & (1u64 << (idx & 63)) != 0 }

    /// false if v was already seen
    #[inline]
    fn vis_first(&mut self, idx: usize) -> bool {
        let w = idx >> 6;
        let b = 1u64 << (idx & 63);
        if self.vis[w] & b != 0 { return false; }
        if self.vis[w] == 0 { self.vis_touched.push(w); }
        self.vis[w] |= b;
        true
    }

    /// builds D_R, its bounding box, the membership bitset and the three lane tables for one
    /// matrix. `p` has to be normalized (coords >= 0).
    fn build(&mut self, p: &[Cell], m: &Mat) {
        self.reset();
        self.s.clear();
        for c in p {
            self.s.push([
                m[0][0]*c[0] + m[0][1]*c[1] + m[0][2]*c[2],
                m[1][0]*c[0] + m[1][1]*c[1] + m[1][2]*c[2],
                m[2][0]*c[0] + m[2][1]*c[1] + m[2][2]*c[2],
            ]);
        }
        // raw differences, and their bounding box
        let mut lo = [i32::MAX; 3];
        let mut hi = [i32::MIN; 3];
        for a in p {
            for c in self.s.iter() {
                let d = [a[0]-c[0], a[1]-c[1], a[2]-c[2]];
                for k in 0..3 {
                    if d[k] < lo[k] { lo[k] = d[k]; }
                    if d[k] > hi[k] { hi[k] = d[k]; }
                }
                self.dr.push(d);
            }
        }
        for k in 0..3 {
            self.span[k] = (hi[k] - lo[k] + 1) as usize;
            assert!(self.span[k] <= STR, "D_R span {} exceeds stride {}", self.span[k], STR);
        }
        self.lo = lo;
        // shift into q-space and fill the tables
        for i in 0..self.dr.len() {
            let d = self.dr[i];
            let q = [d[0]-lo[0], d[1]-lo[1], d[2]-lo[2]];
            self.dr[i] = q;
            self.bit_set(widx(q[0] as usize, q[1] as usize, q[2] as usize));
            // axis 0 keyed on (q1, q2), axis 1 on (q0, q2), axis 2 on (q0, q1)
            let keys = [
                0 * LANE + (q[1] as usize) * STR + q[2] as usize,
                1 * LANE + (q[0] as usize) * STR + q[2] as usize,
                2 * LANE + (q[0] as usize) * STR + q[1] as usize,
            ];
            for k in 0..3 {
                let i2 = keys[k];
                if self.lmin[i2] == i32::MAX && self.lmax[i2] == i32::MIN {
                    self.lane_touched.push(i2);
                }
                if q[k] < self.lmin[i2] { self.lmin[i2] = q[k]; }
                if q[k] > self.lmax[i2] { self.lmax[i2] = q[k]; }
            }
        }
    }

    /// blocked along all six axis directions? `q` is in q-space and inside the box on every axis.
    /// an empty lane can't block (lmax stays i32::MIN, lmin i32::MAX).
    #[inline]
    fn locked(&self, q: [i32; 3]) -> bool {
        let k0 = 0 * LANE + (q[1] as usize) * STR + q[2] as usize;
        if !(self.lmax[k0] > q[0] && self.lmin[k0] < q[0]) { return false; }
        let k1 = 1 * LANE + (q[0] as usize) * STR + q[2] as usize;
        if !(self.lmax[k1] > q[1] && self.lmin[k1] < q[1]) { return false; }
        let k2 = 2 * LANE + (q[0] as usize) * STR + q[1] as usize;
        self.lmax[k2] > q[2] && self.lmin[k2] < q[2]
    }

    /// scans one matrix, calling `hit(v)` for every locking offset v (original frame). returns the
    /// number of hits, stops at the first one if `early`.
    fn scan_matrix<F: FnMut(Cell)>(&mut self, p: &[Cell], m: &Mat, early: bool, mut hit: F) -> u64 {
        self.build(p, m);
        let (sx, sy, sz) = (self.span[0] as i32, self.span[1] as i32, self.span[2] as i32);
        let mut found = 0u64;
        // candidates are the signed-unit neighbours of D_R. a v outside D_R's box on some axis has an
        // empty lane there, so it can't be blocked both ways on that axis. skipping it is the same
        // test done early, not a prefilter.
        let ndr = self.dr.len();
        for i in 0..ndr {
            let d = self.dr[i];
            for e in NEIGH.iter() {
                let q = [d[0]-e[0], d[1]-e[1], d[2]-e[2]];
                if q[0] <= 0 || q[0] >= sx - 1 { continue; }
                if q[1] <= 0 || q[1] >= sy - 1 { continue; }
                if q[2] <= 0 || q[2] >= sz - 1 { continue; }
                let idx = widx(q[0] as usize, q[1] as usize, q[2] as usize);
                if self.bit_get(idx) { continue; }        // overlap
                if !self.vis_first(idx) { continue; }     // already tried
                if self.locked(q) {
                    found += 1;
                    hit([q[0] + self.lo[0], q[1] + self.lo[1], q[2] + self.lo[2]]);
                    if early { return found; }
                }
            }
        }
        found
    }

    /// verdict with early exit: (index into `mats`, v) of the first witness, original frame
    pub fn verdict(&mut self, p: &[Cell], mats: &[Mat]) -> Option<(usize, Cell)> {
        for (mi, m) in mats.iter().enumerate() {
            let mut got: Option<Cell> = None;
            self.scan_matrix(p, m, true, |v| got = Some(v));
            if let Some(v) = got { return Some((mi, v)); }
        }
        None
    }

    /// every witness over every matrix, no early exit, no reduction. returns (locks, count,
    /// identity was a witness, first witness).
    pub fn witnesses(&mut self, p: &[Cell], mats: &[Mat])
        -> (bool, u64, bool, Option<(usize, Cell)>)
    {
        let mut w = 0u64;
        let mut by_identity = false;
        let mut first: Option<(usize, Cell)> = None;
        for (mi, m) in mats.iter().enumerate() {
            let mut one: Option<Cell> = None;
            let got = self.scan_matrix(p, m, false, |v| { if one.is_none() { one = Some(v); } });
            if got > 0 {
                if mi == 0 { by_identity = true; }
                if first.is_none() { first = Some((mi, one.unwrap())); }
                w += got;
            }
        }
        (w > 0, w, by_identity, first)
    }

    /// builds D_R only, for the ablation timing
    pub fn build_only(&mut self, p: &[Cell], m: &Mat) { self.build(p, m); }
}
