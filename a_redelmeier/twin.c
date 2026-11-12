/* twin.c - counts free polycubes that lock with a rotated copy of themselves.
   see the README for the definitions and the checksums.

   shapes come from Redelmeier's algorithm over fixed polycubes. a fixed one is only kept if it's
   the canonical image of its free class: of its 48 images with a sorted box (b0 >= b1 >= b2), the
   one with the smallest sorted cell list.

   the test works on D = { a - M(c) : a, c in P }. the copy M(P) + v overlaps P iff v is in D,
   touches it iff v + e is in D for a unit e, and is blocked along d iff v + t*d is in D for some
   t >= 1, which is just a min/max check on v's line of D.

   build: gcc -O3 -march=native -o twin.exe twin.c
   usage:
     twin count N [debug]                 one n in one process, prints one JSON line
     twin selftest REFFILE                checks n = 1.. against the reference file
     twin mkshards N D SHARDFILE          writes the shard table (search states at depth D)
     twin run N SHARDFILE DONEFILE RESULTS SAMPLES WORKER W REFFILE
*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

/* cpu time (user + kernel), wall time depends on whatever else is running */
static double cpu_seconds(void) {
    FILETIME a, b, k, u;
    GetProcessTimes(GetCurrentProcess(), &a, &b, &k, &u);
    ULARGE_INTEGER kk, uu;
    kk.LowPart = k.dwLowDateTime; kk.HighPart = k.dwHighDateTime;
    uu.LowPart = u.dwLowDateTime; uu.HighPart = u.dwHighDateTime;
    return (double)(kk.QuadPart + uu.QuadPart) * 1e-7;
}

typedef uint64_t u64;
typedef uint8_t u8;

#define MAXN 20
#define KEY(x, y, z) (((x) << 10) | ((y) << 5) | (z))

/* the 48 signed permutation matrices, (M x)_i = s[i] * x[p[i]] */

typedef struct { int p[3], s[3], det; } Mat;
static Mat M48[48];
static int inv48[48], ident;
static int prop24[24], improp24[24];   /* all proper / all improper */
static int prop_red[24], nprop_red;    /* proper, one from each {M, M^-1} */
static int improp_red[24], nimprop_red;
static int cand[4][48], ncand[4];      /* transforms keeping a sorted box sorted, by tie pattern */

static void init_mats(void) {
    static const int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    int k = 0;
    for (int a = 0; a < 6; a++) {
        int inv = 0;
        for (int i = 0; i < 3; i++)
            for (int j = i + 1; j < 3; j++)
                if (perms[a][i] > perms[a][j]) inv++;
        for (int sg = 0; sg < 8; sg++) {
            Mat *m = &M48[k++];
            for (int i = 0; i < 3; i++) { m->p[i] = perms[a][i]; m->s[i] = (sg >> i & 1) ? -1 : 1; }
            m->det = ((inv & 1) ? -1 : 1) * m->s[0] * m->s[1] * m->s[2];
        }
    }
    int np = 0, ni = 0;
    for (int a = 0; a < 48; a++) {
        Mat *m = &M48[a];
        if (m->p[0] == 0 && m->p[1] == 1 && m->p[2] == 2 && m->s[0] == 1 && m->s[1] == 1 && m->s[2] == 1)
            ident = a;
        inv48[a] = -1;
        for (int b = 0; b < 48; b++) {  /* (M M') x = x  iff  p'[p[i]] = i and s[i] s'[p[i]] = 1 */
            Mat *n = &M48[b];
            int ok = 1;
            for (int i = 0; i < 3; i++)
                if (n->p[m->p[i]] != i || m->s[i] * n->s[m->p[i]] != 1) ok = 0;
            if (ok) inv48[a] = b;
        }
        if (inv48[a] < 0) { fprintf(stderr, "no inverse for matrix %d\n", a); exit(2); }
        if (m->det > 0) prop24[np++] = a; else improp24[ni++] = a;
    }
    if (np != 24 || ni != 24) { fprintf(stderr, "matrix count wrong\n"); exit(2); }
    nprop_red = nimprop_red = 0;
    for (int a = 0; a < 48; a++) {
        if (inv48[a] < a) continue;   /* keep a if it's its own inverse or the smaller of the two */
        if (M48[a].det > 0) prop_red[nprop_red++] = a; else improp_red[nimprop_red++] = a;
    }
    if (nprop_red != 17 || nimprop_red != 17) {
        fprintf(stderr, "expected 17 + 17 reduced matrices, got %d + %d\n", nprop_red, nimprop_red);
        exit(2);
    }
    static const int pat[4][3] = {{3,2,1},{2,2,1},{2,1,1},{1,1,1}};
    for (int t = 0; t < 4; t++) {
        ncand[t] = 0;
        for (int a = 0; a < 48; a++) {
            const int *b = pat[t], *p = M48[a].p;
            if (b[p[0]] == b[0] && b[p[1]] == b[1] && b[p[2]] == b[2] && a != ident)
                cand[t][ncand[t]++] = a;
        }
    }
}

/* helpers */

static void isort(int *a, int n) {
    for (int i = 1; i < n; i++) {
        int v = a[i], j = i - 1;
        while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; j--; }
        a[j + 1] = v;
    }
}

static int lexcmp(const int *a, const int *b, int n) {
    for (int i = 0; i < n; i++)
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* keys of M(P) shifted into its box (c[] is P in box b[]), returns the smallest key */
static int image_keys(const Mat *m, int (*c)[3], const int *b, int n, int *out) {
    int mn = 1 << 30;
    for (int j = 0; j < n; j++) {
        int q[3];
        for (int i = 0; i < 3; i++) {
            int v = c[j][m->p[i]];
            q[i] = m->s[i] > 0 ? v : b[m->p[i]] - 1 - v;
        }
        int k = KEY(q[0], q[1], q[2]);
        out[j] = k;
        if (k < mn) mn = k;
    }
    return mn;
}

static u64 splitmix64(u64 x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static u64 hash_ref(const int *ref, int n, u64 seed) {
    u64 h = seed;
    for (int j = 0; j < n; j++) {
        u64 x = (u64)(ref[j] >> 10), y = (u64)((ref[j] >> 5) & 31), z = (u64)(ref[j] & 31);
        h = splitmix64(h ^ (x | (y << 8) | (z << 16)));
    }
    return splitmix64(h ^ (u64)n);
}

#define SEED1 0x243F6A8885A308D3ULL
#define SEED2 0x13198A2E03707344ULL

/* locking test */

#define GMAX (2 * MAXN + 2)
static u8 grid[GMAX * GMAX * GMAX];
static u8 lxmin[GMAX * GMAX], lxmax[GMAX * GMAX], lymin[GMAX * GMAX], lymax[GMAX * GMAX],
          lzmin[GMAX * GMAX], lzmax[GMAX * GMAX];
static int dpts[MAXN * MAXN][3];

typedef struct {
    int found;          /* some witness exists */
    u64 witnesses;      /* number of v (count_all only) */
    int nocontact;      /* some v is free and blocked but not touching (count_all only) */
    int wm, wv[3];      /* first witness: matrix index and v */
} TestOut;

/* test one matrix; count_all = 0 stops at the first witness */
static void test_matrix(int mi, int (*c)[3], const int *b, int n, int count_all, TestOut *o) {
    const Mat *m = &M48[mi];
    int lo[3], G[3];
    for (int i = 0; i < 3; i++) {
        int bp = b[m->p[i]];
        lo[i] = m->s[i] > 0 ? -(bp - 1) : 0;
        G[i] = b[i] + bp - 1 + 2;
    }
    int G01 = G[0] * G[1];
    memset(grid, 0, (size_t)G01 * G[2]);
    memset(lxmin, 255, G[1] * G[2]); memset(lxmax, 0, G[1] * G[2]);
    memset(lymin, 255, G[0] * G[2]); memset(lymax, 0, G[0] * G[2]);
    memset(lzmin, 255, G01);         memset(lzmax, 0, G01);
    int rc[MAXN][3];
    for (int j = 0; j < n; j++)
        for (int i = 0; i < 3; i++) rc[j][i] = m->s[i] * c[j][m->p[i]];
    int nd = 0;
    for (int a = 0; a < n; a++)
        for (int j = 0; j < n; j++) {
            int g0 = c[a][0] - rc[j][0] - lo[0] + 1;
            int g1 = c[a][1] - rc[j][1] - lo[1] + 1;
            int g2 = c[a][2] - rc[j][2] - lo[2] + 1;
            int idx = g0 + G[0] * g1 + G01 * g2;
            if (grid[idx]) continue;
            grid[idx] = 1;
            dpts[nd][0] = g0; dpts[nd][1] = g1; dpts[nd][2] = g2; nd++;
            int lx = g1 + G[1] * g2, ly = g0 + G[0] * g2, lz = g0 + G[0] * g1;
            if (g0 < lxmin[lx]) lxmin[lx] = (u8)g0;
            if (g0 > lxmax[lx]) lxmax[lx] = (u8)g0;
            if (g1 < lymin[ly]) lymin[ly] = (u8)g1;
            if (g1 > lymax[ly]) lymax[ly] = (u8)g1;
            if (g2 < lzmin[lz]) lzmin[lz] = (u8)g2;
            if (g2 > lzmax[lz]) lzmax[lz] = (u8)g2;
        }
    static const int E[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    for (int k = 0; k < nd; k++)
        for (int e = 0; e < 6; e++) {
            int v0 = dpts[k][0] - E[e][0], v1 = dpts[k][1] - E[e][1], v2 = dpts[k][2] - E[e][2];
            int idx = v0 + G[0] * v1 + G01 * v2;
            if (grid[idx]) continue;   /* in D (overlap) or already tested */
            grid[idx] = 2;
            int lx = v1 + G[1] * v2, ly = v0 + G[0] * v2, lz = v0 + G[0] * v1;
            if (lxmin[lx] < v0 && lxmax[lx] > v0 && lymin[ly] < v1 && lymax[ly] > v1 &&
                lzmin[lz] < v2 && lzmax[lz] > v2) {
                if (!o->found) {
                    o->found = 1; o->wm = mi;
                    o->wv[0] = v0 - 1 + lo[0]; o->wv[1] = v1 - 1 + lo[1]; o->wv[2] = v2 - 1 + lo[2];
                }
                if (!count_all) return;
                o->witnesses++;
            }
        }
    if (count_all && !o->nocontact) {
        for (int v2 = 0; v2 < G[2]; v2++)
            for (int v1 = 0; v1 < G[1]; v1++)
                for (int v0 = 0; v0 < G[0]; v0++) {
                    if (grid[v0 + G[0] * v1 + G01 * v2] == 1) continue;
                    int lx = v1 + G[1] * v2, ly = v0 + G[0] * v2, lz = v0 + G[0] * v1;
                    if (lxmin[lx] < v0 && lxmax[lx] > v0 && lymin[ly] < v1 && lymax[ly] > v1 &&
                        lzmin[lz] < v2 && lzmax[lz] > v2) { o->nocontact = 1; return; }
                }
    }
}

/* same test on bitmasks. D is one 64-bit word per (y, z) line, bit g0 = x. blocked along +-x
   means strictly between the lowest and highest bit of the line, along +-y or +-z it's an OR of the
   lines below and above in that column. a witness is a bit blocked all six ways, not in D, and
   next to a bit of D. */
typedef uint64_t row_t;
static int prow_y[MAXN], prow_z[MAXN], nprow;
static row_t prow_m[MAXN];
static const void *prow_owner;   /* the coordinate array the rows were built from */
static row_t Dr[GMAX * GMAX], aboveY[GMAX * GMAX], aboveZ[GMAX * GMAX], belowZ[GMAX];

static void set_shape_rows(int (*c)[3], int n) {
    prow_owner = c;
    nprow = 0;
    for (int j = 0; j < n; j++) {
        int k;
        for (k = 0; k < nprow; k++)
            if (prow_y[k] == c[j][1] && prow_z[k] == c[j][2]) break;
        if (k == nprow) { prow_y[k] = c[j][1]; prow_z[k] = c[j][2]; prow_m[k] = 0; nprow++; }
        prow_m[k] |= (row_t)1 << c[j][0];
    }
}

static inline row_t between(row_t d) {
    if (!d) return 0;
    int l = __builtin_ctzll(d), h = 63 - __builtin_clzll(d);
    if (h - l < 2) return 0;
    return (((row_t)1 << h) - 1) & ~(((row_t)2 << l) - 1);
}

static void test_matrix_bits(int mi, int (*c)[3], const int *b, int n, int count_all, TestOut *o) {
    if (prow_owner != (const void *)c) { fprintf(stderr, "FATAL: line masks belong to another shape\n"); exit(3); }
    const Mat *m = &M48[mi];
    int lo[3], G[3];
    for (int i = 0; i < 3; i++) {
        int bp = b[m->p[i]];
        lo[i] = m->s[i] > 0 ? -(bp - 1) : 0;
        G[i] = b[i] + bp - 1 + 2;
    }
    int G1 = G[1], G2 = G[2], NR = G1 * G2;
    memset(Dr, 0, NR * sizeof(row_t));
    for (int j = 0; j < n; j++) {
        int r0 = m->s[0] * c[j][m->p[0]], r1 = m->s[1] * c[j][m->p[1]], r2 = m->s[2] * c[j][m->p[2]];
        int s = 1 - lo[0] - r0, oy = 1 - lo[1] - r1, oz = 1 - lo[2] - r2;
        for (int k = 0; k < nprow; k++)
            Dr[(prow_y[k] + oy) + G1 * (prow_z[k] + oz)] |= prow_m[k] << s;
    }
    for (int gz = 0; gz < G2; gz++) {           /* OR of the lines above, along y */
        row_t acc = 0;
        for (int gy = G1 - 1; gy >= 0; gy--) { aboveY[gy + G1 * gz] = acc; acc |= Dr[gy + G1 * gz]; }
    }
    for (int gy = 0; gy < G1; gy++) {           /* OR of the lines above, along z */
        row_t acc = 0;
        for (int gz = G2 - 1; gz >= 0; gz--) { aboveZ[gy + G1 * gz] = acc; acc |= Dr[gy + G1 * gz]; }
        belowZ[gy] = 0;
    }
    for (int gz = 0; gz < G2; gz++) {
        row_t belowY = 0;
        for (int gy = 0; gy < G1; gy++) {
            int i = gy + G1 * gz;
            row_t d = Dr[i];
            row_t blocked = between(d) & belowY & aboveY[i] & belowZ[gy] & aboveZ[i] & ~d;
            belowY |= d; belowZ[gy] |= d;
            if (!blocked) continue;
            if (count_all) o->nocontact = 1;
            /* blocked bits are interior, so the four neighbour lines exist */
            row_t w = blocked & ((d << 1) | (d >> 1) | Dr[i - 1] | Dr[i + 1] | Dr[i - G1] | Dr[i + G1]);
            if (!w) continue;
            if (!o->found) {
                o->found = 1; o->wm = mi;
                o->wv[0] = __builtin_ctzll(w) - 1 + lo[0]; o->wv[1] = gy - 1 + lo[1]; o->wv[2] = gz - 1 + lo[2];
            }
            if (!count_all) return;
            o->witnesses += (u64)__builtin_popcountll(w);
        }
    }
}

static int locks(const int *mats, int nm, int (*c)[3], const int *b, int n, int count_all, TestOut *o) {
    for (int t = 0; t < nm; t++) {
        test_matrix_bits(mats[t], c, b, n, count_all, o);
        if (o->found && !count_all) return 1;
    }
    return o->found;
}

/* counters */

typedef struct { u64 c, h1, h2; } Acc;
typedef struct {
    Acc all, twin, mirror;
    u64 fixed_equiv, achiral, leaves, translate, sampled, wcount, wh1, nocontact;
} Counts;

static int N, DEBUG, NOTEST;
static Counts K;
static FILE *samples;

static void acc_add(Acc *a, u64 h1, u64 h2) { a->c++; a->h1 += h1; a->h2 += h2; }

static void dump_sample(const int *ref, int n, int twin, int mirror) {
    int c[MAXN][3], b[3] = {0, 0, 0};
    for (int j = 0; j < n; j++) {
        c[j][0] = ref[j] >> 10; c[j][1] = (ref[j] >> 5) & 31; c[j][2] = ref[j] & 31;
        for (int i = 0; i < 3; i++) if (c[j][i] + 1 > b[i]) b[i] = c[j][i] + 1;
    }
    set_shape_rows(c, n);
    TestOut o; memset(&o, 0, sizeof o);
    int t2 = locks(prop24, 24, c, b, n, 0, &o);
    int m2 = t2 || locks(improp24, 24, c, b, n, 0, &o);
    if (t2 != twin || m2 != mirror) {
        fprintf(stderr, "FATAL: verdict on REF coordinates differs from the enumerated copy\n");
        exit(3);
    }
    if (!samples) return;
    K.sampled++;
    fprintf(samples, "{\"n\": %d, \"ref\": [", n);
    for (int j = 0; j < n; j++) fprintf(samples, "%s[%d, %d, %d]", j ? ", " : "", c[j][0], c[j][1], c[j][2]);
    fprintf(samples, "], \"twin\": %d, \"mirror\": %d, \"witness\": ", twin, mirror);
    if (o.found) {
        const Mat *m = &M48[o.wm];
        int R[3][3] = {{0}};
        for (int i = 0; i < 3; i++) R[i][m->p[i]] = m->s[i];
        fprintf(samples, "{\"R\": [[%d, %d, %d], [%d, %d, %d], [%d, %d, %d]], \"v\": [%d, %d, %d]}",
                R[0][0], R[0][1], R[0][2], R[1][0], R[1][1], R[1][2], R[2][0], R[2][1], R[2][2],
                o.wv[0], o.wv[1], o.wv[2]);
    } else {
        fprintf(samples, "null");
    }
    fprintf(samples, "}\n");
}

/* P in normalized coords c, sorted keys P, box b (sorted); |Aut| and chirality already known */
static void emit(int (*c)[3], const int *b, int n, int aut, int achiral) {
    int ref[MAXN], img[MAXN], have = 0;
    for (int a = 0; a < 48; a++) {
        int mn = image_keys(&M48[a], c, b, n, img);
        if (have && mn > ref[0]) continue;
        isort(img, n);
        if (!have || lexcmp(img, ref, n) < 0) { memcpy(ref, img, n * sizeof(int)); have = 1; }
    }
    if (48 % aut) { fprintf(stderr, "FATAL: |Aut| = %d does not divide 48\n", aut); exit(3); }
    u64 h1 = hash_ref(ref, n, SEED1), h2 = hash_ref(ref, n, SEED2);
    acc_add(&K.all, h1, h2);
    K.fixed_equiv += (u64)(48 / aut);
    if (achiral) K.achiral++;

    if (NOTEST) return;   /* profiling: everything except the locking test */
    /* identity first, gives the translate-only count for free */
    set_shape_rows(c, n);
    TestOut o; memset(&o, 0, sizeof o);
    test_matrix_bits(ident, c, b, n, DEBUG, &o);
    if (o.found) K.translate++;
    if (DEBUG) {
        for (int t = 0; t < 24; t++)
            if (prop24[t] != ident) test_matrix_bits(prop24[t], c, b, n, 1, &o);
        K.wcount += o.witnesses;
        K.wh1 += o.witnesses * h1;
        if (o.nocontact) K.nocontact++;
    } else {
        for (int t = 0; t < nprop_red && !o.found; t++)
            if (prop_red[t] != ident) test_matrix_bits(prop_red[t], c, b, n, 0, &o);
    }
    int twin = o.found;
    TestOut o2; memset(&o2, 0, sizeof o2);
    int mirror = twin || (DEBUG ? locks(improp24, 24, c, b, n, 0, &o2)
                                : locks(improp_red, nimprop_red, c, b, n, 0, &o2));
    if (twin) acc_add(&K.twin, h1, h2);
    if (mirror) acc_add(&K.mirror, h1, h2);
    if ((h2 & 0xFFFF) == 0) dump_sample(ref, n, twin, mirror);
}

/* Redelmeier enumeration */

static int GX, GY, GXY;
static u8 *seen, *seen_base;
static int8_t *gx, *gy, *gz;
static int origin;
static int off[6];
static int cell[MAXN][3];
static int mn_[MAXN + 1][3], mx_[MAXN + 1][3];

static void init_grid(int n) {
    GX = 2 * n + 1; GY = 2 * n + 1; GXY = GX * GY;
    int GZ = n + 2, total = GXY * GZ;
    free(seen); free(seen_base); free(gx); free(gy); free(gz);
    seen = malloc(total); seen_base = malloc(total);
    gx = malloc(total); gy = malloc(total); gz = malloc(total);
    for (int z = -1; z <= n; z++)
        for (int y = -n; y <= n; y++)
            for (int x = -n; x <= n; x++) {
                int idx = (x + n) + GX * (y + n) + GXY * (z + 1);
                gx[idx] = (int8_t)x; gy[idx] = (int8_t)y; gz[idx] = (int8_t)z;
                int inside = x >= -(n - 1) && x <= n - 1 && y >= -(n - 1) && y <= n - 1 && z >= 0 && z <= n - 1;
                int legal = inside && (z > 0 || (z == 0 && y > 0) || (z == 0 && y == 0 && x >= 0));
                seen_base[idx] = (u8)!legal;
            }
    origin = n + GX * n + GXY * 1;
    seen_base[origin] = 1;
    memcpy(seen, seen_base, total);
    off[0] = 1; off[1] = -1; off[2] = GX; off[3] = -GX; off[4] = GXY; off[5] = -GXY;
}

static void leaf(void) {
    K.leaves++;
    const int *lo = mn_[N], *hi = mx_[N];
    int b[3] = {hi[0] - lo[0] + 1, hi[1] - lo[1] + 1, hi[2] - lo[2] + 1};
    if (b[0] < b[1] || b[1] < b[2]) return;
    int c[MAXN][3], P[MAXN], img[MAXN];
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < 3; i++) c[j][i] = cell[j][i] - lo[i];
        P[j] = KEY(c[j][0], c[j][1], c[j][2]);
    }
    isort(P, N);
    int pat = b[0] == b[1] ? (b[1] == b[2] ? 3 : 1) : (b[1] == b[2] ? 2 : 0);
    int aut = 1, achiral = 0;
    for (int t = 0; t < ncand[pat]; t++) {
        const Mat *m = &M48[cand[pat][t]];
        int mn = image_keys(m, c, b, N, img);
        if (mn < P[0]) return;
        if (mn > P[0]) continue;
        isort(img, N);
        int r = lexcmp(img, P, N);
        if (r < 0) return;
        if (r == 0) { aut++; if (m->det < 0) achiral = 1; }
    }
    emit(c, b, N, aut, achiral);
}

static void set_cell(int depth, int idx) {
    cell[depth][0] = gx[idx]; cell[depth][1] = gy[idx]; cell[depth][2] = gz[idx];
    for (int i = 0; i < 3; i++) {
        int v = cell[depth][i];
        mn_[depth + 1][i] = depth == 0 || v < mn_[depth][i] ? v : mn_[depth][i];
        mx_[depth + 1][i] = depth == 0 || v > mx_[depth][i] ? v : mx_[depth][i];
    }
}

static void rec(const int *U, int m, int depth) {
    int newU[6 * MAXN + 8], added[6];
    while (m > 0) {
        int cidx = U[--m];
        set_cell(depth, cidx);
        if (depth + 1 == N) { leaf(); continue; }
        memcpy(newU, U, m * sizeof(int));
        int k = m, na = 0;
        for (int e = 0; e < 6; e++) {
            int d = cidx + off[e];
            if (!seen[d]) { seen[d] = 1; newU[k++] = d; added[na++] = d; }
        }
        rec(newU, k, depth + 1);
        for (int i = 0; i < na; i++) seen[added[i]] = 0;
    }
}

/* shard generation: record every search state at depth D */
static FILE *shard_out;
static int SD;
static u64 nshards;

static void gen(const int *U, int m, int depth) {
    int newU[6 * MAXN + 8], added[6];
    while (m > 0) {
        int cidx = U[--m];
        set_cell(depth, cidx);
        memcpy(newU, U, m * sizeof(int));
        int k = m, na = 0;
        for (int e = 0; e < 6; e++) {
            int d = cidx + off[e];
            if (!seen[d]) { seen[d] = 1; newU[k++] = d; added[na++] = d; }
        }
        if (depth + 1 == SD) {
            fprintf(shard_out, "%llu %d", (unsigned long long)nshards++, SD);
            for (int j = 0; j < SD; j++) fprintf(shard_out, " %d %d %d", cell[j][0], cell[j][1], cell[j][2]);
            fprintf(shard_out, " %d", k);
            for (int j = 0; j < k; j++) fprintf(shard_out, " %d %d %d", gx[newU[j]], gy[newU[j]], gz[newU[j]]);
            fprintf(shard_out, "\n");
        } else {
            gen(newU, k, depth + 1);
        }
        for (int i = 0; i < na; i++) seen[added[i]] = 0;
    }
}

static int idx_of(int x, int y, int z) { return (x + N) + GX * (y + N) + GXY * (z + 1); }

/* modes */

static void count_all_n(int n, int debug) {
    N = n; DEBUG = debug;
    memset(&K, 0, sizeof K);
    init_grid(n);
    int U[1] = {origin};
    rec(U, 1, 0);
}

static void print_counts(FILE *f, const char *prog, int n, const char *shard, u64 total, double secs) {
    fprintf(f, "{\"prog\": \"%s\", \"n\": %d, \"shard\": \"%s\", \"shards_total\": %llu, ", prog, n, shard,
            (unsigned long long)total);
    const Acc *a[3] = {&K.all, &K.twin, &K.mirror};
    const char *nm[3] = {"all", "twin", "mirror"};
    for (int i = 0; i < 3; i++)
        fprintf(f, "\"%s\": {\"count\": %llu, \"h1\": \"%016llx\", \"h2\": \"%016llx\"}, ", nm[i],
                (unsigned long long)a[i]->c, (unsigned long long)a[i]->h1, (unsigned long long)a[i]->h2);
    fprintf(f, "\"fixed_equiv\": %llu, \"achiral\": %llu, \"leaves\": %llu, \"translate\": %llu, \"sampled\": %llu, ",
            (unsigned long long)K.fixed_equiv, (unsigned long long)K.achiral, (unsigned long long)K.leaves,
            (unsigned long long)K.translate, (unsigned long long)K.sampled);
    if (DEBUG)
        fprintf(f, "\"wcount\": %llu, \"wh1\": \"%016llx\", \"nocontact\": %llu, ", (unsigned long long)K.wcount,
                (unsigned long long)K.wh1, (unsigned long long)K.nocontact);
    fprintf(f, "\"seconds\": %.3f}\n", secs);
    fflush(f);
}

/* known values from OEIS A038119, A001931, A007743 */
static const u64 A038119[16] = {0, 1, 1, 2, 7, 23, 112, 607, 3811, 25413, 178083, 1279537, 9371094,
                                69513546, 520878101, 3934285874ULL};
static const u64 A001931[16] = {0, 1, 3, 15, 86, 534, 3481, 23502, 162913, 1152870, 8294738, 60494549,
                                446205905, 3322769321ULL, 24946773111ULL, 188625900446ULL};
static const u64 A007743[16] = {0, 1, 1, 2, 6, 17, 58, 191, 700, 2515, 9623, 36552, 143761, 564443,
                                2259905, 9057278};

/* reference file: lines "n all_c all_h1 all_h2 twin_c twin_h1 twin_h2 mirror_c mirror_h1 mirror_h2" */
static int selftest(const char *reffile) {
    FILE *f = fopen(reffile, "r");
    if (!f) { fprintf(stderr, "selftest: cannot open %s\n", reffile); return 0; }
    int ok = 1, lines = 0, n;
    unsigned long long v[9];
    FILE *save = samples; samples = NULL;
    while (fscanf(f, "%d %llu %llx %llx %llu %llx %llx %llu %llx %llx", &n, &v[0], &v[1], &v[2], &v[3], &v[4],
                  &v[5], &v[6], &v[7], &v[8]) == 10) {
        lines++;
        count_all_n(n, 0);
        const Acc *a[3] = {&K.all, &K.twin, &K.mirror};
        for (int i = 0; i < 3; i++)
            if (a[i]->c != v[3 * i] || a[i]->h1 != v[3 * i + 1] || a[i]->h2 != v[3 * i + 2]) ok = 0;
        if (n < 16 && (K.all.c != A038119[n] || K.fixed_equiv != A001931[n] || K.achiral != A007743[n] ||
                       K.leaves != A001931[n])) ok = 0;
        fprintf(stderr, "selftest n=%d %s\n", n, ok ? "ok" : "MISMATCH");
        if (!ok) break;
    }
    fclose(f);
    samples = save;
    if (lines < 9) { fprintf(stderr, "selftest: only %d reference lines\n", lines); ok = 0; }
    return ok;
}

int main(int argc, char **argv) {
    init_mats();
    if (argc >= 3 && !strcmp(argv[1], "count")) {
        int n = atoi(argv[2]);
        if (n < 1 || n > MAXN) return 2;
        double t0 = cpu_seconds();
        NOTEST = argc >= 4 && !strcmp(argv[3], "notest");
        count_all_n(n, argc >= 4 && !strcmp(argv[3], "debug"));
        print_counts(stdout, "A", n, "all", 1, cpu_seconds() - t0);
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "classify")) {
        /* input lines: n x1 y1 z1 ... ; output: twin twin_fast mirror w24 translate aut achiral h1 h2 */
        FILE *f = fopen(argv[2], "r");
        if (!f) return 2;
        int n;
        while (fscanf(f, "%d", &n) == 1) {
            int raw[MAXN][3], lo[3] = {99, 99, 99}, hi[3] = {-99, -99, -99};
            for (int j = 0; j < n; j++)
                for (int i = 0; i < 3; i++) {
                    if (fscanf(f, "%d", &raw[j][i]) != 1) return 2;
                    if (raw[j][i] < lo[i]) lo[i] = raw[j][i];
                    if (raw[j][i] > hi[i]) hi[i] = raw[j][i];
                }
            int c[MAXN][3], b[3], P[MAXN], img[MAXN], ref[MAXN];
            for (int i = 0; i < 3; i++) b[i] = hi[i] - lo[i] + 1;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < 3; i++) c[j][i] = raw[j][i] - lo[i];
                P[j] = KEY(c[j][0], c[j][1], c[j][2]);
            }
            isort(P, n);
            int aut = 0, achiral = 0, have = 0;
            for (int a = 0; a < 48; a++) {
                image_keys(&M48[a], c, b, n, img);
                isort(img, n);
                if (!lexcmp(img, P, n)) { aut++; if (M48[a].det < 0) achiral = 1; }
                if (!have || lexcmp(img, ref, n) < 0) { memcpy(ref, img, n * sizeof(int)); have = 1; }
            }
            TestOut o; memset(&o, 0, sizeof o);
            test_matrix(ident, c, b, n, 1, &o);
            int translate = o.found;
            for (int t = 0; t < 24; t++)
                if (prop24[t] != ident) test_matrix(prop24[t], c, b, n, 1, &o);
            set_shape_rows(c, n);
            TestOut ob; memset(&ob, 0, sizeof ob);   /* bitmask version, full count */
            for (int t = 0; t < 24; t++) test_matrix_bits(prop24[t], c, b, n, 1, &ob);
            TestOut of; memset(&of, 0, sizeof of);
            int twin_fast = locks(prop_red, nprop_red, c, b, n, 0, &of);
            TestOut om; memset(&om, 0, sizeof om);
            int mirror = o.found || locks(improp24, 24, c, b, n, 0, &om);
            printf("%d %d %d %llu %d %d %d %016llx %016llx %llu %d %d\n", o.found, twin_fast, mirror,
                   (unsigned long long)o.witnesses, translate, aut, achiral,
                   (unsigned long long)hash_ref(ref, n, SEED1), (unsigned long long)hash_ref(ref, n, SEED2),
                   (unsigned long long)ob.witnesses, ob.nocontact, o.nocontact);
        }
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "selftest")) {
        int ok = selftest(argv[2]);
        printf("selftest %s\n", ok ? "PASSED" : "FAILED");
        return ok ? 0 : 1;
    }
    if (argc >= 5 && !strcmp(argv[1], "mkshards")) {
        N = atoi(argv[2]); SD = atoi(argv[3]);
        if (SD < 1 || SD >= N) { fprintf(stderr, "need 1 <= D < N\n"); return 2; }
        init_grid(N);
        shard_out = fopen(argv[4], "w");
        if (!shard_out) return 2;
        int U[1] = {origin};
        gen(U, 1, 0);
        fclose(shard_out);
        printf("%llu shards\n", (unsigned long long)nshards);
        return 0;
    }
    if (argc >= 10 && !strcmp(argv[1], "run")) {
        int n = atoi(argv[2]);
        const char *shardfile = argv[3], *donefile = argv[4], *results = argv[5], *sampfile = argv[6];
        int worker = atoi(argv[7]), W = atoi(argv[8]);
        if (!selftest(argv[9])) { fprintf(stderr, "selftest FAILED, refusing to run\n"); return 1; }
        /* done shard ids */
        u64 total = 0;
        FILE *sf = fopen(shardfile, "r");
        if (!sf) return 2;
        for (int ch; (ch = fgetc(sf)) != EOF;) if (ch == '\n') total++;
        rewind(sf);
        u8 *done = calloc(total + 1, 1);
        FILE *df = fopen(donefile, "r");
        if (df) {
            unsigned long long id;
            while (fscanf(df, "%llu", &id) == 1) if (id < total) done[id] = 1;
            fclose(df);
        }
        FILE *rf = fopen(results, "a"), *spf = fopen(sampfile, "a");
        if (!rf || !spf) return 2;
        N = n; DEBUG = 0;
        init_grid(n);
        char *line = malloc(1 << 16);
        while (fgets(line, 1 << 16, sf)) {
            char *p = line, *q;
            unsigned long long id = strtoull(p, &q, 10); p = q;
            int D = (int)strtol(p, &q, 10); p = q;
            if (id % (unsigned)W != (unsigned)worker || done[id]) continue;
            memcpy(seen, seen_base, (size_t)GXY * (n + 2));
            for (int j = 0; j < D; j++) {
                int x = (int)strtol(p, &q, 10); p = q;
                int y = (int)strtol(p, &q, 10); p = q;
                int z = (int)strtol(p, &q, 10); p = q;
                int idx = idx_of(x, y, z);
                set_cell(j, idx);
                for (int e = 0; e < 6; e++) seen[idx + off[e]] = 1;
                seen[idx] = 1;
            }
            int k = (int)strtol(p, &q, 10); p = q;
            int U[6 * MAXN + 8];
            for (int j = 0; j < k; j++) {
                int x = (int)strtol(p, &q, 10); p = q;
                int y = (int)strtol(p, &q, 10); p = q;
                int z = (int)strtol(p, &q, 10); p = q;
                U[j] = idx_of(x, y, z);
            }
            memset(&K, 0, sizeof K);
            samples = spf;
            double t0 = cpu_seconds();
            rec(U, k, D);
            char sid[32];
            snprintf(sid, sizeof sid, "%llu", id);
            fflush(spf);
            print_counts(rf, "A", n, sid, total, cpu_seconds() - t0);
        }
        fclose(rf); fclose(spf); fclose(sf);
        return 0;
    }
    fprintf(stderr, "usage: see the comment at the top of twin.c\n");
    return 2;
}
