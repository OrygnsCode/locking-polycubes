# counts free polycubes with n cells that lock with a rotated copy of themselves, for small n.
#
# P locks with a copy Q = R(P) + v (R one of the 24 rotations, identity included, v an integer
# vector) if
#   1. Q doesn't overlap P,
#   2. Q shares at least one face with P, and
#   3. Q is blocked in all six axis directions: moving it some whole number of steps that way
#      makes it overlap P.
# polycubes are counted up to rotation and reflection (free polycubes, A038119).
#
#   python lock_small.py [nmax]        default nmax = 7, prints n, A038119(n), a(n)
#
# it just checks the definition cell by cell, no attempt at speed. n = 8 takes under a minute and
# each n after that about eight times longer.
import sys
from itertools import permutations, product

UNITS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


def signed_permutations():
    # all 48 symmetries of the cube as (perm, signs, determinant)
    out = []
    for perm in permutations(range(3)):
        inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
        for signs in product((1, -1), repeat=3):
            det = (-1) ** inversions * signs[0] * signs[1] * signs[2]
            out.append((perm, signs, det))
    return out


ALL48 = signed_permutations()
ROTATIONS = [m for m in ALL48 if m[2] == 1]


def image(m, cell):
    perm, signs, _ = m
    return tuple(signs[i] * cell[perm[i]] for i in range(3))


def canonical(cells):
    # smallest translated image under the 48 symmetries, so one name per free polycube
    best = None
    for m in ALL48:
        img = [image(m, c) for c in cells]
        low = [min(c[i] for c in img) for i in range(3)]
        key = tuple(sorted(tuple(c[i] - low[i] for i in range(3)) for c in img))
        if best is None or key < best:
            best = key
    return best


def grow(level):
    # all free polycubes with one more cell
    bigger = set()
    for shape in level:
        cells = set(shape)
        for c in shape:
            for u in UNITS:
                d = (c[0] + u[0], c[1] + u[1], c[2] + u[2])
                if d not in cells:
                    bigger.add(canonical(list(cells) + [d]))
    return bigger


def locks(shape):
    P = set(shape)
    # for each axis line: the smallest and largest coordinate P has on it
    extent = [{}, {}, {}]
    for c in shape:
        for a in range(3):
            line = c[:a] + c[a + 1:]
            lo, hi = extent[a].get(line, (c[a], c[a]))
            extent[a][line] = (min(lo, c[a]), max(hi, c[a]))
    for m in ROTATIONS:
        turned = [image(m, c) for c in shape]
        # condition 2: some cell of Q is a face neighbour of some cell of P
        offsets = {(p[0] + u[0] - t[0], p[1] + u[1] - t[1], p[2] + u[2] - t[2])
                   for p in shape for u in UNITS for t in turned}
        for v in offsets:
            Q = [(t[0] + v[0], t[1] + v[1], t[2] + v[2]) for t in turned]
            if any(q in P for q in Q):                      # condition 1
                continue
            ok = True
            for a in range(3):                              # condition 3, both directions of axis a
                ahead = behind = False
                for q in Q:
                    e = extent[a].get(q[:a] + q[a + 1:])
                    if e is not None:
                        if q[a] < e[1]:
                            ahead = True                    # a cell of P lies further along +a
                        if q[a] > e[0]:
                            behind = True
                if not (ahead and behind):
                    ok = False
                    break
            if ok:
                return True
    return False


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    level = {canonical([(0, 0, 0)])}
    for n in range(1, nmax + 1):
        if n > 1:
            level = grow(level)
        print(n, len(level), sum(locks(s) for s in level), flush=True)


if __name__ == "__main__":
    main()
