# draws lock6.svg: the first 6-cell polycube (blue) locked with its copy (orange), seen from
# two opposite sides. standard library only, run it from this folder.

P = [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (0, 1, 2), (1, 1, 0)]
Q = [(-z + 1, -y + 1, -x + 1) for x, y, z in P]   # half turn (x,y,z) -> (-z,-y,-x), then +(1,1,1)

UNITS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


# check what the README says before drawing it
p, q = set(P), set(Q)
assert q == {(-1, 0, 1), (-1, 1, 1), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 1)}
assert not p & q
assert sum(add(c, u) in q for c in p for u in UNITS) == 8          # shared faces
assert all({add(c, u) for c in q} & p for u in UNITS)              # every step overlaps
# each piece is a C of five cubes plus one, and each C's gap holds a cube of the other piece
p_c, q_c = p - {(1, 1, 0)}, q - {(1, 0, 0)}
assert all(c[0] == 0 for c in p_c) and (0, 1, 1) in q
assert all(c[2] == 1 for c in q_c) and (0, 0, 1) in p
# the two C's stop five of the six steps; +y is stopped only by the two extra cubes
assert [u for u in UNITS if not {add(c, u) for c in q_c} & p_c] == [(0, 1, 0)]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def view(eye, ox, oy, size=46):
    # polygons for the pair seen from eye, as (points, fill), back to front
    right = cross((-eye[0], -eye[1], -eye[2]), (0, 0, 1))
    up = cross(right, (-eye[0], -eye[1], -eye[2]))
    rn, un = dot(right, right) ** 0.5, dot(up, up) ** 0.5

    def screen(v):
        return (ox + size * dot(v, right) / rn, oy - size * dot(v, up) / un)

    # x then y must turn counterclockwise on the page, or the drawing is mirrored
    (ax, ay), (bx, by), (cx, cy) = screen((0, 0, 0)), screen((1, 0, 0)), screen((0, 1, 0))
    assert (bx - ax) * (ay - cy) - (ay - by) * (cx - ax) > 0

    shade = {0: 0.80, 1: 0.65, 2: 1.0}       # side faces darker than tops
    colour = {"P": (66, 118, 190), "Q": (226, 128, 60)}
    both = {c: "P" for c in p} | {c: "Q" for c in q}
    polys = []
    for c in sorted(both, key=lambda c: dot(c, eye)):
        for u in UNITS:
            if dot(u, eye) <= 0 or add(c, u) in both:
                continue
            k = [i for i in range(3) if u[i]][0]
            a, b = [i for i in range(3) if i != k]
            base = list(c)
            base[k] += 1 if u[k] > 0 else 0
            corners = []
            for da, db in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                v = list(base)
                v[a] += da
                v[b] += db
                corners.append(screen(v))
            r, g, bl = (round(x * shade[k]) for x in colour[both[c]])
            polys.append((corners, "#%02x%02x%02x" % (r, g, bl)))
    return polys


# the two views side by side, trimmed to the drawing with a small margin
left, right = view((1, -1, 1), 0, 0), view((-1, 1, 1), 0, 0)


def box(polys):
    xs = [x for c, _ in polys for x, _ in c]
    ys = [y for c, _ in polys for _, y in c]
    return min(xs), max(xs), min(ys), max(ys)


m, gap = 16, 50
l, r = box(left), box(right)
h = max(l[3] - l[2], r[3] - r[2])
polys = [([(x - l[0] + m, y - l[2] + m) for x, y in c], f) for c, f in left]
polys += [([(x - r[0] + m + (l[1] - l[0]) + gap, y - r[2] + m) for x, y in c], f) for c, f in right]
W, H = round(2 * m + (l[1] - l[0]) + gap + (r[1] - r[0])), round(2 * m + h)
out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H)]
for corners, fill in polys:
    pts = " ".join("%.1f,%.1f" % xy for xy in corners)
    out.append('<polygon points="%s" fill="%s" stroke="#222222" stroke-width="1.3" '
               'stroke-linejoin="round"/>' % (pts, fill))
out.append("</svg>")
with open("lock6.svg", "w", newline="\n") as f:
    f.write("\n".join(out) + "\n")
