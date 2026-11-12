"""Runs twin.exe over every shard of one n, in parallel. Windows only (uses job objects).

    python driver.py N D W OUTDIR [MIN_COMMIT_GB]

N is the number of cells, D the shard depth (each search state at depth D is a shard), W the
number of workers. build twin.exe next to this file first. result files are named after the
first 8 hex digits of twin.exe's sha256, so a rebuilt exe never mixes its lines with an old one's.

if a run stops for any reason, run the same command again: shards that already have a complete
result line are skipped. workers run a copy of the exe from OUTDIR at below-normal priority with
no console window, each capped at 512 MB, and each one runs the selftest before its first shard.
the driver stops its own workers if the output drive gets under MIN_FREE_GB, and won't start
with less than MIN_COMMIT_GB of free commit (default 4).
"""
import ctypes, hashlib, json, msvcrt, os, shutil, subprocess, sys, time
from ctypes import wintypes

HERE = os.path.dirname(os.path.abspath(__file__))
REFLINES = os.path.join(HERE, "..", "data", "selftest_lines.txt")
CREATE_NO_WINDOW, BELOW_NORMAL = 0x08000000, 0x00004000
MIN_FREE_GB = 5
MEM_CAP = 512 * 1024 * 1024


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    with open(os.path.join(OUT, "driver.log"), "a", newline="") as f:
        f.write(line + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


class MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [("dwLength", wintypes.DWORD), ("dwMemoryLoad", wintypes.DWORD)] + [
        (n, ctypes.c_uint64) for n in ("ullTotalPhys", "ullAvailPhys", "ullTotalPageFile", "ullAvailPageFile",
                                       "ullTotalVirtual", "ullAvailVirtual", "ullAvailExtendedVirtual")]


def free_commit_gb():
    m = MEMORYSTATUSEX()
    m.dwLength = ctypes.sizeof(m)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPageFile / 2**30


def free_gb(path):
    return shutil.disk_usage(os.path.splitdrive(os.path.abspath(path))[0] + "\\").free / 2**30


class JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64), ("PerJobUserTimeLimit", ctypes.c_int64),
                ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD),
                ("SchedulingClass", wintypes.DWORD)]


class IO_COUNTERS(ctypes.Structure):
    _fields_ = [(n, ctypes.c_uint64) for n in ("ReadOperationCount", "WriteOperationCount",
                "OtherOperationCount", "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("BasicLimitInformation", JOBOBJECT_BASIC_LIMIT_INFORMATION), ("IoInfo", IO_COUNTERS),
                ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]


def capped_job():
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.CreateJobObjectW.restype = wintypes.HANDLE
    job = k32.CreateJobObjectW(None, None)
    info = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
    info.BasicLimitInformation.LimitFlags = 0x100  # JOB_OBJECT_LIMIT_PROCESS_MEMORY
    info.ProcessMemoryLimit = MEM_CAP
    if not k32.SetInformationJobObject(job, 9, ctypes.byref(info), ctypes.sizeof(info)):
        raise OSError(ctypes.get_last_error(), "SetInformationJobObject")
    return k32, job


def fail(msg):
    log("REFUSED: " + msg)
    sys.exit(1)


def image_path(pid):
    # full path of the exe running as pid, or None
    k32 = ctypes.windll.kernel32
    k32.OpenProcess.restype = wintypes.HANDLE
    h = k32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
    if not h:
        return None
    try:
        code = wintypes.DWORD()
        if not k32.GetExitCodeProcess(h, ctypes.byref(code)) or code.value != 259:  # STILL_ACTIVE
            return None
        buf = ctypes.create_unicode_buffer(1024)
        size = wintypes.DWORD(1024)
        if not k32.QueryFullProcessImageNameW(h, 0, buf, ctypes.byref(size)):
            return None
        return buf.value
    finally:
        k32.CloseHandle(h)


def terminate_partial_lines():
    # a worker killed mid-write can leave a line with no newline, add one so the next line isn't
    # glued on
    for name in os.listdir(OUT):
        if name.startswith(("results_", "samples_")) and name.endswith(".jsonl"):
            path = os.path.join(OUT, name)
            if os.path.getsize(path) == 0:
                continue
            with open(path, "rb+") as f:
                f.seek(-1, 2)
                if f.read(1) != b"\n":
                    f.write(b"\n")
                    log("closed a partial last line in " + name)


def done_ids():
    ids = set()
    for name in os.listdir(OUT):
        if name.startswith("results_%s_w" % TAG) and name.endswith(".jsonl"):
            for line in open(os.path.join(OUT, name)):
                if not line.endswith("\n"):
                    continue
                try:
                    ids.add(int(json.loads(line)["shard"]))
                except (ValueError, KeyError):
                    pass
    return ids


def main():
    exe_src = os.path.join(HERE, "twin.exe")
    src = os.path.join(HERE, "twin.c")
    if not os.path.exists(exe_src):
        fail("build twin.exe next to driver.py first (see the README)")
    exe = os.path.join(OUT, "twin_%s.exe" % TAG)
    # one driver per OUTDIR: hold a lock the whole time, and refuse if old workers are still running
    lockf = open(os.path.join(OUT, "driver.lock"), "a+")
    try:
        msvcrt.locking(lockf.fileno(), msvcrt.LK_NBLCK, 1)
    except OSError:
        fail("another driver holds driver.lock in this OUTDIR")
    globals()["_lock"] = lockf
    pids_path = os.path.join(OUT, "pids.txt")
    if os.path.exists(pids_path):
        for line in open(pids_path):
            parts = line.split()
            if parts and parts[0] == "worker":
                pid = int(parts[2])
                img = image_path(pid)
                if img and os.path.normcase(os.path.abspath(img)) == os.path.normcase(exe):
                    fail("worker %s (pid %d) from an earlier launch is still running %s" % (parts[1], pid, img))
    manifest_path = os.path.join(OUT, "manifest.json")
    if os.path.exists(manifest_path):
        man = json.load(open(manifest_path))
        if man["n"] != N or man["depth"] != D or man["tag"] != TAG:
            fail("manifest is for a different n, depth or build")
        if sha256(exe) != man["exe_sha256"]:
            fail("the exe in OUTDIR changed since the manifest was written")
        if sha256(REFLINES) != man["reflines_sha256"]:
            fail("reference lines changed since the manifest was written (update the manifest deliberately)")
    else:
        shutil.copy2(exe_src, exe)
        man = {"n": N, "depth": D, "tag": TAG, "exe_sha256": sha256(exe), "source_sha256": sha256(src),
               "reflines_sha256": sha256(REFLINES), "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    shards = os.path.join(OUT, "shards_n%d_d%d.txt" % (N, D))
    if not os.path.exists(shards):
        tmp = shards + ".tmp"
        subprocess.run([exe, "mkshards", str(N), str(D), tmp], check=True, creationflags=CREATE_NO_WINDOW)
        os.replace(tmp, shards)
        log("shard table written")
    total = sum(1 for _ in open(shards))
    shards_hash = sha256(shards)
    if "shards_sha256" in man and (man["shards_total"] != total or man["shards_sha256"] != shards_hash):
        fail("shard table changed since the manifest was written")
    man["shards_total"] = total
    man["shards_sha256"] = shards_hash
    with open(manifest_path, "w", newline="") as f:
        json.dump(man, f, indent=1)
        f.write("\n")
    terminate_partial_lines()
    done = done_ids()
    with open(os.path.join(OUT, "done.txt"), "w", newline="") as f:
        f.writelines("%d\n" % i for i in sorted(done))
    log("launch: n=%d depth=%d shards=%d done=%d workers=%d free=%.0f GB" % (N, D, total, len(done), W, free_gb(OUT)))
    if len(done) == total:
        log("nothing to do")
        return
    if free_gb(OUT) < MIN_FREE_GB:
        fail("free space below %d GB" % MIN_FREE_GB)
    log("free commit %.1f GB, floor %.1f GB" % (free_commit_gb(), MIN_COMMIT))
    if free_commit_gb() < MIN_COMMIT:
        fail("free commit %.1f GB is below %.1f GB" % (free_commit_gb(), MIN_COMMIT))
    k32, job = capped_job()
    k32.OpenProcess.restype = wintypes.HANDLE
    procs = []
    for k in range(W):
        lf = open(os.path.join(OUT, "worker_w%d.log" % k), "a")
        p = subprocess.Popen([exe, "run", str(N), shards, os.path.join(OUT, "done.txt"),
                              os.path.join(OUT, "results_%s_w%d.jsonl" % (TAG, k)), os.path.join(OUT, "samples_%s_w%d.jsonl" % (TAG, k)),
                              str(k), str(W), REFLINES],
                             stdout=lf, stderr=subprocess.STDOUT, creationflags=CREATE_NO_WINDOW | BELOW_NORMAL)
        h = k32.OpenProcess(0x1F0FFF, False, p.pid)
        if not h or not k32.AssignProcessToJobObject(job, h):
            log("warning: could not put worker %d (pid %d) in the memory-capped job" % (k, p.pid))
        procs.append(p)
    with open(os.path.join(OUT, "pids.txt"), "w", newline="") as f:
        f.write("driver %d\n" % os.getpid() + "".join("worker %d %d\n" % (k, p.pid) for k, p in enumerate(procs)))
    log("workers: " + " ".join(str(p.pid) for p in procs))
    while any(p.poll() is None for p in procs):
        time.sleep(30)
        n_done = len(done_ids())
        log("progress %d/%d shards, free %.0f GB, running %d" % (n_done, total, free_gb(OUT),
                                                                    sum(p.poll() is None for p in procs)))
        if free_gb(OUT) < MIN_FREE_GB:
            log("free space below %d GB: stopping own workers" % MIN_FREE_GB)
            for p in procs:
                if p.poll() is None:
                    p.terminate()
    codes = [p.returncode for p in procs]
    got = len(done_ids())
    log("workers exited %s; shards with results %d/%d" % (codes, got, total))
    if got == total and all(c == 0 for c in codes):
        log("# DRIVER COMPLETE %d/%d" % (got, total))
    else:
        log("# DRIVER INCOMPLETE %d/%d exit codes %s" % (got, total, codes))


if __name__ == "__main__":
    N, D, W, OUT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), os.path.abspath(sys.argv[4])
    MIN_COMMIT = float(sys.argv[5]) if len(sys.argv) > 5 else 4.0
    TAG = sha256(os.path.join(HERE, "twin.exe"))[:8] if os.path.exists(os.path.join(HERE, "twin.exe")) else "none"
    os.makedirs(OUT, exist_ok=True)
    main()
