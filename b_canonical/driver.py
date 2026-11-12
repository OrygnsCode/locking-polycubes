"""Runs one n of the Rust program as shards on several cores, and picks up where it stopped.

    python driver.py <exe> <table> <nmax> <wupto> <outdir> [workers] [min_commit_gb] [min_disk_gb]

<table> comes from "<exe> shardtable <d> <table>". each shard runs as its own process. run
the same command again on the same <outdir> to resume. when every shard is done it writes
<outdir>/B_collected.jsonl, the input for a tally.

what it guards against:
  - <outdir> is tied to the table: its sha256 goes in <outdir>/TABLE_SHA256 on the first
    launch and is checked after that, so resuming with a different table can't quietly skip
    the wrong shards.
  - only shards with a done marker get collected, and the marker is written after the
    results, so a half-written file is never summed. a marked shard is never rerun.
  - a failing shard's stderr is kept in <outdir>/fail_<id>.err (the binary aborts on panic,
    so that message is all you get).
  - all children share one job object with a memory cap, to bound a fault.
  - a lock file stops two drivers on the same directory. free disk and free commit are
    checked first, and if free commit drops under the floor it waits instead of starting
    more.
  - TWINB_SHARD_TIMEOUT (default 7200 s) turns a stuck shard into a logged failure that the
    next launch retries.

children run at below-normal priority with no console window, and the driver only ever kills
process ids it started itself.
"""
import ctypes, hashlib, json, os, shutil, subprocess, sys, time

if len(sys.argv) < 6:
    sys.exit(__doc__)
exe, table, nmax, wupto, outdir = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
# windows won't start a process from a bare relative path to the exe (WinError 2), even
# though open() finds the file, so make it absolute
exe = os.path.abspath(exe)
workers = int(sys.argv[6]) if len(sys.argv) > 6 else 4
min_commit_gb = float(sys.argv[7]) if len(sys.argv) > 7 else 5.0
min_disk_gb = float(sys.argv[8]) if len(sys.argv) > 8 else 5.0

JOB_MEMORY_CAP = 2 << 30          # 2 GB across all children
# a shard only writes a few MB, the disk floor is just so a long run can't fill the drive
# so one stuck shard can't hold up the whole run forever. raise it if a real shard needs longer
SHARD_TIMEOUT_S = float(os.environ.get("TWINB_SHARD_TIMEOUT", "7200"))

os.makedirs(outdir, exist_ok=True)
LOG = open(os.path.join(outdir, "driver.log"), "a", buffering=1, newline="")


def log(msg):
    line = "%s %s" % (time.strftime("%H:%M:%S"), msg)
    LOG.write(line + "\n")
    print(line, flush=True)


# box state

class MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]


def free_commit_gb():
    st = MEMORYSTATUSEX()
    st.dwLength = ctypes.sizeof(st)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return st.ullAvailPageFile / (1 << 30)


# job object

class IO_COUNTERS(ctypes.Structure):
    _fields_ = [("ReadOperationCount", ctypes.c_ulonglong),
                ("WriteOperationCount", ctypes.c_ulonglong),
                ("OtherOperationCount", ctypes.c_ulonglong),
                ("ReadTransferCount", ctypes.c_ulonglong),
                ("WriteTransferCount", ctypes.c_ulonglong),
                ("OtherTransferCount", ctypes.c_ulonglong)]


class JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong),
                ("PerJobUserTimeLimit", ctypes.c_longlong),
                ("LimitFlags", ctypes.c_ulong),
                ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t),
                ("ActiveProcessLimit", ctypes.c_ulong),
                ("Affinity", ctypes.c_size_t),
                ("PriorityClass", ctypes.c_ulong),
                ("SchedulingClass", ctypes.c_ulong)]


class JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("BasicLimitInformation", JOBOBJECT_BASIC_LIMIT_INFORMATION),
                ("IoInfo", IO_COUNTERS),
                ("ProcessMemoryLimit", ctypes.c_size_t),
                ("JobMemoryLimit", ctypes.c_size_t),
                ("PeakProcessMemoryUsed", ctypes.c_size_t),
                ("PeakJobMemoryUsed", ctypes.c_size_t)]


JOB_OBJECT_LIMIT_JOB_MEMORY = 0x00000200
JobObjectExtendedLimitInformation = 9


def make_job(cap_bytes):
    k32 = ctypes.windll.kernel32
    # HANDLEs are pointer-sized, ctypes' default c_int would truncate them
    k32.CreateJobObjectW.restype = ctypes.c_void_p
    k32.OpenProcess.restype = ctypes.c_void_p
    h = k32.CreateJobObjectW(None, None)
    if not h:
        log("WARNING: CreateJobObject failed, running without a memory cap")
        return None
    info = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
    info.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_JOB_MEMORY
    info.JobMemoryLimit = cap_bytes
    ok = k32.SetInformationJobObject(ctypes.c_void_p(h), JobObjectExtendedLimitInformation,
                                     ctypes.byref(info), ctypes.sizeof(info))
    if not ok:
        log("WARNING: SetInformationJobObject failed, running without a memory cap")
        return None
    return h


def assign(job, pid):
    if job is None:
        return
    k32 = ctypes.windll.kernel32
    PROCESS_SET_QUOTA, PROCESS_TERMINATE = 0x0100, 0x0001
    hp = k32.OpenProcess(PROCESS_SET_QUOTA | PROCESS_TERMINATE, False, pid)
    if hp:
        k32.AssignProcessToJobObject(ctypes.c_void_p(job), ctypes.c_void_p(hp))
        k32.CloseHandle(ctypes.c_void_p(hp))


# preflight

# free space on the output drive
drive = (os.path.splitdrive(os.path.abspath(outdir))[0] or "C:").upper()
floor = min_disk_gb
free_gb = shutil.disk_usage(drive + os.sep).free / (1 << 30)
if free_gb < floor:
    sys.exit("free space on %s is %.1f GB, under the %.1f GB floor"
             % (drive, free_gb, floor))

# bind this directory to this table
tsha = hashlib.sha256(open(table, "rb").read()).hexdigest()
shafile = os.path.join(outdir, "TABLE_SHA256")
if os.path.exists(shafile):
    was = open(shafile).read().strip().split()[0]
    if was != tsha:
        sys.exit("TABLE MISMATCH: %s was built for table %s..., this table is %s...; "
                 "resuming would skip the wrong subtrees" % (outdir, was[:16], tsha[:16]))
else:
    with open(shafile, "w", newline="") as fh:
        fh.write("%s  %s\n" % (tsha, os.path.basename(table)))

# one driver per directory
lockpath = os.path.join(outdir, "driver.lock")
try:
    lockfd = os.open(lockpath, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(lockfd, ("pid %d at %s\n" % (os.getpid(), time.strftime("%H:%M:%S"))).encode())
    os.close(lockfd)
except FileExistsError:
    sys.exit("%s exists: another driver is working %s, or one died. "
             "Check, then delete it." % (lockpath, outdir))

nshards = sum(1 for l in open(table) if l.strip() and not l.startswith("#"))
log("table %s (%s...) has %d shards, nmax=%d wupto=%d workers=%d commit floor=%.1f GB "
    "shard timeout %.0f s  disk floor %.0f GB  free disk on %s %.0f GB"
    % (table, tsha[:16], nshards, nmax, wupto, workers, min_commit_gb,
       SHARD_TIMEOUT_S, floor, drive, free_gb))

todo = [i for i in range(nshards)
        if not os.path.exists(os.path.join(outdir, "done_%d" % i))]
log("%d of %d shards still to do" % (len(todo), nshards))

FLAGS = 0
if os.name == "nt":
    FLAGS = subprocess.BELOW_NORMAL_PRIORITY_CLASS | subprocess.CREATE_NO_WINDOW
job = make_job(JOB_MEMORY_CAP)
log("job object memory cap %.1f GB across all children" % (JOB_MEMORY_CAP / (1 << 30)))

running = {}
queue = list(todo)
done, failed = 0, []
launch_cpu = 0.0          # CPU of the shards THIS launch ran
t0 = time.time()
try:
    while queue or running:
        while queue and len(running) < workers:
            fc = free_commit_gb()
            if fc < min_commit_gb:
                if running:
                    break
                log("HOLD: free commit %.1f GB under the %.1f GB floor, no children alive"
                    % (fc, min_commit_gb))
                time.sleep(60)
                continue
            sid = queue.pop(0)
            errf = open(os.path.join(outdir, "err_%d.tmp" % sid), "wb")
            p = subprocess.Popen(
                [exe, "shard", str(nmax), str(wupto), table, str(sid), outdir],
                stdout=subprocess.DEVNULL, stderr=errf, creationflags=FLAGS)
            assign(job, p.pid)
            running[p] = [sid, errf, time.time(), False]
        time.sleep(0.2)
        for p in list(running):
            rc = p.poll()
            if rc is None:
                # terminate once, not on every poll while it dies
                if not running[p][3] and time.time() - running[p][2] > SHARD_TIMEOUT_S:
                    running[p][3] = True
                    el = time.time() - running[p][2]
                    log("shard %d exceeded %.0f s, terminating it (my own PID %d only)"
                        % (running[p][0], SHARD_TIMEOUT_S, p.pid))
                    # a timeout isn't a crash: write it down, so a run that can never
                    # finish at this timeout says so instead of looping
                    with open(os.path.join(outdir, "timeout_%d" % running[p][0]),
                              "w", newline="") as tf:
                        tf.write("%d exceeded %.0f s (ran %.0f s)\n"
                                 % (running[p][0], SHARD_TIMEOUT_S, el))
                    p.terminate()
                continue
            sid, errf = running.pop(p)[:2]
            errf.close()
            tmp = os.path.join(outdir, "err_%d.tmp" % sid)
            marked = os.path.exists(os.path.join(outdir, "done_%d" % sid))
            if rc != 0 or not marked:
                # keep the reason
                keep = os.path.join(outdir, "fail_%d.err" % sid)
                try:
                    os.replace(tmp, keep)
                except OSError:
                    pass
                failed.append(sid)
                tail = ""
                try:
                    tail = open(keep, errors="replace").read()[-300:].replace("\n", " | ")
                except OSError:
                    pass
                log("shard %d FAILED rc=%s marker=%s :: %s" % (sid, rc, marked, tail))
            else:
                try:
                    os.remove(tmp)
                except OSError:
                    pass
                done += 1
                # it finished, so forget any earlier timeout
                try:
                    os.remove(os.path.join(outdir, "timeout_%d" % sid))
                except OSError:
                    pass
                try:
                    launch_cpu += float(open(os.path.join(outdir, "done_%d" % sid))
                                        .read().split()[1])
                except (OSError, IndexError, ValueError):
                    pass
                if done % 100 == 0 or not queue:
                    log("%d/%d done, %d running, %d queued, %.0f s elapsed"
                        % (done, len(todo), len(running), len(queue), time.time() - t0))
except KeyboardInterrupt:
    log("interrupted, terminating only my own %d children" % len(running))
    for p in list(running):
        p.terminate()
    raise
finally:
    try:
        os.remove(lockpath)
    except OSError:
        pass

timed_out = sorted(int(f.split("_")[1]) for f in os.listdir(outdir)
                   if f.startswith("timeout_"))
log("finished: %d done, %d failed %s" % (done, len(failed), failed[:20] if failed else ""))
if timed_out:
    log("TIMED OUT, and these will fail again at this timeout: %s. "
        "Raise TWINB_SHARD_TIMEOUT above the largest shard's CPU time and rerun."
        % timed_out[:20])

collected = os.path.join(outdir, "B_collected.jsonl")
have, lines, cpu = [], [], 0.0
for i in range(nshards):
    if not os.path.exists(os.path.join(outdir, "done_%d" % i)):
        continue
    rp = os.path.join(outdir, "B_results_%d.jsonl" % i)
    if not os.path.exists(rp):
        log("shard %d has a marker but no results file" % i)
        continue
    have.append(i)
    first = True
    for line in open(rp):
        line = line.strip()
        if not line:
            continue
        lines.append(line)
        if first:
            cpu += json.loads(line)["seconds"]
            first = False
with open(collected, "w", newline="") as fh:
    fh.write("\n".join(lines) + "\n")
log("collected %d of %d shards into %s, %d lines, %.1f CPU-s in total"
    % (len(have), nshards, collected, len(lines), cpu))
if len(have) != nshards:
    log("INCOMPLETE: %d shards missing, rerun this driver to finish them"
        % (nshards - len(have)))
    sys.exit(1)
wall = time.time() - t0
if done == nshards:
    log("COMPLETE %d/%d, wall %.0f s, %.1f CPU-s this launch, parallel speedup %.2fx on %d workers"
        % (len(have), nshards, wall, launch_cpu, launch_cpu / max(wall, 1e-9), workers))
else:
    # this launch skipped finished shards, so only report what it actually ran
    log("COMPLETE %d/%d (this launch ran %d of them), wall %.0f s, %.1f CPU-s this launch, "
        "parallel speedup %.2fx on %d workers -- NOT a figure for the whole run"
        % (len(have), nshards, done, wall, launch_cpu,
           launch_cpu / max(wall, 1e-9), workers))
