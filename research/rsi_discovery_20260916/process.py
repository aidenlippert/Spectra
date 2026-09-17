"""Serial, deadline-limited children; kill only our own process group."""
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def run(command, directory, timeout, stdin=None, memory_mb=None, limit_files=True):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1")
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = "1"

    def limits():
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        resource.setrlimit(resource.RLIMIT_CPU, (int(timeout) + 1, int(timeout) + 2))
        if limit_files:
            resource.setrlimit(resource.RLIMIT_FSIZE, (128 * 1024**2, 128 * 1024**2))
        if memory_mb and sys.platform != "darwin":
            resource.setrlimit(resource.RLIMIT_AS, (memory_mb * 1024**2,) * 2)

    start = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    timed_out = False
    with (directory / "stdout").open("wb") as out, (directory / "stderr").open("wb") as err:
        p = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.PIPE if stdin else subprocess.DEVNULL,
                             stdout=out, stderr=err, start_new_session=True, preexec_fn=limits)
        try:
            p.communicate(stdin.encode() if stdin else None, timeout=timeout)
        except BaseException as error:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
            if isinstance(error, subprocess.TimeoutExpired):
                timed_out = True
            else:
                raise
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {"exit_code": p.returncode, "timeout": timed_out,
            "wall_seconds": time.monotonic() - start,
            "cpu_seconds": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
            # ru_maxrss is a process-family high-water mark, not a per-attempt peak.
            "process_family_peak_rss_bytes": after.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
            "rss_scope": "children high-water mark since controller start; not additive",
            "stdout": str(directory / "stdout"), "stderr": str(directory / "stderr")}
