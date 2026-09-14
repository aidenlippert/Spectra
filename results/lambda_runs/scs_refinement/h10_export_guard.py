"""Give the existing H10 exact exporter a 2100s whole-process ceiling."""
import os,signal,time,json
from pathlib import Path
solver_pid=16627;timer_pid=16626;cap=2100

def info(pid):
    text=Path(f'/proc/{pid}/stat').read_text();parts=text[text.rindex(')')+2:].split()
    return parts[0],int(parts[1]),int(parts[19])
def interrupted(signum,frame):raise SystemExit(signum)
for sig in (signal.SIGTERM,signal.SIGHUP,signal.SIGINT):signal.signal(sig,interrupted)
state,parent,start_ticks=info(solver_pid)
assert parent==timer_pid
assert b'spin_invariant_discovery' in Path(f'/proc/{solver_pid}/cmdline').read_bytes()
assert Path(f'/proc/{timer_pid}/cmdline').read_bytes().startswith(b'timeout\x001500s\x00')
start=start_ticks/os.sysconf('SC_CLK_TCK');outcome='running'
try:
    os.kill(timer_pid,signal.SIGSTOP)
    print(json.dumps({'action':'export_guard_started','solver_pid':solver_pid,'whole_process_cap_seconds':cap}),flush=True)
    while True:
        try:state,parent,ticks=info(solver_pid)
        except FileNotFoundError:outcome='completed';break
        if ticks!=start_ticks or state=='Z':outcome='completed';break
        uptime=float(Path('/proc/uptime').read_text().split()[0])
        if uptime-start>=cap:
            os.kill(solver_pid,signal.SIGTERM);outcome='whole_process_cap_reached';break
        time.sleep(2)
finally:
    try:os.kill(timer_pid,signal.SIGCONT)
    except ProcessLookupError:pass
    print(json.dumps({'action':'export_guard_finished','outcome':outcome}),flush=True)
