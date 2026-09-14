"""Bounded computational profiling; all counts describe this implementation.

The profiled task includes fixture construction, all seven methods, both
studies, summaries and JSON serialization. It excludes imports, proof writing,
earlier research attempts and verification reruns; those exclusions are explicit.
"""
from pathlib import Path
import cProfile
import json
import pstats
from time import perf_counter, process_time
import tracemalloc
from experiments.v5_run import run
from experiments.v5_complementarity import run as retention

ROOT=Path(__file__).resolve().parents[1]


def profile():
    profiler=cProfile.Profile()
    tracemalloc.start()
    wall=perf_counter(); cpu=process_time()
    profiler.enable()
    result=run()
    comp=retention()
    payload=json.dumps({'transfer':result,'retention':comp})
    profiler.disable()
    cpu=process_time()-cpu; wall=perf_counter()-wall
    _,peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
    stats=pstats.Stats(profiler)
    calls=[]
    for (path,line,function),(primitive,total,own,cumulative,_) in stats.stats.items():
        if '/experiments/v5_' in path or path.endswith('/fractions.py'):
            calls.append({'file':path.split('/')[-1],'line':line,'function':function,
                          'primitive_calls':primitive,'total_calls':total,
                          'own_seconds':own,'cumulative_seconds':cumulative})
    out={'transfer_worlds':64,'retention_worlds':128,'elapsed_seconds_with_profiling':wall,
         'cpu_seconds_with_profiling':cpu,'peak_traced_python_bytes':peak,
         'serialized_result_bytes':len(payload.encode()),'profile_total_calls':stats.total_calls,
         'profile_primitive_calls':stats.prim_calls,'selected_function_counts':sorted(calls,key=lambda x:(x['file'],x['line'])),
         'scope':'fixture generation, acquisition, reuse, all baselines, summaries and serialization',
         'exclusions':'imports, native allocator memory, proof construction, previous research and verification reruns',
         'not_claimed':'hardware energy, physical elapsed time, algorithm-independent bit complexity, or speedup over stateful Bayes'}
    (ROOT/'results/v5/computational_accounting.json').write_text(json.dumps(out,indent=2)+'\n')
    return {k:v for k,v in out.items() if k!='selected_function_counts'}


if __name__=='__main__': print(json.dumps(profile(),indent=2))
