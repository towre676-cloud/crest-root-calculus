"""Time the complete certificate and the unrestricted incidence on its positive chart."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import sympy as sp

HERE = Path(__file__).resolve().parent


def address_limit():
    import resource
    cap = 1024**3
    resource.setrlimit(resource.RLIMIT_AS,(cap,cap))


def measure(seconds):
    env = {**os.environ, "OPENBLAS_NUM_THREADS":"1"}
    source = str(HERE/"quartic.json")
    start = time.perf_counter()
    result = subprocess.run([sys.executable,"-B",str(HERE/"branch_fibres.py"),source],
                            capture_output=True,text=True,check=True,env=env)
    answer = json.loads(result.stdout)
    stats = {"python":platform.python_version(),"sympy":sp.__version__,
             "positive_status":answer["status"],"certificate_seconds":round(time.perf_counter()-start,3),
             "raw_chart":2,"raw_wall_budget_seconds":seconds,"raw_address_limit_MiB":None}
    kwargs = {}
    if sys.platform.startswith("linux"):
        kwargs["preexec_fn"] = address_limit
        stats["raw_address_limit_MiB"] = 1024
    start = time.perf_counter()
    command = [sys.executable,"-B",str(HERE.parent/"degree-four-incidence/incidence.py"),source,"2","--solve"]
    try:
        raw = subprocess.run(command,timeout=seconds,env=env,capture_output=True,text=True,**kwargs)
        stats.update(raw_status="finished" if raw.returncode==0 else "error",
                     raw_returncode=raw.returncode,raw_output=(raw.stdout if raw.returncode==0 else raw.stderr)[-1800:])
    except subprocess.TimeoutExpired:
        stats["raw_status"] = "wall budget reached"
    stats["raw_seconds"] = round(time.perf_counter()-start,3)
    return stats


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=int,default=90)
    args=parser.parse_args()
    if args.seconds<1:
        parser.error("Choose a positive wall-clock budget")
    print(json.dumps(measure(args.seconds),indent=2))
