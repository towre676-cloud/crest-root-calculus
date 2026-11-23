"""Post-study compiled-PARI comparison, preserving the original frozen study.

This protocol is frozen AFTER viewing the original GP comparisons. It reuses
ALL twenty inputs, with a contemporaneous rerun of Crest as the anchor. It is
not an independently selected confirmation panel.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import resource
import signal
import subprocess
import time

from benchmark import digest, output_digest, write_json, cpu_usage
from model import expand_case

HERE=Path(__file__).resolve().parent
METHODS=["crest-cpp","pari-native-linear","pari-native-doubling"]
SOURCES=["model.py","solver.cpp","pari_native.c","supplement.py","benchmark.py"]


def invoke(case, method, iterations, protocol, expected):
    a,b=expand_case(case)
    if method=="crest-cpp":
        command=[protocol["crest_path"]]
    else:
        command=[protocol["pari_native_path"],method.rsplit("-",1)[1]]
    source=f"{len(a)} {iterations}\n"+" ".join(map(str,a))+"\n"+" ".join(map(str,b))+"\n"
    def limits():
        os.setsid()
        resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
        os.sched_setaffinity(0,{protocol["worker_cpu"]})
    start_cpu,start=cpu_usage(),time.perf_counter()
    process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                             text=True,preexec_fn=limits)
    try:
        stdout,stderr=process.communicate(source,timeout=protocol["per_process_wall_limit_seconds"])
    except subprocess.TimeoutExpired:
        os.killpg(process.pid,signal.SIGKILL)
        stdout,stderr=process.communicate()
        return {"status":"timeout","iterations":iterations,
                "limit_seconds":protocol["per_process_wall_limit_seconds"],
                "process_wall_seconds":time.perf_counter()-start,"process_cpu_seconds":cpu_usage()-start_cpu,
                "stderr":stderr[-3000:]}
    wall,cpu=time.perf_counter()-start,cpu_usage()-start_cpu
    try:
        if process.returncode or "***" in stderr:
            raise RuntimeError("Native solver failed")
        raw=json.loads(stdout)
        values=[int(v) for v in raw.pop("coefficients")]
        actual_period=next(m for m in range(1,len(a)+1) if len(a)%m==0
                          and all(a[i]==a[i%m] and b[i]==b[i%m] for i in range(len(a))))
        status="complete" if values==expected and raw["period"]==actual_period else "wrong_answer"
        result={"status":status,**raw,"kernel_cpu_seconds":raw["cpu_seconds"],
                "process_wall_seconds":wall,"process_cpu_seconds":cpu,
                "output_sha256":output_digest(values)}
        if status!="complete":
            result["actual_coefficients"]=[str(v) for v in values]
            result["expected_period"]=actual_period
        return result
    except Exception as error:
        return {"status":"error","iterations":iterations,"reason":str(error),
                "process_wall_seconds":wall,"process_cpu_seconds":cpu,
                "stdout":stdout[-3000:],"stderr":stderr[-3000:]}


def run(args):
    target=args.output.resolve()
    protocol_path=target.with_name(target.stem+"-protocol.json")
    if target.exists() or protocol_path.exists():
        raise RuntimeError("Refusing to overwrite a recorded supplement")
    panel=json.loads((HERE/"panel.json").read_text())
    reference_path=HERE/"results-reference.json"
    reference={key:list(map(int,values)) for key,values in json.loads(reference_path.read_text()).items()}
    cases={c["id"]:c for c in panel["cases"]}
    schedule=[];rng=random.Random(76128)
    for replicate in range(1,args.replicates+1):
        order=[(c["id"],m) for c in panel["cases"] for m in METHODS];rng.shuffle(order)
        schedule.extend({"case":c,"method":m,"replicate":replicate} for c,m in order)
    protocol={"schema":1,"created_utc":datetime.now(timezone.utc).isoformat(),
              "original_protocol_sha256":digest(HERE/"results-protocol.json"),
              "original_results_sha256":digest(HERE/"results.json"),
              "panel_sha256":digest(HERE/"panel.json"),"reference_sha256":digest(reference_path),
              "source_sha256":{name:digest(HERE/name) for name in SOURCES},
              "crest_path":str(args.crest.resolve()),"crest_sha256":digest(args.crest),
              "pari_native_path":str(args.pari_native.resolve()),"pari_native_sha256":digest(args.pari_native),
              "pari_version":"2.19.1","pari_build":json.loads((HERE/"pari-build.json").read_text()),
              "platform":platform.platform(),"worker_cpu":args.cpu,"methods":METHODS,
              "schedule_seed":76128,"replicates":args.replicates,"schedule":schedule,
              "per_process_wall_limit_seconds":15,"minimum_batch_cpu_seconds":0.1,
              "maximum_batch_iterations":1048576,"address_space_bytes":4*1024**3,
              "design_status":"Supplement designed after seeing the original 240-run study, to test the native-versus-GP-language confound. The entire unchanged 20-input panel is retained; no speed-based input selection. Contemporaneous Crest anchor included. Not an independently selected confirmation panel.",
              "timing":"Identical common expanded-array input. Validation, period discovery and polynomial construction repeated inside every timed solve. Native CLOCK_PROCESS_CPUTIME_ID; batches at least 100 ms. Parsing, stack initialization and serialization outside kernel CPU but inside cold process wall/CPU. Every attempted calibration retained.",
              "compile_command":json.loads((HERE/"pari-native-build.json").read_text())["compile_command"]}
    write_json(protocol_path,protocol)
    raw={"schema":1,"completed":False,"protocol_sha256":digest(protocol_path),
         "reference_file":"results-reference.json","runs":[]}
    write_json(target,raw)
    for i,planned in enumerate(schedule):
        case=cases[planned["case"]];method=planned["method"]
        row=dict(planned)
        row["cold"]=invoke(case,method,1,protocol,reference[case["id"]])
        row["batches"]=[];row["status"]=row["cold"]["status"]
        if row["status"]=="complete":
            t=row["cold"]["kernel_cpu_seconds"]
            count=1 if t>=0.1 else min(1048576,max(2,math.ceil(0.12/max(t,0.0001))))
            while True:
                batch=invoke(case,method,count,protocol,reference[case["id"]]);row["batches"].append(batch)
                if batch["status"]!="complete":row["status"]=batch["status"];break
                if batch["kernel_cpu_seconds"]>=0.1:
                    row["kernel_seconds_per_solve"]=batch["kernel_cpu_seconds"]/count;break
                if count>=1048576:row["status"]="timer_floor";break
                count=min(1048576,2*count)
        raw["runs"].append(row);write_json(target,raw)
        print(f"{i+1}/{len(schedule)} {case['id']} {method} r{planned['replicate']} {row['status']}",flush=True)
    for name,sha in protocol["source_sha256"].items():
        if digest(HERE/name)!=sha:raise ArithmeticError("Supplement source changed: "+name)
    for name in ("crest","pari_native"):
        if digest(protocol[name+"_path"])!=protocol[name+"_sha256"]:raise ArithmeticError("Executable changed")
    raw["dependencies_unchanged"]=True;raw["completed"]=True
    raw["finished_utc"]=datetime.now(timezone.utc).isoformat();write_json(target,raw)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crest",type=Path,required=True)
    parser.add_argument("--pari-native",type=Path,required=True)
    parser.add_argument("--cpu",type=int,default=min(os.sched_getaffinity(0)))
    parser.add_argument("--replicates",type=int,default=3)
    parser.add_argument("--output",type=Path,default=HERE/"native-results.json")
    run(parser.parse_args())
