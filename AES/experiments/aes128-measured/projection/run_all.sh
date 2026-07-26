#!/bin/bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
# run_all.sh -- pinned timing runs on the benchmark machine 
# usage: ./run_all.sh <offline-core> <online-core>
set -u
D="${WORKDIR:-$(pwd)}"
cd $D
OC=${1:-90}
NC=${2:-30}
LOG=$D/logs
mkdir -p $LOG

sibl(){ cat /sys/devices/system/cpu/cpu$1/topology/thread_siblings_list; }
echo "=== run_all start $(date -u) host $(hostname) ===" | tee $LOG/run_all.log
echo "offline core $OC (siblings $(sibl $OC)), online core $NC (siblings $(sibl $NC))" | tee -a $LOG/run_all.log
cat /proc/loadavg | tee -a $LOG/run_all.log
lscpu | grep -E "Model name|MHz|NUMA node[01]" | tee -a $LOG/run_all.log
gcc --version | head -1 | tee -a $LOG/run_all.log

# background per-CPU idleness sampler for our cores + siblings (node-wide /proc/stat)
python3 - "$OC" "$NC" > $LOG/coreidle.log 2>&1 <<'EOF' &
import sys,time
oc,nc=int(sys.argv[1]),int(sys.argv[2])
def sib(c):
    s=open(f"/sys/devices/system/cpu/cpu{c}/topology/thread_siblings_list").read().strip()
    return [int(x) for x in s.replace('-',',').split(',')]
cores=sorted(set(sib(oc)+sib(nc)))
def snap():
    d={}
    for l in open("/proc/stat"):
        if l.startswith("cpu") and not l.startswith("cpu "):
            p=l.split(); d[int(p[0][3:])]=(int(p[4]), sum(map(int,p[1:9])))
    return d
a=snap()
while True:
    time.sleep(10)
    b=snap()
    la=open("/proc/loadavg").read().split()[0]
    tot=0.0
    for c in sorted(a):
        dt=b[c][1]-a[c][1]
        if dt: tot+=1-(b[c][0]-a[c][0])/dt
    line=time.strftime("%H:%M:%S")+f" loadavg {la} busy_cores_total {tot:6.2f} |"
    for c in cores:
        dt=b[c][1]-a[c][1]; u=1-(b[c][0]-a[c][0])/dt if dt else 0
        line+=f" cpu{c} util {u:4.2f}"
    print(line,flush=True)
    a=b
EOF
SAMPLER=$!
echo "sampler pid $SAMPLER" | tee -a $LOG/run_all.log

# ---- ONLINE (background, on core NC) ----
( echo "=== online3 start $(date -u) core $NC ==="; cat /proc/loadavg;
  taskset -c $NC ./online3 0xA15BE2C4 20 16 5 18;
  echo "=== online3 end $(date -u) ==="; cat /proc/loadavg ) > $LOG/online3.log 2>&1 &
ONPID=$!
echo "online pid $ONPID" | tee -a $LOG/run_all.log

# ---- OFFLINE on core OC ----
{
 echo "=== offline gate $(date -u) ==="; taskset -c $OC ./offline3 gate
 echo "=== A full (stage2) ==="; cat /proc/loadavg; taskset -c $OC ./offline3 bench A 16 5 20 2
 echo "=== A gen-only (stage1) ==="; taskset -c $OC ./offline3 bench A 16 5 20 1
 echo "=== B full (stage2) ==="; cat /proc/loadavg; taskset -c $OC ./offline3 bench B 16 5 20 2
 echo "=== B gen-only (stage1) ==="; taskset -c $OC ./offline3 bench B 16 5 20 1
 echo "=== C full (stage2) ==="; cat /proc/loadavg; taskset -c $OC ./offline3 bench C 16 5 20 2
 echo "=== C gen-only (stage1) ==="; taskset -c $OC ./offline3 bench C 16 5 20 1
 echo "=== offline done $(date -u) ==="; cat /proc/loadavg
} > $LOG/offline3.log 2>&1

# ---- DENOMINATORS on core OC (after offline finished) ----
{
 echo "=== denom $(date -u) core $OC ==="; cat /proc/loadavg
 taskset -c $OC ./denom 100000000 5
 echo "=== denom done $(date -u) ==="; cat /proc/loadavg
} > $LOG/denom.log 2>&1

wait $ONPID
kill $SAMPLER 2>/dev/null
echo "=== run_all done $(date -u) ===" | tee -a $LOG/run_all.log
cat /proc/loadavg | tee -a $LOG/run_all.log
