#!/bin/bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
# aggregate.sh <sweep_dir> : per-key table + aggregates from attack logs
d=$1
printf "%-5s %-4s %12s %8s %8s %12s %10s\n" seed rec lookups hits rawhits ks_tries wall_s
for f in $(ls $d/key_*.log | sort -t_ -k2 -n); do
  seed=$(basename $f .log | sed 's/key_//')
  rec=$(grep -m1 "key recovered" $f | grep -q YES && echo YES || echo NO)
  lk=$(grep -m1 "lookups" $f | awk '{print $NF}')
  hits=$(grep -m1 "bloom/self hits" $f | awk '{print $NF}')
  raw=$(grep -m1 "raw probe hits" $f | awk '{print $NF}'); raw=${raw:-"-"}
  ks=$(grep -m1 "key-sched tries" $f | awk '{print $NF}')
  wall=$(grep -m1 "^  time " $f | awk '{print $3}')
  printf "%-5s %-4s %12s %8s %8s %12s %10s\n" $seed $rec $lk $hits $raw $ks $wall
done
echo "---- aggregates ----"
nrec=$(grep -l "key recovered   : YES" $d/key_*.log | wc -l); n=$(ls $d/key_*.log | wc -l)
echo "recovered: $nrec / $n"
for field in "lookups" "^  time " "bloom/self hits" "key-sched tries" "raw probe hits"; do
  case "$field" in "^  time ") lbl="wall_s"; col=3;; *) lbl=$(echo "$field"|tr -d '^'); col=NF;; esac
  grep -h -m1 "$field" $d/key_*.log | awk -v col=$col '{v=(col=="NF")?$NF:$(col); print v}' | sort -g | \
  awk -v lbl="$lbl" '{a[NR]=$1; s+=$1} END{if(!NR){exit} mean=s/NR; med=(NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2; printf "%-18s n=%d mean=%.4g median=%.4g min=%.4g max=%.4g total=%.6g\n", lbl, NR, mean, med, a[1], a[NR], s}'
done
