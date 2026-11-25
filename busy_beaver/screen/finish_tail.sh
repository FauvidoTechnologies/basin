#!/bin/sh
# Runs chunks 198 and 199 as 10 pieces each (instead of one core each), concatenates them in order,
# then waits for all 200 chunks and writes full/DONE.
cd "$(dirname "$0")"
N=$(( $(stat -c %s rest.int64) / 8 ))
for i in 198 199; do
  a=$((N*i/200)); b=$((N*(i+1)/200))
  for k in $(seq 0 9); do
    echo "$((a+(b-a)*k/10)) $((a+(b-a)*(k+1)/10)) ${i}_$k"
  done
done | xargs -P 20 -L 1 sh -c './screen rec ../bbdata/timelimit.bin rest.int64 $0 $1 26 full/p_$2.part && mv full/p_$2.part full/p_$2.bin'
for i in 198 199; do
  a=$((N*i/200)); b=$((N*(i+1)/200))
  cat $(for k in $(seq 0 9); do echo full/p_${i}_$k.bin; done) > full/c_$i.tmp
  [ $(stat -c %s full/c_$i.tmp) -eq $(( (b-a)*257*4 )) ] && mv full/c_$i.tmp full/c_$i.bin && rm -f full/p_${i}_*.bin || echo "size mismatch chunk $i"
done
while [ $(ls full/c_*.bin | wc -l) -lt 200 ]; do sleep 30; done
echo "finished $(date)" > full/DONE
