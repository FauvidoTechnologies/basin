#!/bin/sh
# Resumable full screen: 200 chunks of rest.int64, 20 at a time; a chunk is done when full/c_<i>.bin
# exists (written to a .part file first, renamed when complete).
cd "$(dirname "$0")"
N=$(( $(stat -c %s rest.int64) / 8 ))
for i in $(seq 0 199); do
  [ -f full/c_$i.bin ] || echo "$((N*i/200)) $((N*(i+1)/200)) $i"
done | xargs -P $(nproc) -L 1 sh -c './screen rec ../bbdata/timelimit.bin rest.int64 $0 $1 26 full/c_$2.part && mv full/c_$2.part full/c_$2.bin'
echo "finished $(date)" > full/DONE
