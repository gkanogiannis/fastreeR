#!/bin/bash

fastreeR=bin/fastreeR.py
vcf2dis=bin/VCF2Dis_multi

MEM=4096

OUTDIR=benchmark_logs && mkdir -p $OUTDIR
LOGDIR=benchmark_logs && mkdir -p $LOGDIR

for t in 64 32 16 8 4 2; do
for v in `ls vcf/subsets/M*.vcf.gz`; do

  NAME=`basename ${v} .vcf.gz`
  OUT=${OUTDIR}/${NAME}_t$t.fastreer.out
  TXT=${OUTDIR}/${NAME}_t$t.fastreer.txt
  LOG=${LOGDIR}/${NAME}_t$t.fastreer.log
  echo "Running fastreeR on $v with $t cores ..."
  /usr/bin/time -f "Time: %E\nMemory: %M KB" -o $LOG \
    ${fastreeR} --mem ${MEM} VCF2DIST -i $v -o $OUT --threads $t \
      1>&2 | tee ${TXT}

  OUT=${OUTDIR}/${NAME}_t$t.vcf2dis.out
  TXT=${OUTDIR}/${NAME}_t$t.vcf2dis.txt
  LOG=${LOGDIR}/${NAME}_t$t.vcf2dis.log
  echo "Running VCF2DIS on $v with $t cores ..."
  /usr/bin/time -f "Time: %E\nMemory: %M KB" -o $LOG \
    ${vcf2dis} -InPut $v -OutPut $OUT -Thread $t\
      1>&2 | tee ${TXT}
done
done
