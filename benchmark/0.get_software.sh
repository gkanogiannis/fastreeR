#!/bin/bash

vcf2dis_VER=1.54
fastreer_VER=2.0.0

mkdir -p ${PWD}/software && pushd ${PWD}/software

wget https://github.com/hewm2008/VCF2Dis/archive/refs/tags/v${vcf2dis_VER}.tar.gz
tar zxvf v${vcf2dis_VER}.tar.gz && rm v${vcf2dis_VER}.tar.gz
chmod +x VCF2Dis-${vcf2dis_VER}/bin/*

wget https://github.com/gkanogiannis/fastreeR/archive/refs/tags/${fastreer_VER}.tar.gz
tar zxvf ${fastreer_VER}.tar.gz && rm ${fastreer_VER}.tar.gz
chmod +x fastreeR-${fastreer_VER}/fastreeR.py

popd

mkdir -p ${PWD}/bin && pushd ${PWD}/bin
ln -sfn ../software/VCF2Dis-${vcf2dis_VER}/bin/* .
ln -sfn ../software/fastreeR-${fastreer_VER}/fastreeR.py .
ln -sfn ../software/fastreeR-${fastreer_VER}/inst .
popd
