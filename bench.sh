#!/bin/bash
# Created on 2026.05.01
# Copyright Youcef Lemsafer

echo '| TF limit | Digits | Time(s) |'
echo '|----------|--------|---------|'

indices=(12221 22221 42221 62221 82221)
digits=(49999 99999 199999 299999 399999)

for i in {0..4};
do
    for bits in 33 34;
    do
        out=$(./bin-cu/cutrialdive --mode smarandache -s ${indices[i]} --tf-bits ${bits} | tee -a bench-run-output.txt)
        echo -n "| 2^${bits} | ${digits[i]} | "
        echo -n $(echo ${out} | sed -e 's@.*\[Factoring took \([0-9]\+\.[0-9]\+\).*@\1@g')
        echo " |"
    done
done
