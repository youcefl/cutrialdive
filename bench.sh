#!/bin/bash
# Created on 2026.05.01
# Copyright Youcef Lemsafer

echo '| TF limit | Digits | Time(s) |'
echo '|----------|--------|---------|'
for i in 12221 22221 42221 62221 82221;
do
    out=$(./bin-cu/cutrialdive -i ${i} | tee -a bench-run-output.txt)
    echo -n "| "
    echo -n $(echo ${out} | sed -e 's@.*Trial factoring to \([0-9^]\+\).*@\1@g')
    echo -n " | "
    echo -n $(echo ${out} | sed -e 's@.* \([0-9]\+\) digits.*@\1@g')
    echo -n " | "
    echo -n $(echo ${out} | sed -e 's@.*\[Factoring took \([0-9]\+\.[0-9]\+\).*@\1@g')
    echo " |"
done
