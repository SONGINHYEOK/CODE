
#단일 프로세스 압축해제


#!/bin/bash

StartTime=$(date +%s)

for i in {1..8}
do
        tar -zxvf conf_420.tar.gz
done

EndTime=$(date +%s)

echo "It takes $(($EndTime - $StartTime)) seconds to complete this task."



#range 범위로 설정해서 병렬


#!/bin/bash

for i in {1..8}
do
    find . -type f -maxdepth 1 -name conf_$i.tar.gz | xargs -I {} -P 8 tar -zxvf {}
done