#!/bin/bash

for num in {2..21}
do
    # 创建临时文件来存储修改后的内容
    sed "s/_dir_/sterics_group$((num-1))/" tmp > temp_input.txt
    sed -i "s/_num_/${num}/g" temp_input.txt 
    cp temp_input.txt sterics_group${num}/input.txt
done

# 清理临时文件
rm temp_input.txt
