#!/usr/bin/env bash
log=${snakemake_log[0]}
# 获取Snakemake传递的参数
extra=${snakemake_params[extra]}
threads=${snakemake[threads]}

# 获取输入文件数量
input_files=("${snakemake_input[sample]}")
n=${#input_files[@]}

# 检查输入文件数量是否正确
if [ $n -ne 1 ] && [ $n -ne 2 ]; then
    echo "错误: input->sample必须有1个(单端)或2个(双端)元素。" >&2
    exit 1
fi

# 设置输入文件参数
if [ $n -eq 1 ]; then
    reads="--in1 ${input_files[0]}"
else
    reads="--in1 ${input_files[0]} --in2 ${input_files[1]}"
fi

# 设置输出文件参数
trimmed=""
if [ -n "${snakemake_output[trimmed]}" ]; then
    if [ $n -eq 1 ]; then
        trimmed="--out1 ${snakemake_output[trimmed]}"
    else
        trimmed_files=("${snakemake_output[trimmed]}")
        trimmed="--out1 ${trimmed_files[0]} --out2 ${trimmed_files[1]}"
    fi
fi



# 设置统计输出
html="--html ${snakemake_output[html]}"
json="--json ${snakemake_output[json]}"

# 执行fastp命令
(fastp --thread "$threads" \
    $extra \
    $reads \
    $trimmed \
    $json \
    $html) 2>&1 > "$log"