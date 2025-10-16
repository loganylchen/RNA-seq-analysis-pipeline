#!/usr/bin/env bash

# 获取Snakemake参数
species=$(echo "${snakemake_params[species]}" | tr '[:upper:]' '[:lower:]')
build=${snakemake_params[build]}
release=${snakemake_params[release]}
gtf_release=$release
output_file=${snakemake_output[0]}
log=${snakemake_log[0]}

# 检查输出文件格式
if [[ "$output_file" == *.gz ]]; then
    out_gz=true
    out_fmt=$(basename "$output_file" | sed 's/\.[^.]*$//' | sed 's/.*\.//')
else
    out_gz=false
    out_fmt=$(basename "$output_file" | sed 's/.*\.//')
fi

# 设置branch
branch=""
if [ "$build" = "GRCh37" ]; then
    if [ "$release" -ge 81 ]; then
        branch="grch37/"
    fi
    if [ "$release" -gt 87 ]; then
        gtf_release=87
    fi
elif [ -n "${snakemake_params[branch]}" ]; then
    branch="${snakemake_params[branch]}/"
fi

# 设置flavor
flavor=${snakemake_params[flavor]:-""}
if [ -n "$flavor" ]; then
    flavor="${flavor}."
fi

# 检查输出格式并设置后缀
if [ "$out_fmt" = "gtf" ]; then
    suffix="gtf.gz"
elif [ "$out_fmt" = "gff3" ]; then
    suffix="gff3.gz"
else
    echo "错误：无效的格式。目前只支持 'gtf[.gz]' 和 'gff3[.gz]'。" >&2
    exit 1
fi

# 构建URL
url=${snakemake_params[url]:-"ftp://ftp.ensembl.org/pub"}
url="${url}/${branch}release-${release}/${out_fmt}/${species}/${species^}.${build}.${gtf_release}.${flavor}${suffix}"

# 下载文件
if $out_gz; then
    if ! curl -L "$url" > "$output_file" 2>> "$log"; then
        echo "无法从Ensembl下载注释数据。" >&2
        echo "请检查物种、构建版本和发布版本的组合是否实际提供。" >&2
        exit 1
    fi
else
    if ! (curl -L "$url" | gzip -d > "$output_file") 2>> "$log"; then
        echo "无法从Ensembl下载注释数据。" >&2
        echo "请检查物种、构建版本和发布版本的组合是否实际提供。" >&2
        exit 1
    fi
fi