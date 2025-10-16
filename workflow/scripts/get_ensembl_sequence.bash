#!/usr/bin/env bash

# 获取Snakemake参数
species=$(echo "${snakemake_params[species]}" | tr '[:upper:]' '[:lower:]')
release=${snakemake_params[release]}
build=${snakemake_params[build]}
datatype=${snakemake_params[datatype]:-""}
threads=${snakemake[threads]}
output_file=${snakemake_output[0]}
log=${snakemake_log[0]}

# 设置branch
branch=""
if [ "$release" -ge 81 ] && [ "$build" = "GRCh37" ]; then
    branch="grch37/"
elif [ -n "${snakemake_params[branch]}" ]; then
    branch="${snakemake_params[branch]}/"
fi

# 设置spec
if [ "$release" -gt 75 ]; then
    spec="$build"
else
    spec="${build}.${release}"
fi

# 检查是否需要解压
if [[ "$output_file" == *.gz ]]; then
    decompress=""
else
    decompress="| gzip -d"
fi

# 设置URL基础部分
url="https://ftp.ensembl.org/pub"
url_prefix="${url}/${branch}release-${release}/fasta/${species}/${datatype}/${species^}.${spec}"


tmpdir=$(mktemp -d -p ./)
trap 'rm -rf "$tmpdir"' EXIT


# 根据数据类型设置后缀
declare -a suffixes
if [ "$datatype" = "dna" ]; then
    suffixes=("dna.primary_assembly.fa.gz" "dna.toplevel.fa.gz")
elif [ "$datatype" = "cdna" ]; then
    suffixes=("cdna.all.fa.gz")
elif [ "$datatype" = "cds" ]; then
    suffixes=("cds.all.fa.gz")
elif [ "$datatype" = "ncrna" ]; then
    suffixes=("ncrna.fa.gz")
elif [ "$datatype" = "pep" ]; then
    suffixes=("pep.all.fa.gz")
else
    echo "错误：无效的数据类型，必须是 dna、cdna、cds、ncrna 或 pep 之一" >&2
    exit 1
fi

success=false
for suffix in "${suffixes[@]}"; do
    url_https="${url_prefix}.${suffix}"
    url_ftp="${url_https/https:\/\//ftp:\/\/}"
    
    # 检查文件是否存在
    if curl --location --head "$url_https" 2>/dev/null | grep -q 'Content-Length'; then
        (lftp -c "pget -n ${threads} ${url_https} -" $decompress >> "$output_file") 2>&1 | tee -a "$log"
        success=true

    elif curl --location --head "$url_ftp" 2>/dev/null | grep -q 'Content-Length'; then
        (lftp -c "pget -n ${threads} ${url_ftp} -" $decompress >> "$output_file") 2>&1 | tee -a "$log"
        success=true
    fi
done

if ! $success; then

    echo "请检查以上URL是否当前可用（可能是临时服务器问题）。" >&2
    echo "此外，请检查物种、构建版本和发布版本的组合是否实际提供。" >&2
    exit 1
fi