#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/cheney/Projects/PeleC"

dirs=(
  "Source"
  "Submodules/PelePhysics/Source"
  "Submodules/PelePhysics/Submodules/amrex"
  "Submodules/PelePhysics/Submodules/sundials"
)

# 需要统计的扩展名
exts=(cpp c C H h hpp cuh cu cc)

grand_files=0
grand_lines=0

for d in "${dirs[@]}"; do
  dir_path="$ROOT/$d"
  if [[ ! -d "$dir_path" ]]; then
    echo "Directory not found: $dir_path"
    continue
  fi

  # 构造 find 条件
  find_expr=()
  for ext in "${exts[@]}"; do
    find_expr+=(-name "*.${ext}")
    find_expr+=(-o)
  done
  # 去掉最后一个 -o
  unset 'find_expr[${#find_expr[@]}-1]'

  # 查找文件
  mapfile -t files < <(find "$dir_path" -type f \( "${find_expr[@]}" \))

  file_count=${#files[@]}
  if [[ $file_count -eq 0 ]]; then
    total_lines=0
  else
    # 统计总行数
    total_lines=$(wc -l "${files[@]}" | awk 'END{print $1}')
  fi

  echo "Directory: $d"
  echo "  Files: $file_count"
  echo "  Lines: $total_lines"

  grand_files=$((grand_files + file_count))
  grand_lines=$((grand_lines + total_lines))
done

echo "Total"
echo "  Files: $grand_files"
echo "  Lines: $grand_lines"
