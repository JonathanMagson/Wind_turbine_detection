#!/usr/bin/env bash
set -e

venv=$1
filename=$2
metaname=$3
t_NFA=$4
t_shadow=$5
t_hub=$6


if [ ! -d "$venv" ]; then
  echo "Virtualenv not found" > demo_failure.txt
  exit 0
else
  source $venv/bin/activate
fi

main_IPOL_OH_nogdal.py --filename $filename --metaname $metaname --t_NFA $t_NFA --t_shadow $t_shadow --t_hub $t_hub
