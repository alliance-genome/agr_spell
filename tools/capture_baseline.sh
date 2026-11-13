#!/usr/bin/env bash
# Send validated queries one at a time to the live spell++ engine on sgd-spell and
# record its two response lines. Usage: capture_baseline.sh QUERY_IDS_TSV [QID ...]
set -euo pipefail
cd "$(dirname "$0")/.."
QFILE=$1; shift
ids=$(cut -f1 data/raw/dataset_list.txt | paste -sd, -)
valid=$(cut -f1 data/raw/gene_list.txt | sort -u)
{
  printf 'DS=%q\n' "$ids"
  echo 'QS=()'
  while IFS=$'\t' read -r qid genes; do
    if [ $# -gt 0 ] && ! printf '%s\n' "$@" | grep -qx "$qid"; then continue; fi
    # Safety: non-empty, numeric, every id known to the engine.
    [[ $genes =~ ^[0-9]+(,[0-9]+)*$ ]] || { echo "bad query $qid: '$genes'" >&2; exit 1; }
    for g in ${genes//,/ }; do grep -qx "$g" <<<"$valid" || { echo "unknown gene $g in $qid" >&2; exit 1; }; done
    printf 'QS+=(%q)\n' "$qid:$genes"
  done < "$QFILE"
  cat <<'REMOTE'
for item in "${QS[@]}"; do
  qid=${item%%:*}; genes=${item#*:}
  exec 3<>/dev/tcp/127.0.0.1/5011 || { echo "$qid	ERROR	connect"; exit 1; }
  printf '%s\t%s\n' "$DS" "$genes" >&3
  l1=; l2=
  IFS= read -r -t 120 l1 <&3 || true
  IFS= read -r -t 120 l2 <&3 || true
  exec 3<&-
  printf '%s\t%s\t%s\n' "$qid" "$l1" "$l2"
  [ -n "$l2" ] || { echo "$qid	ERROR	empty-response" ; exit 1; }
  sleep 1
done
REMOTE
} | ssh -o BatchMode=yes sgd-spell 'bash -s' 2>/dev/null
