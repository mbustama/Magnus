#!/bin/bash
S=/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad
LOG=$S/nb28_par.log
run() {  # $1 copy dir, $2 extra env
  cd $S/par_$1/notebooks && env -u MAGNUS_PAPER_CACHE_ONLY MAGNUS_PAPER_REDO=1 MAGNUS_PAPER_REDO_ONLY="$(cat redo_only.txt)" \
    MAGNUS_PAPER_FIGDIR=$S/par_$1/notebooks/figs $2 jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=-1 28_magnus_paper_figures.ipynb --output 28_redo.ipynb > $S/$1.exec.log 2>&1
  rc=$?; [ $rc -eq 0 ] || touch $S/FAILED_$1
  echo "$(date +%T) $1 exit $rc; mem avail $(free -m | awk '/Mem/ {print $7}') MB" >> $LOG
}
echo "$(date +%T) parallel pass: four copies, value-only sections, no timing; mem avail $(free -m | awk '/Mem/ {print $7}') MB" > $LOG
for w in 0 1 2 3; do run nb28_w$w MAGNUS_PAPER_NO_TIMING=1 & done
wait
if ls $S/FAILED_* > /dev/null 2>&1; then echo "$(date +%T) FAILED: $(ls $S/FAILED_* | xargs -n1 basename); stopping before the solo pass" >> $LOG; exit 2; fi
echo "$(date +%T) parallel pass done; waiting for the usual quiet level (1-min and 5-min load <= 1.0, three minutes)" >> $LOG
ok=0; deadline=$(( $(date +%s) + 3*3600 ))
while [ $ok -lt 3 ]; do
  read l1 l5 rest < /proc/loadavg
  if awk "BEGIN{exit !($l1 <= 1.0 && $l5 <= 1.0)}"; then ok=$((ok+1)); else ok=0; fi
  if [ $(date +%s) -gt $deadline ]; then echo "$(date +%T) NOT QUIET within 3 h; solo pass not run" >> $LOG; exit 3; fi
  [ $ok -lt 3 ] && sleep 60
done
echo "$(date +%T) solo pass: the timed prob_vs sections, one process" >> $LOG
run nb28_solo ""
if [ -e $S/FAILED_nb28_solo ]; then echo "$(date +%T) solo pass FAILED; no merge" >> $LOG; exit 2; fi
python $S/merge_nb28.py >> $LOG 2>&1
python $S/diff_cache.py /home/mbustamante/Research/magnus/notebooks/paper_figure_cache.json $S/paper_figure_cache_recomputed.json > $S/nb28_cache_diff.txt 2>&1
echo "$(date +%T) diff written; TIMING-SKIPPED markers: $(grep -l TIMING-SKIPPED $S/par_nb28_w*/notebooks/28_redo.ipynb 2>/dev/null | wc -l)" >> $LOG
