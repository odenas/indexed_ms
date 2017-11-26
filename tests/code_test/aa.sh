#!/bin/env bash
EXE=/home/brt/Documents/projects/matching_statistics/indexed_ms/fast_ms/bin/matching_stats.opt.x
SPATH=datasets/testing/rnd_20_10.s
TPATH=datasets/testing/rnd_20_10.t

# lazy/single_rank doesn't work with maxrep

$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0 2>/dev/null
echo $?
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0 2>/dev/null
echo $?
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 0 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0 2>/dev/null
echo $?
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 0 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0 2>/dev/null
echo $?
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 0 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 0
$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 0 -double_rank 1 2>/dev/null
echo $?
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 0 -use_maxrep_vanilla 1 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 0 -double_rank 1
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 0
#$EXE -s_path $SPATH -t_path $TPATH -lca_parents 1 -lazy_wl 1 -rank_fail 1 -use_maxrep_rc 1 -use_maxrep_vanilla 1 -double_rank 1


