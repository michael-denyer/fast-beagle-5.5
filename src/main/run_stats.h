/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) main/RunStats.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef MAIN_RUN_STATS_H
#define MAIN_RUN_STATS_H

#include <stdint.h>
#include <stdio.h>

#include "main/run_outputs.h"
#include "phase/fixed_phase_data.h"
#include "phase/phase_data.h"
#include "vcf/sliding_window.h"

/* The run's progress report, printed to standard output and to <out>.log in
 * Beagle's format. The banner and command line name this program instead of
 * the jar, and the cumulative statistics add the CPU time and max memory. */
typedef struct {
    const par *par;
    FILE *log;
    const char *log_path;   /* borrowed from the run outputs */
    int64_t start_nanos;
    int64_t total_phase_nanos;
    int64_t total_impute_nanos;
} run_stats;

/* A monotonic clock reading for timing the steps reported here. */
int64_t run_stats_nanos(void);

/* new RunStats(par) and printStartInfo(): opens the run's log and prints the
 * banner, start time and command line. From here until run_stats_close, an
 * error message also goes to the log. */
void run_stats_open(run_stats *rs, const par *p, const run_outputs *out, const char *program);
/* printSampleSummary */
void run_stats_sample_summary(run_stats *rs, int n_ref_samples, int n_targ_samples);
/* printWindowUpdate */
void run_stats_window_update(run_stats *rs, const window *w, const fixed_phase_data *fpd);
/* printStage1Info: call after the iteration pd->it has run, before it advances. */
void run_stats_stage1(run_stats *rs, const phase_data *pd, int64_t nanos);
/* printStage2Info */
void run_stats_stage2(run_stats *rs, int64_t nanos);
/* imputationNanos and printImputationUpdate */
void run_stats_imputation(run_stats *rs, int64_t nanos);
/* printSummaryAndClose */
void run_stats_close(run_stats *rs, int64_t n_targ_markers, int64_t n_markers);

#endif
