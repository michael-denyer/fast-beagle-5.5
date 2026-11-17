/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef MAIN_RUN_OUTPUTS_H
#define MAIN_RUN_OUTPUTS_H

#include <stdio.h>

#include "main/par.h"

/* fast-beagle's output destinations. Paths are immutable and borrowed by
 * their format writers. */
typedef struct run_outputs run_outputs;
typedef enum {
    RUN_OUTPUT_VCF, RUN_OUTPUT_LOG, RUN_OUTPUT_BGEN, RUN_OUTPUT_INFO,
    RUN_OUTPUT_SAMPLE, RUN_OUTPUT_TBI, RUN_OUTPUT_COUNT
} run_output_file;

/* Check destinations against inputs and each other during parameter
 * validation, before opening any writer. */
void run_outputs_check(const par *p);
run_outputs *run_outputs_new(const par *p);
/* Paths include an existing tabix index scheduled for invalidation. */
const char *run_outputs_path(const run_outputs *out, run_output_file file);
/* Open and check the actual destination before truncation. VCF open also
 * invalidates the old index. Stream opens return NULL on error and set
 * conflict to the borrowed path that aliases this output, if any. */
int run_outputs_open_vcf(const run_outputs *out);
FILE *run_outputs_open_stream(const run_outputs *out, run_output_file file, const char **conflict);
FILE *run_outputs_open_trace(const run_outputs *out, const char *path, const char **conflict);
void run_outputs_free(run_outputs *out);

#endif
