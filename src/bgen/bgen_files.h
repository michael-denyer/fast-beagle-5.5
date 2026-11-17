/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef BGEN_BGEN_FILES_H
#define BGEN_BGEN_FILES_H

#include "main/run_outputs.h"

/* Removal of partial BGEN output at process exit. A run creates at most the
 * three members (.bgen, .info, .sample); exit removes the ones it created
 * unless bgen_files_complete ran first. Any thread's exit() can run the
 * removal, so creation and removal are serialized, and removal is terminal:
 * later opens fail. */

/* Once, before starting any thread that can exit. */
void bgen_files_register_cleanup(void);
/* Open a checked destination for writing and borrow its path until exit or
 * completion. An alias is fatal. An open error or a call after cleanup or
 * completion returns NULL. */
FILE *bgen_files_open(const run_outputs *out, run_output_file file);
/* After every member is written and closed: exit keeps them. */
void bgen_files_complete(void);

#endif
