/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef BLBUTIL_TRACE_H
#define BLBUTIL_TRACE_H

#include <stdbool.h>
#include <stdint.h>

struct run_outputs;

/* Trace seams for comparison with the Java trace build (java/trace.patch).
 * Enabled by trace=<dir>; each seam writes <dir>/<seam>.txt. */
void trace_init(const char *dir, const struct run_outputs *outputs);
bool trace_on(void);
void trace_line(const char *seam, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
void trace_close(void);

/* FNV-1a over 32-bit words, for the digests in trace lines. */
#define TRACE_FNV_BASIS 14695981039346656037ULL
static inline uint64_t trace_fold(uint64_t h, uint32_t v) {
    return (h ^ v) * 1099511628211ULL;
}

#endif
