/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef MAIN_VCF_INDEX_H
#define MAIN_VCF_INDEX_H

#include <stddef.h>
#include <stdint.h>

/* fast-beagle only: the tabix index of a BGZF VCF, byte-identical to
 * tabix -p vcf. Records are collected with their uncompressed end offsets as
 * the VCF is written; virtual offsets are read from the closed file, as a
 * multithreaded BGZF writer cannot report them. */
typedef struct vcf_index vcf_index;

/* The index of vcf_path, written to tbi_path, whose header lines are
 * u_header_end uncompressed bytes. */
vcf_index *vcf_index_new(const char *vcf_path, const char *tbi_path, uint64_t u_header_end);
/* A record whose line ends at uncompressed offset u_end, just past its
 * newline. line[0, info_end) runs through the INFO column and line[info_end]
 * is a tab; the line is changed and restored. */
void vcf_index_add(vcf_index *x, char *line, size_t info_end, uint64_t u_end);
/* After the BGZF file of u_total uncompressed bytes is closed: writes
 * tbi_path and frees x. */
void vcf_index_write(vcf_index *x, uint64_t u_total);

#endif
