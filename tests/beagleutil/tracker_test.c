/* Composite haplotype transitions, through the interface used by all callers. */
#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "beagleutil/comp_hap_queue.h"

static void change(comp_hap_change c, int index, int old_hap, int start, int end) {
    assert(c.index == index);
    assert(c.old_hap == old_hap);
    assert(c.start_step == start);
    assert(c.end_step == end);
}

int main(void) {
    comp_hap_tracker t;
    comp_hap_tracker_init(&t, 3, 61);
    change(comp_hap_tracker_observe(&t, 10, 0, INT_MAX), 0, -1, 0, 0);
    change(comp_hap_tracker_observe(&t, 20, 0, INT_MAX), 1, -1, 0, 0);
    change(comp_hap_tracker_observe(&t, 30, 0, INT_MAX), 2, -1, 0, 0);
    change(comp_hap_tracker_observe(&t, 10, 4, INT_MAX), -1, -1, 0, 0);
    /* Refreshing the old head exposes 30 before 20 under Java's tie rules. */
    change(comp_hap_tracker_observe(&t, 40, 4, INT_MAX), 2, 30, 0, 2);
    change(comp_hap_tracker_observe(&t, 30, 5, INT_MAX), 1, 20, 0, 2);
    change(comp_hap_tracker_observe(&t, 40, 6, INT_MAX), -1, -1, 0, 0);
    change(comp_hap_tracker_observe(&t, 50, 8, INT_MAX), 0, 10, 0, 6);
    change(comp_hap_tracker_observe(&t, 60, 9, INT_MAX), 1, 30, 2, 7);
    assert(comp_hap_tracker_size(&t) == 3);
    const comp_hap_segment *s = comp_hap_tracker_segment(&t, 1);
    assert(s->hap == 60 && s->start_step == 7);
    puts("PASS lazy refresh, Java ties, eviction and retired step ranges");

    comp_hap_tracker_clear(&t);
    change(comp_hap_tracker_observe(&t, 10, 0, 3), 0, -1, 0, 0);
    change(comp_hap_tracker_observe(&t, 20, 2, 3), 1, -1, 0, 0);
    change(comp_hap_tracker_observe(&t, 30, 3, 3), 0, 10, 0, 1);
    assert(comp_hap_tracker_size(&t) == 2);
    change(comp_hap_tracker_observe(&t, 10, 4, 3), 2, -1, 0, 0);
    puts("PASS stale retirement before capacity, and removed haplotype reentry");

    comp_hap_tracker_clear(&t);
    comp_hap_tracker_observe(&t, 10, 0, INT_MAX);
    change(comp_hap_tracker_observe(&t, 20, 100, INT_MAX), 1, -1, 0, 0);
    comp_hap_tracker_clear(&t);
    comp_hap_tracker_observe(&t, 10, 0, 3);
    comp_hap_tracker_observe(&t, 10, 0, 3);
    assert(comp_hap_tracker_size(&t) == 1);
    comp_hap_tracker_clear(&t);
    assert(comp_hap_tracker_seed(&t, 10) == 0);
    assert(comp_hap_tracker_seed(&t, 10) == 1);
    assert(comp_hap_tracker_size(&t) == 2);
    assert(comp_hap_tracker_segment(&t, 0)->hap == 10);
    assert(comp_hap_tracker_segment(&t, 1)->hap == 10);
    comp_hap_tracker_free(&t);
    puts("PASS disabled staleness, clear, distinct and repeated fallback draws");

    /* An offered key equal to its parent's stays below it. */
    comp_hap_tracker_init(&t, 2, 6);
    comp_hap_tracker_observe(&t, 0, 2, INT_MAX);
    comp_hap_tracker_observe(&t, 3, 2, INT_MAX);
    change(comp_hap_tracker_observe(&t, 2, 2, INT_MAX), 0, 0, 0, 2);
    comp_hap_tracker_free(&t);
    /* A polled heap sifts its last node down past an equal child only when
     * that child is the smaller, taking the left of two equal children. */
    comp_hap_tracker_init(&t, 4, 6);
    comp_hap_tracker_observe(&t, 1, 0, INT_MAX);
    comp_hap_tracker_observe(&t, 2, 1, INT_MAX);
    comp_hap_tracker_observe(&t, 4, 1, INT_MAX);
    comp_hap_tracker_observe(&t, 5, 2, INT_MAX);
    comp_hap_tracker_observe(&t, 1, 2, INT_MAX);
    change(comp_hap_tracker_observe(&t, 0, 3, INT_MAX), 1, 2, 0, 2);
    comp_hap_tracker_free(&t);
    puts("PASS java.util.PriorityQueue order for equal keys on offer and poll");

    comp_hap_tracker_init(&t, 1, 21);
    comp_hap_tracker_observe(&t, 10, 0, 3);
    change(comp_hap_tracker_observe(&t, 20, 0, 3), 0, 10, 0, 0);
    comp_hap_tracker_free(&t);
    puts("PASS zero-length retirement at one step");
    return 0;
}
