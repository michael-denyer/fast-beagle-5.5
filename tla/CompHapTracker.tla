---------------------------- MODULE CompHapTracker ----------------------------
(***************************************************************************)
(* The composite-haplotype tracker of src/beagleutil/comp_hap_queue.c:     *)
(* comp_hap_tracker_observe with update_head_of_q and add_segment,         *)
(* comp_hap_tracker_seed, comp_hap_tracker_clear, and the binary heap      *)
(* comp_hap_queue_offer and comp_hap_queue_poll it keeps its segments in.  *)
(* The heap is the array, so the order in which equal keys leave is the    *)
(* C's. Haplotypes are 1 .. H and steps 0 .. MaxStep. The environment is   *)
(* the caller, as in basic_phase_states.c, low_freq_phase_states.c and     *)
(* imp_ls.c: observations at non-decreasing steps, seeding only into an    *)
(* empty tracker and then no observation until clear, and clear before     *)
(* each new target, up to MaxEvents calls. MinSteps > MaxStep stands for   *)
(* INT_MAX, which ImpStates uses. Each step is one pass of the             *)
(* update_head_of_q loop, the retire-or-add decision, the map update, or   *)
(* one seed or clear.                                                      *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets

CONSTANTS H, MaxStates, MaxStep, MinSteps, MaxEvents

ASSUME H \in Nat \ {0} /\ MaxStates \in Nat \ {0} /\ MaxStep \in Nat /\ MinSteps \in Nat \ {0}

VARIABLES
    segs,      \* t->segs[0 .. ), as records; index i is comp_hap_index i - 1
    heap,      \* t->q.es[0 .. size), as indices into segs
    map,       \* t->last_ibs_step
    lastStep,  \* the caller's step so far
    phase,     \* "obs", or "seeded" from the first seed to the next clear
    pc,        \* "idle", "repair" (update_head_of_q), "decide" or "set"
    cur,       \* the last observation, and whether its hap was tracked; a seed or
               \* clear since then resets tracked
    before,    \* the slots when it started
    change,    \* what it returned
    evictOk,   \* every retired segment had the least true step in the queue
    events     \* caller calls so far

vars == <<segs, heap, map, lastStep, phase, pc, cur, before, change, evictOk, events>>

NIL == -103
Haps == 1 .. H
NoChange == [index |-> -1, old |-> -1, start |-> 0, end |-> 0]
Seg == [hap : Haps, start : 0 .. MaxStep, key : {NIL} \cup 0 .. MaxStep]

Key(sg) == [i \in 1 .. Len(sg) |-> sg[i].key]

(* comp_hap_queue_offer: sift up, stopping at a parent with an equal key. *)
RECURSIVE Up(_, _, _, _)
Up(h, p, x, k) ==
    IF p > 1 /\ k[x] < k[h[p \div 2]]
    THEN Up([h EXCEPT ![p] = h[p \div 2]], p \div 2, x, k)
    ELSE [h EXCEPT ![p] = x]
Offer(h, x, k) == Up(Append(h, x), Len(h) + 1, x, k)

(* comp_hap_queue_poll: the last element sifts down from the root, taking
   the left child on equal children and stopping at an equal child. *)
RECURSIVE Down(_, _, _, _, _)
Down(h, p, x, n, k) ==
    IF 2 * p > n THEN [h EXCEPT ![p] = x]
    ELSE LET c == IF 2 * p + 1 <= n /\ k[h[2 * p + 1]] < k[h[2 * p]] THEN 2 * p + 1 ELSE 2 * p
         IN  IF k[x] <= k[h[c]] THEN [h EXCEPT ![p] = x]
             ELSE Down([h EXCEPT ![p] = h[c]], c, x, n, k)
Poll(h, k) ==
    LET n == Len(h) - 1
    IN  IF n = 0 THEN <<>> ELSE Down(SubSeq(h, 1, n), 1, h[n + 1], n, k)

(* add_segment's index, q.size *)
NextIndex == Len(heap) + 1

AddSeg(sg, hp, step) ==
    LET i == NextIndex
        s == [hap |-> hp, start |-> 0, key |-> step]
    IN  IF i <= Len(sg) THEN [sg EXCEPT ![i] = s] ELSE Append(sg, s)

Queued == {heap[p] : p \in 1 .. Len(heap)}
TrueStep(i) == map[segs[i].hap]
Min(X) == CHOOSE x \in X : \A y \in X : x <= y
Slots == [i \in 1 .. Len(heap) |-> <<segs[i].hap, segs[i].start>>]

Untracked(h, st) == map[h] = NIL

TypeOK ==
    /\ segs \in Seq(Seg) /\ Len(segs) <= MaxStates
    /\ heap \in Seq(1 .. MaxStates)
    /\ map \in [Haps -> {NIL} \cup 0 .. MaxStep]
    /\ lastStep \in 0 .. MaxStep
    /\ phase \in {"obs", "seeded"}
    /\ pc \in {"idle", "repair", "decide", "set"}
    /\ evictOk \in BOOLEAN
    /\ events \in 0 .. MaxEvents

Init ==
    /\ segs = <<>> /\ heap = <<>>
    /\ map = [h \in Haps |-> NIL]
    /\ lastStep = 0 /\ phase = "obs" /\ pc = "idle"
    /\ cur = [hap |-> 1, step |-> 0, tracked |-> FALSE]
    /\ before = <<>> /\ change = NoChange
    /\ evictOk = TRUE /\ events = 0

(* comp_hap_tracker_observe(t, h, st): the repair and the decision run only
   for a haplotype the map does not hold. *)
Observe(h, st) ==
    /\ pc = "idle" /\ phase = "obs" /\ events < MaxEvents
    /\ st >= lastStep
    /\ cur' = [hap |-> h, step |-> st, tracked |-> map[h] # NIL]
    /\ before' = Slots /\ change' = NoChange
    /\ lastStep' = st /\ events' = events + 1
    /\ pc' = IF Untracked(h, st) THEN "repair" ELSE "set"
    /\ UNCHANGED <<segs, heap, map, phase, evictOk>>

(* One pass of update_head_of_q: a head whose key is behind its haplotype's
   step leaves, takes that step, and goes back in. *)
Repair ==
    /\ pc = "repair"
    /\ IF heap = <<>> \/ segs[heap[1]].key = TrueStep(heap[1])
         THEN pc' = "decide" /\ UNCHANGED <<segs, heap>>
         ELSE LET hd == heap[1]
                  sg == [segs EXCEPT ![hd].key = TrueStep(hd)]
              IN  /\ segs' = sg
                  /\ heap' = Offer(Poll(heap, Key(segs)), hd, Key(sg))
                  /\ UNCHANGED pc
    /\ UNCHANGED <<map, lastStep, phase, cur, before, change, evictOk, events>>

(* The queue is full, or its head is stale: the head is retired and its
   segment takes the haplotype from the midpoint step. Otherwise a new
   segment is added. *)
Decide ==
    /\ pc = "decide"
    /\ LET full == Len(heap) = MaxStates
           stale == heap # <<>> /\ cur.step - segs[heap[1]].key >= MinSteps
       IN  IF full \/ stale
           THEN LET hd == heap[1]
                    old == segs[hd]
                    mid == (old.key + cur.step) \div 2
                    sg == [segs EXCEPT ![hd] = [hap |-> cur.hap, start |-> mid, key |-> cur.step]]
                IN  /\ evictOk' = (evictOk /\ TrueStep(hd) = Min({TrueStep(i) : i \in Queued}))
                    /\ map' = [map EXCEPT ![old.hap] = NIL]
                    /\ segs' = sg
                    /\ heap' = Offer(Poll(heap, Key(segs)), hd, Key(sg))
                    /\ change' = [index |-> hd - 1, old |-> old.hap, start |-> old.start, end |-> mid]
           ELSE LET sg == AddSeg(segs, cur.hap, cur.step)
                IN  /\ segs' = sg
                    /\ heap' = Offer(heap, NextIndex, Key(sg))
                    /\ change' = [NoChange EXCEPT !.index = NextIndex - 1]
                    /\ UNCHANGED <<map, evictOk>>
    /\ pc' = "set"
    /\ UNCHANGED <<lastStep, phase, cur, before, events>>

(* t->last_ibs_step[hap] = step *)
Set ==
    /\ pc = "set"
    /\ map' = [map EXCEPT ![cur.hap] = cur.step]
    /\ pc' = "idle"
    /\ UNCHANGED <<segs, heap, lastStep, phase, cur, before, change, evictOk, events>>

(* comp_hap_tracker_seed, from fill_q_with_random_haps after observations
   left the tracker empty; draws may repeat. *)
Seed(h) ==
    /\ pc = "idle" /\ events < MaxEvents
    /\ Len(heap) < MaxStates
    /\ phase = "seeded" \/ heap = <<>>
    /\ LET sg == AddSeg(segs, h, 0)
       IN  segs' = sg /\ heap' = Offer(heap, NextIndex, Key(sg))
    /\ phase' = "seeded" /\ events' = events + 1
    /\ cur' = [cur EXCEPT !.tracked = FALSE]
    /\ UNCHANGED <<map, lastStep, pc, before, change, evictOk>>

(* comp_hap_tracker_clear: NIL for the haplotype of each segment in use. *)
Clear ==
    /\ pc = "idle" /\ events < MaxEvents
    /\ map' = [h \in Haps |-> IF \E i \in 1 .. Len(heap) : segs[i].hap = h THEN NIL ELSE map[h]]
    /\ heap' = <<>> /\ phase' = "obs" /\ lastStep' = 0 /\ events' = events + 1
    /\ cur' = [cur EXCEPT !.tracked = FALSE]
    /\ UNCHANGED <<segs, pc, before, change, evictOk>>

Done == pc = "idle" /\ events = MaxEvents

Stop == Done /\ UNCHANGED vars

Next ==
    \/ \E h \in Haps, st \in 0 .. MaxStep : Observe(h, st)
    \/ \E h \in Haps : Seed(h)
    \/ Repair \/ Decide \/ Set \/ Clear
    \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Repair \/ Decide \/ Set)

-----------------------------------------------------------------------------
(* At most max_states segments, which are segs[0 .. size), each queued once. *)
SizeBound == Len(heap) <= MaxStates
SlotsArePrefix == Queued = 1 .. Len(heap) /\ Cardinality(Queued) = Len(heap)

(* The heap order on the stored keys. *)
HeapOrdered == \A p \in 2 .. Len(heap) : segs[heap[p \div 2]].key <= segs[heap[p]].key

(* Outside seeding, the map holds a step exactly for the queued haplotypes,
   each queued once. *)
NilIffSlot ==
    (phase = "obs" /\ pc = "idle") =>
        /\ \A h \in Haps : map[h] # NIL <=> \E i \in Queued : segs[i].hap = h
        /\ \A i, j \in Queued : i # j => segs[i].hap # segs[j].hap

(* A stored key never runs ahead of its haplotype's last IBS step. *)
KeyBehindStep ==
    (phase = "obs" /\ pc # "set") => \A i \in Queued : segs[i].key <= TrueStep(i)

(* The lazily repaired head is a segment with the least true step. *)
EvictsLeast == evictOk

(* Observing a tracked haplotype returns no change and moves no segment. *)
TrackedNoChange == (pc = "idle" /\ cur.tracked) => change = NoChange /\ Slots = before

(* update_head_of_q returns. *)
RepairEnds == pc = "repair" ~> pc = "idle"
=============================================================================
