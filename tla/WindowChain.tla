------------------------------ MODULE WindowChain ------------------------------
(***************************************************************************)
(* The chain of windows the reader thread builds in                        *)
(* src/vcf/sliding_window.c: read_windows, next_targ_window (Java's        *)
(* TargSlidingWindow.Reader), next_ref_window                              *)
(* (RefTargSlidingWindow.Reader), overlap_start and check_advances. The    *)
(* hand-over of each window to the caller is tla/SlidingWindow.tla; here a *)
(* window is one step. Mode "targ" has no reference; "ref" and "impute"    *)
(* read the reference with impute=false and impute=true. Init chooses the  *)
(* target records T and reference records R, each a marker [c              *)
(* (chromosome), p (position), a (allele set)], sorted by chromosome and   *)
(* position. In "targ" mode records at one position stand for markers with *)
(* different alleles, so a is 1 and only c and p are used. The genetic map *)
(* is the identity, one cM per base, so window=W and overlap=O are in      *)
(* bases; this hides the rounding of genetic_map_base_pos and flat         *)
(* stretches of a real map. WM is window-markers and OM is                 *)
(* overlap_markers. A window that check_advances rejects is not made (pc = *)
(* "stalled"); empty_window_error and Java's empty-window exception end    *)
(* the chain at pc = "empty".                                              *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets

CONSTANTS Mode, NT, NR, NC, P, A, W, O, WM

ASSUME /\ Mode \in {"targ", "ref", "impute"}
       /\ NT \in Nat /\ NR \in Nat /\ (Mode = "targ" => NR = 0 /\ A = 1)
       /\ {NC, P, A, WM} \subseteq Nat \ {0}
       /\ O \in Nat /\ W > O /\ 10 * W >= 11 * O

OM == WM \div 4

Markers == [c : 1..NC, p : 1..P, a : 1..A]

Before(m, n) == m.c < n.c \/ (m.c = n.c /\ m.p <= n.p)

(* Input files sorted by chromosome and position; a reference or target VCF
   holds a marker once. *)
Inputs(n) == {s \in [1..n -> Markers] :
                /\ \A j \in 1..n - 1 : Before(s[j], s[j + 1])
                /\ Mode # "targ" => \A j, k \in 1..n : j # k => s[j] # s[k]}

Entries == [r : 1..NR, in : BOOLEAN]

VARIABLES
    T, R,     \* the target and reference records
    pc,       \* "start", "stop", "stalled", "empty" or "nodata"
    next,     \* the index of sw->next_targ; NT + 1 for NULL
    rnext,    \* the index of sw->next_ref; NR + 1 for NULL
    nRead,    \* sw->n_read
    tov,      \* sw->targ_overlap, as indices into T
    rov,      \* sw->ref_overlap, as indices into R with in_targ
    twin,     \* the target records of the last window made
    rwin,     \* the reference records of the last window made
    inLen,    \* how many records of the last window were carried into it
    seen,     \* the target records of every window made
    done,     \* sw->done
    fixed     \* the last window, made or rejected, left the reader where it found it

vars == <<T, R, pc, next, rnext, nRead, tov, rov, twin, rwin, inLen, seen, done, fixed>>

Range(s) == {s[j] : j \in 1..Len(s)}

TypeOK ==
    /\ T \in Inputs(NT) /\ R \in Inputs(NR)
    /\ pc \in {"start", "stop", "stalled", "empty", "nodata"}
    /\ next \in 1..NT + 1 /\ rnext \in 1..NR + 1
    /\ nRead \in 0..NT + NR + 2
    /\ tov \in Seq(1..NT) /\ rov \in Seq(Entries)
    /\ twin \in Seq(1..NT) /\ rwin \in Seq(Entries)
    /\ inLen \in Nat /\ seen \subseteq 1..NT
    /\ done \in BOOLEAN /\ fixed \in BOOLEAN

(* read_windows reads the first record of each input. *)
Init ==
    /\ T \in Inputs(NT) /\ R \in Inputs(NR)
    /\ next = 1 /\ rnext = 1
    /\ nRead = IF Mode = "targ" THEN 1 ELSE 2
    /\ pc = IF NT = 0 \/ (Mode # "targ" /\ NR = 0) THEN "nodata" ELSE "start"
    /\ tov = <<>> /\ rov = <<>> /\ twin = <<>> /\ rwin = <<>>
    /\ inLen = 0 /\ seen = {} /\ done = FALSE /\ fixed = FALSE

(* overlap_start over the positions ps of a window's markers: binary search
   of the last overlap_markers markers for the position overlap= before
   endGen, then back to the first marker at that position. Indices are
   0-based, as in the C. *)
RECURSIVE FirstWithPos(_, _)
FirstWithPos(ps, j) == IF j > 0 /\ ps[j] = ps[j + 1] THEN FirstWithPos(ps, j - 1) ELSE j

RECURSIVE Search(_, _, _, _)
Search(ps, key, low, high) ==
    IF low > high THEN FirstWithPos(ps, IF high > 0 THEN high ELSE 0)
    ELSE LET mid == (low + high) \div 2
         IN  IF ps[mid + 1] < key THEN Search(ps, key, mid + 1, high)
             ELSE IF ps[mid + 1] > key THEN Search(ps, key, low, mid - 1)
             ELSE FirstWithPos(ps, mid)

OverlapStart(ps, endGen) ==
    Search(ps, endGen - O, IF Len(ps) - OM > 0 THEN Len(ps) - OM ELSE 0, Len(ps) - 1)

(* The read loop of next_targ_window. *)
RECURSIVE TLoop(_, _, _, _)
TLoop(recs, k, chrom, endPos) ==
    IF k <= NT /\ T[k].c = chrom /\ T[k].p < endPos /\ Len(recs) < WM
    THEN TLoop(Append(recs, k), k + 1, chrom, endPos)
    ELSE <<recs, k>>

(* next_targ_window, then check_advances and the carried overlap. The
   stages are bound with \E x \in {e}, not LET, because TLC evaluates a LET
   again each time it checks the fairness condition; with LET the reference
   modes ran ten times slower. *)
TargWindow ==
    /\ pc = "start" /\ ~done /\ Mode = "targ"
    /\ LET chrom == T[next].c
           endPos == T[next].p + (IF tov = <<>> THEN W ELSE W - O)
       IN  \E lp \in {TLoop(tov, next, chrom, endPos)} :
           LET recs == lp[1]
               k == lp[2]
               nr == nRead + (k - next)
               last == k > NT
               chromEnd == last \/ T[k].c # chrom
           IN  IF recs = <<>>
               THEN pc' = "empty" /\ UNCHANGED <<next, nRead, tov, twin, inLen, seen, done, fixed>>
               ELSE \E ovStart \in {IF chromEnd THEN Len(recs)
                                     ELSE OverlapStart([j \in 1..Len(recs) |-> T[recs[j]].p], T[recs[Len(recs)]].p)} :
                    LET carry == SubSeq(recs, ovStart + 1, Len(recs))
                    IN  /\ fixed' = (k = next /\ carry = tov)
                        /\ IF ~last /\ ovStart = 0 /\ nr = nRead
                           THEN pc' = "stalled" /\ UNCHANGED <<next, nRead, tov, twin, inLen, seen, done>>
                           ELSE /\ next' = k /\ nRead' = nr
                                /\ tov' = carry
                                /\ twin' = recs /\ inLen' = Len(tov) /\ seen' = seen \cup Range(recs)
                                /\ done' = last
                                /\ UNCHANGED pc
    /\ UNCHANGED <<T, R, rnext, rov, rwin>>

Impute == Mode = "impute"

(* advanceRefItToChrom: ref_has_next reads ahead without counting a read. *)
RECURSIVE Advance(_, _, _)
Advance(rk, nr, chrom) == IF R[rk].c # chrom /\ rk < NR THEN Advance(rk + 1, nr + 1, chrom) ELSE <<rk, nr>>

(* The inner loop of next_ref_window over the reference records before
   target record k. *)
RECURSIVE Skip(_, _, _, _, _)
Skip(ref, rk, nr, k, chrom) ==
    IF rk <= NR /\ R[rk].c = chrom /\ (R[rk].p < T[k].p \/ (R[rk].p = T[k].p /\ R[rk] # T[k]))
    THEN Skip(IF Impute THEN Append(ref, [r |-> rk, in |-> FALSE]) ELSE ref, rk + 1, nr + 1, k, chrom)
    ELSE [ref |-> ref, rk |-> rk, nr |-> nr]

(* The target loop of next_ref_window; s is the inner loop's result for
   target record k. *)
RECURSIVE RLoop(_, _, _, _, _, _, _), RMatch(_, _, _, _, _)
RLoop(targ, ref, k, rk, nr, chrom, endPos) ==
    IF k <= NT /\ T[k].c = chrom /\ T[k].p < endPos /\ Len(ref) < WM
    THEN RMatch(Skip(ref, rk, nr, k, chrom), targ, k, chrom, endPos)
    ELSE [targ |-> targ, ref |-> ref, k |-> k, rk |-> rk, nr |-> nr]

RMatch(s, targ, k, chrom, endPos) ==
    IF s.rk <= NR /\ R[s.rk] = T[k]
    THEN RLoop(Append(targ, k), Append(s.ref, [r |-> s.rk, in |-> TRUE]), k + 1, s.rk + 1, s.nr + 2, chrom, endPos)
    ELSE RLoop(targ, s.ref, k + 1, s.rk, s.nr + 1, chrom, endPos)

(* The impute=true loop over the reference records after the last target. *)
RECURSIVE Trail(_, _, _, _, _)
Trail(ref, rk, nr, chrom, endPos) ==
    IF Impute /\ rk <= NR /\ R[rk].c = chrom /\ R[rk].p < endPos /\ Len(ref) < WM
    THEN Trail(Append(ref, [r |-> rk, in |-> FALSE]), rk + 1, nr + 1, chrom, endPos)
    ELSE [ref |-> ref, rk |-> rk, nr |-> nr]

(* next_ref_window, then check_advances and the carried overlaps.
   MarkerIndices.targOverlapStart is the first target marker at or after the
   reference overlap start. *)
RefWindow ==
    /\ pc = "start" /\ ~done /\ Mode # "targ"
    /\ LET chrom == T[next].c
       IN  \E adv \in {Advance(rnext, nRead, chrom)} :
           LET rk0 == adv[1]
               nr0 == adv[2]
               endPos == R[rk0].p + (IF rov = <<>> THEN W ELSE W - O)
           IN  \E lp \in {RLoop(tov, rov, next, rk0, nr0, chrom, endPos)} :
               \E tr \in {Trail(lp.ref, lp.rk, lp.nr, chrom, endPos)} :
               LET targ == lp.targ
                   ref == tr.ref
                   k == lp.k
                   rk == tr.rk
                   nr == tr.nr
                   last == k > NT \/ rk > NR
                   chromEnd == last \/ R[ref[1].r].c # R[rk].c
               IN  IF targ = <<>> \/ ref = <<>>
                   THEN pc' = "empty" /\ UNCHANGED <<next, rnext, nRead, tov, rov, twin, rwin, inLen, done, fixed>>
                   ELSE \E ovStart \in {IF chromEnd THEN Len(ref)
                                         ELSE OverlapStart([j \in 1..Len(ref) |-> R[ref[j].r].p], endPos - 1)} :
                        LET tStart == Cardinality({j \in 1..ovStart : ref[j].in})
                            tcarry == SubSeq(targ, tStart + 1, Len(targ))
                            rcarry == SubSeq(ref, ovStart + 1, Len(ref))
                        IN  /\ fixed' = (k = next /\ rk = rk0 /\ tcarry = tov /\ rcarry = rov)
                            /\ IF ~last /\ ovStart = 0 /\ nr = nr0
                               THEN pc' = "stalled" /\ UNCHANGED <<next, rnext, nRead, tov, rov, twin, rwin, inLen, done>>
                               ELSE /\ next' = k /\ rnext' = rk /\ nRead' = nr
                                    /\ tov' = tcarry /\ rov' = rcarry
                                    /\ twin' = targ /\ rwin' = ref /\ inLen' = Len(rov)
                                    /\ done' = last
                                    /\ UNCHANGED pc
    /\ UNCHANGED <<T, R, seen>>

(* read_windows hands over NULL after the last window. *)
Finish == pc = "start" /\ done /\ pc' = "stop" /\ UNCHANGED <<T, R, next, rnext, nRead, tov, rov, twin, rwin, inLen, seen, done, fixed>>

Stop == pc # "start" /\ UNCHANGED vars

Next == TargWindow \/ RefWindow \/ Finish \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(TargWindow \/ RefWindow \/ Finish)

-----------------------------------------------------------------------------
Made == twin # <<>>

IsSuffix(s, t) == Len(s) <= Len(t) /\ s = SubSeq(t, Len(t) - Len(s) + 1, Len(t))

(* A window is a run of consecutive target records. *)
WinInterval == \A j \in 1..Len(twin) - 1 : twin[j + 1] = twin[j] + 1

(* The windows made so far cover a prefix of the target records, and all of
   them once the chain ends. *)
Covered == /\ \E m \in 0..NT : seen = 1..m
           /\ pc = "stop" => seen = 1..NT

(* A window lies on one chromosome. *)
OneChrom == \A a, b \in Range(twin) : T[a].c = T[b].c

(* The overlap carried into the next window is a suffix of the window made. *)
CarrySuffix == pc = "start" /\ Made => IsSuffix(tov, twin) /\ IsSuffix(rov, rwin)

(* Nothing is carried across a chromosome boundary or past the end. *)
BoundaryEmpty == pc = "start" /\ Made /\ (next > NT \/ T[next].c # T[twin[1]].c) => tov = <<>>

(* done is set by the window holding the last target record and no other. *)
DoneAtLast == Made => (done <=> NT \in Range(twin))

IsIn(e) == e.in

(* The in-target reference records pair, in order, with target records of
   the same marker: in the window made and in what it carries. *)
Pairs(rs, ts) == LET m == SelectSeq(rs, IsIn)
                 IN  Len(m) = Len(ts) /\ \A j \in 1..Len(ts) : R[m[j].r] = T[ts[j]]

Paired == Made => Pairs(rwin, twin) /\ Pairs(rov, tov)

(* The reference records a window read lie on one chromosome, and so do the
   ones carried into it: a window spans two chromosomes only when it starts
   with an overlap from the chromosome before (Window.chromIndex). *)
ChromsBesideCarry ==
    Mode # "targ" =>
    /\ \A a, b \in inLen + 1..Len(rwin) : R[rwin[a].r].c = R[rwin[b].r].c
    /\ \A a, b \in 1..inLen : R[rwin[a].r].c = R[rwin[b].r].c

(* check_advances rejects only a window that would repeat forever. *)
StallIsReal == pc = "stalled" => fixed

(* No window that would repeat forever is made, unless it is the last. *)
NoSilentRepeat == pc = "start" /\ Made /\ ~done => ~fixed

(* The chain ends: the last window, an error, or check_advances. *)
Terminates == <>(pc # "start")
=============================================================================
