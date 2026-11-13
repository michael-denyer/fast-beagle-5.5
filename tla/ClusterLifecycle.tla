--------------------------- MODULE ClusterLifecycle ---------------------------
(***************************************************************************)
(* The cluster types of one sample in src/phase/sample_phase.c:            *)
(* set_clusters, set_clust_type, sample_phase_mark_unphased_het_as_phased, *)
(* sample_phase_mark_masked_het_as_phased, mask_run and                    *)
(* sample_phase_mask_trailing_unphased_hets. No other file writes          *)
(* clust_type or clust_type_cnt. Clusters are 0 .. N-1. Init chooses the   *)
(* types set_clusters can produce: no MASKED_HET, and a PHASED_HET before  *)
(* every UNPHASED_HET, because PbwtPhaser.indices (pbwt_phaser.c) never    *)
(* lists a sample's first heterozygote as unphased. The masking pass runs  *)
(* at the start of a sample in each iteration whose threshold is below 50  *)
(* (phase_baum2.c begin_sample), at most MaxPasses times. Between passes   *)
(* the HMM marks a cluster phased after testing its type (fwd_step,        *)
(* impute_alleles), when the environment says the likelihood ratio passed. *)
(* Whether a run of three spans at most 3000 bp is an environment choice   *)
(* made at the start of each pass, per first het of the run. Each step is  *)
(* one pass of the masking loop with its mask_run, or one mark.            *)
(***************************************************************************)
EXTENDS Naturals, Sequences, FiniteSets

CONSTANTS N, MaxPasses

ASSUME N \in Nat /\ MaxPasses \in Nat

VARIABLES
    type,      \* sp->clust_type
    cnt,       \* sp->clust_type_cnt
    exited,    \* set_clust_type reached util_exit
    pc,        \* "idle", "mask" (inside the masking loop) or "masked"
    passes,    \* masking passes started
    mc,        \* the loop's c
    run,       \* the loop's clusters list
    snap,      \* clust_type when the current masking pass started
    fits       \* per unphased het, whether a 3-run starting there spans <= 3000 bp

vars == <<type, cnt, exited, pc, passes, mc, run, snap, fits>>

Cl == 0 .. N - 1
Types == {"MISSING", "MASKED", "HOM", "PHASED", "UNPH"}
Allowed == {<<"UNPH", "PHASED">>, <<"MASKED", "PHASED">>, <<"UNPH", "MASKED">>}

Count(t, x) == Cardinality({c \in Cl : t[c] = x})
Counts(t) == [x \in Types |-> Count(t, x)]

(* set_clust_type: util_exit unless the cluster has type from. *)
SetType(st, c, from, to) ==
    IF st.type[c] # from THEN [st EXCEPT !.exited = TRUE]
    ELSE [type |-> [st.type EXCEPT ![c] = to],
          cnt |-> [st.cnt EXCEPT ![from] = @ - 1, ![to] = @ + 1],
          exited |-> st.exited]

RECURSIVE MaskAll(_, _)
MaskAll(st, cs) ==
    IF cs = <<>> THEN st
    ELSE MaskAll(SetType(st, Head(cs), "UNPH", "MASKED"), Tail(cs))

(* mask_run on clusters, with fit the outcome of its bp test. *)
MaskRun(st, cs, fit) ==
    LET lastMasked == Len(cs) - 2
    IN  IF lastMasked = 0 THEN MaskAll(st, <<cs[1]>>)
        ELSE IF lastMasked > 0 /\ fit THEN MaskAll(st, SubSeq(cs, 1, lastMasked + 1))
        ELSE st

Cur == [type |-> type, cnt |-> cnt, exited |-> exited]

(* The masking a pass must do, by segments: a segment is a maximal interval
   of clusters with no PHASED_HET in it. In a segment with two or three
   UNPHASED_HETs every one but the last is masked, for three only when the
   bp test passed for the segment's first one. *)
Segment(t, c) == {d \in Cl : \A e \in Cl : (e >= c /\ e <= d) \/ (e >= d /\ e <= c) => t[e] # "PHASED"}
Unph(t, c) == {d \in Segment(t, c) : t[d] = "UNPH"}
Max(S) == CHOOSE x \in S : \A y \in S : y <= x
Min(S) == CHOOSE x \in S : \A y \in S : x <= y
MaskRef(t, f) ==
    [c \in Cl |->
        IF t[c] = "UNPH"
           /\ Cardinality(Unph(t, c)) \in {2, 3}
           /\ c # Max(Unph(t, c))
           /\ (Cardinality(Unph(t, c)) = 2 \/ f[Min(Unph(t, c))])
        THEN "MASKED" ELSE t[c]]

InitTypes ==
    {t \in [Cl -> Types \ {"MASKED"}] :
        \A c \in Cl : t[c] = "UNPH" => \E d \in Cl : d < c /\ t[d] = "PHASED"}

TypeOK ==
    /\ type \in [Cl -> Types]
    /\ cnt \in [Types -> Nat]
    /\ exited \in BOOLEAN
    /\ pc \in {"idle", "mask", "masked"}
    /\ passes \in 0 .. MaxPasses
    /\ mc \in 0 .. N
    /\ run \in Seq(Cl)
    /\ fits \in [Cl -> BOOLEAN]

Init ==
    /\ type \in InitTypes
    /\ cnt = Counts(type)
    /\ exited = FALSE
    /\ pc = "idle"
    /\ passes = 0
    /\ mc = 0
    /\ run = <<>>
    /\ snap = type
    /\ fits = [c \in Cl |-> FALSE]

(* fwd_step: the cluster is an unphased het and phase_het marks it. *)
MarkUnphased(c) ==
    /\ pc = "idle" /\ ~exited /\ type[c] = "UNPH"
    /\ LET st == SetType(Cur, c, "UNPH", "PHASED")
       IN type' = st.type /\ cnt' = st.cnt /\ exited' = st.exited
    /\ UNCHANGED <<pc, passes, mc, run, snap, fits>>

(* impute_alleles: the cluster is a masked het and impute_masked_het marks it. *)
MarkMasked(c) ==
    /\ pc = "idle" /\ ~exited /\ type[c] = "MASKED"
    /\ LET st == SetType(Cur, c, "MASKED", "PHASED")
       IN type' = st.type /\ cnt' = st.cnt /\ exited' = st.exited
    /\ UNCHANGED <<pc, passes, mc, run, snap, fits>>

(* begin_sample calls sample_phase_mask_trailing_unphased_hets. *)
MaskBegin ==
    /\ pc = "idle" /\ ~exited /\ passes < MaxPasses
    /\ pc' = "mask" /\ passes' = passes + 1 /\ mc' = 0 /\ run' = <<>>
    /\ snap' = type
    /\ fits' \in {f \in [Cl -> BOOLEAN] : \A c \in Cl : f[c] => type[c] = "UNPH"}
    /\ UNCHANGED <<type, cnt, exited>>

(* One pass of the loop over clusters. *)
MaskStep ==
    /\ pc = "mask" /\ ~exited /\ mc < N
    /\ IF type[mc] = "PHASED"
         THEN /\ IF 2 <= Len(run) /\ Len(run) <= 3
                   THEN LET st == MaskRun(Cur, run, fits[run[1]])
                        IN type' = st.type /\ cnt' = st.cnt /\ exited' = st.exited
                   ELSE UNCHANGED <<type, cnt, exited>>
              /\ run' = <<>>
         ELSE /\ run' = IF type[mc] = "UNPH" THEN Append(run, mc) ELSE run
              /\ UNCHANGED <<type, cnt, exited>>
    /\ mc' = mc + 1
    /\ UNCHANGED <<pc, passes, snap, fits>>

(* The mask_run after the loop. *)
MaskEnd ==
    /\ pc = "mask" /\ ~exited /\ mc = N
    /\ IF 2 <= Len(run) /\ Len(run) <= 3
         THEN LET st == MaskRun(Cur, run, fits[run[1]])
              IN type' = st.type /\ cnt' = st.cnt /\ exited' = st.exited
         ELSE UNCHANGED <<type, cnt, exited>>
    /\ pc' = "masked"
    /\ UNCHANGED <<passes, mc, run, snap, fits>>

(* sample_phase_mask_trailing_unphased_hets returns. *)
MaskReturn ==
    /\ pc = "masked"
    /\ pc' = "idle"
    /\ UNCHANGED <<type, cnt, exited, passes, mc, run, snap, fits>>

Done == exited \/ (pc = "idle" /\ passes = MaxPasses)

Stop == Done /\ UNCHANGED vars

Next ==
    \/ \E c \in Cl : MarkUnphased(c) \/ MarkMasked(c)
    \/ MaskBegin \/ MaskStep \/ MaskEnd \/ MaskReturn
    \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(MaskBegin \/ MaskStep \/ MaskEnd \/ MaskReturn)

-----------------------------------------------------------------------------
(* clust_type_cnt holds the count of each type. *)
CountsMatch == cnt = Counts(type)

(* No set_clust_type finds a type other than the one its caller expects. *)
NoExit == ~exited

(* A masking pass masks exactly the leading hets of each short run. *)
MaskMatchesRef == pc = "masked" => type = MaskRef(snap, fits)

(* Only UNPHASED to PHASED, MASKED to PHASED and UNPHASED to MASKED; MISSING,
   HOMOZYGOUS and PHASED never change. *)
Transitions ==
    [][\A c \in Cl : type'[c] = type[c] \/ <<type[c], type'[c]>> \in Allowed]_vars

(* Every masking pass returns. *)
Terminates == <>Done
=============================================================================
