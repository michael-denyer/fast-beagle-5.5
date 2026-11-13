------------------------------- MODULE HmmPass -------------------------------
(***************************************************************************)
(* One phase_baum2 (pb) phasing S samples in turn, from                    *)
(* src/phase/phase_baum2.c: begin_sample, bwd_begin, bwd_step with         *)
(* bwd_end_step, fwd_begin, fwd_step with fwd_end_step, phase_het,         *)
(* swap_haps, impute_alleles and end_sample. Each sample has clusters      *)
(* 0 .. N-1. Init chooses the run's threshold regime ("burnin" with        *)
(* pb->burnin set and lr infinite, "finite", or "one" on the last          *)
(* iteration) and each sample's cluster types when begin_sample has        *)
(* masked them. A PHASED_HET precedes every UNPHASED_HET and MASKED_HET,   *)
(* since PbwtPhaser.indices never lists a sample's first heterozygote, and *)
(* no MASKED_HET exists in burn-in, since masking needs lr < 50            *)
(* (PhasingSchedule checks that burn-in runs use the infinite threshold).  *)
(* At each unphased het the environment chooses num < den and whether the  *)
(* likelihood-ratio test passes; at each masked het, whether               *)
(* impute_masked_het's test passes. At lr = 1.0 both tests always pass     *)
(* (one side of each comparison holds), NaN aside. In burn-in and between  *)
(* the tests are free choices, so phase_het's !pb->burnin guard is checked *)
(* on its own, although with finite num and den the comparison is false at *)
(* infinity too. The backward values are not modeled, only the indices the *)
(* backward pass stores them under and the forward pass reads them from.   *)
(* Each step is one bwd_step with bwd_end_step, one fwd_step with          *)
(* fwd_end_step, or a sample's begin or end.                               *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets

CONSTANTS N, S

ASSUME N \in Nat \ {0} /\ S \in Nat \ {0}

VARIABLES
    regime,    \* "burnin", "finite" or "one"
    types0,    \* each sample's clust_type when the run starts
    types,     \* each sample's clust_type
    pc,        \* "begin", "bwd", "fwd" or "done"
    s,         \* the sample being phased
    snap,      \* clust_type of sample s at begin_sample
    ran,       \* samples whose begin_sample returned true
    swapHaps,  \* pb->swap_haps
    nSwaps,    \* pb->n_swaps
    rateSwaps, \* rate->n_swaps
    rateUnph,  \* rate->n_unph_hets
    outs,      \* per sample, the num < den outcome at each unphased het, in order
    bc, bMiss, bUnph,  \* bwd_pass c, miss_index, unph_het_index
    fc, fMiss, fUnph,  \* fwd_pass c, miss_index, unph_het_index
    hetW, missW,       \* per cluster, the index bwd stored its values under; -1 for none
    hetR, missR,       \* per cluster, the index fwd read them from; -1 for none
    swaps      \* per cluster, how many swap_haps calls covered it in sample s

vars == <<regime, types0, types, pc, s, snap, ran, swapHaps, nSwaps, rateSwaps, rateUnph,
          outs, bc, bMiss, bUnph, fc, fMiss, fUnph, hetW, missW, hetR, missR, swaps>>

Cl == 0 .. N - 1
Samples == 1 .. S
Types == {"MISSING", "MASKED", "HOM", "PHASED", "UNPH"}
None == -1
NoIdx == [k \in Cl |-> None]

MissOrMasked(t) == t \in {"MISSING", "MASKED"}
CountUnph(t, K) == Cardinality({k \in K : t[k] = "UNPH"})
CountMM(t, K) == Cardinality({k \in K : MissOrMasked(t[k])})

(* mc->unph_het_clusters: the unphased hets at begin_sample, in order. *)
RECURSIVE UnphFrom(_, _)
UnphFrom(t, k) == IF k >= N THEN <<>>
                  ELSE IF t[k] = "UNPH" THEN <<k>> \o UnphFrom(t, k + 1)
                  ELSE UnphFrom(t, k + 1)
UnphClusters(t) == UnphFrom(t, 0)

(* The environment's verdict on a likelihood-ratio test, by regime. *)
Ratio(e) == IF regime = "one" THEN TRUE ELSE e

(* Reference definitions over the environment. *)
Rank(t, k) == CountUnph(t, 0 .. k - 1)
MissRank(t, k) == CountMM(t, 0 .. k - 1)
Flips(o) == Cardinality({i \in 1 .. Len(o) : o[i] # (IF i = 1 THEN FALSE ELSE o[i - 1])})
LastOr(o) == IF o = <<>> THEN FALSE ELSE o[Len(o)]
Max(X) == CHOOSE x \in X : \A y \in X : y <= x
(* How many times cluster k is swapped once the hets before fc are phased: the
   nearest unphased het at or before k decides, if it has been phased. *)
SwapRef(t, o, k) ==
    LET H == {j \in 0 .. k : t[j] = "UNPH"}
    IN  IF H = {} THEN 0
        ELSE LET j == Max(H)
             IN  IF j < fc /\ o[Rank(t, j) + 1] THEN 1 ELSE 0
SumRan(f(_)) ==
    LET sum[i \in 0 .. S] == IF i = 0 THEN 0 ELSE sum[i - 1] + (IF i \in ran THEN f(i) ELSE 0)
    IN  sum[S]

InitTypes(r) ==
    {t \in [Cl -> Types] :
        \A c \in Cl :
            /\ t[c] \in {"UNPH", "MASKED"} => \E d \in Cl : d < c /\ t[d] = "PHASED"
            /\ r = "burnin" => t[c] # "MASKED"}

TypeOK ==
    /\ regime \in {"burnin", "finite", "one"}
    /\ types0 \in [Samples -> [Cl -> Types]] /\ types \in [Samples -> [Cl -> Types]]
    /\ snap \in [Cl -> Types]
    /\ pc \in {"begin", "bwd", "fwd", "done"}
    /\ s \in 1 .. S + 1
    /\ ran \subseteq Samples
    /\ swapHaps \in BOOLEAN /\ nSwaps \in Nat
    /\ rateSwaps \in Nat /\ rateUnph \in Nat
    /\ outs \in [Samples -> Seq(BOOLEAN)]
    /\ bc \in -1 .. N /\ bMiss \in -1 .. N /\ bUnph \in -1 .. N
    /\ fc \in 0 .. N /\ fMiss \in 0 .. N /\ fUnph \in 0 .. N
    /\ hetW \in [Cl -> -1 .. N] /\ missW \in [Cl -> -1 .. N]
    /\ hetR \in [Cl -> -1 .. N] /\ missR \in [Cl -> -1 .. N]
    /\ swaps \in [Cl -> Nat]

Init ==
    /\ regime \in {"burnin", "finite", "one"}
    /\ types0 \in [Samples -> InitTypes(regime)]
    /\ types = types0
    /\ pc = "begin" /\ s = 1
    /\ snap = types0[1]
    /\ ran = {}
    /\ swapHaps = FALSE /\ nSwaps = 0
    /\ rateSwaps = 0 /\ rateUnph = 0
    /\ outs = [x \in Samples |-> <<>>]
    /\ bc = None /\ bMiss = None /\ bUnph = None
    /\ fc = 0 /\ fMiss = 0 /\ fUnph = 0
    /\ hetW = NoIdx /\ missW = NoIdx /\ hetR = NoIdx /\ missR = NoIdx
    /\ swaps = [k \in Cl |-> 0]

(* begin_sample, and bwd_begin when it returns true. *)
Begin ==
    /\ pc = "begin" /\ s <= S
    /\ LET t == types[s]
           nUnph == CountUnph(t, Cl)
           nMM == CountMM(t, Cl)
       IN  IF nMM = 0 /\ nUnph = 0
           THEN /\ s' = s + 1
                /\ UNCHANGED <<pc, snap, ran, swapHaps, nSwaps, bc, bMiss, bUnph,
                               hetW, missW, hetR, missR, swaps>>
           ELSE /\ pc' = "bwd" /\ snap' = t /\ ran' = ran \cup {s}
                /\ nSwaps' = 0 /\ swapHaps' = FALSE
                /\ hetR' = NoIdx /\ missR' = NoIdx /\ swaps' = [k \in Cl |-> 0]
                /\ bUnph' = nUnph - 1
                /\ IF MissOrMasked(t[N - 1])
                     THEN missW' = [NoIdx EXCEPT ![N - 1] = nMM - 1] /\ bMiss' = nMM - 2
                     ELSE missW' = NoIdx /\ bMiss' = nMM - 1
                /\ hetW' = NoIdx
                /\ bc' = N - 2
                /\ UNCHANGED s
    /\ UNCHANGED <<regime, types0, types, rateSwaps, rateUnph, outs, fc, fMiss, fUnph>>

(* bwd_step and bwd_end_step: the values reach cluster bc, which stores them
   if it is missing or masked, and cluster bc + 1 if it is an unphased het. *)
BwdStep ==
    /\ pc = "bwd" /\ bc >= 0
    /\ LET t == types[s]
           mm == MissOrMasked(t[bc])
           het == t[bc + 1] = "UNPH"
       IN  /\ missW' = IF mm THEN [missW EXCEPT ![bc] = bMiss] ELSE missW
           /\ bMiss' = IF mm THEN bMiss - 1 ELSE bMiss
           /\ hetW' = IF het THEN [hetW EXCEPT ![bc + 1] = bUnph] ELSE hetW
           /\ bUnph' = IF het THEN bUnph - 1 ELSE bUnph
    /\ bc' = bc - 1
    /\ UNCHANGED <<regime, types0, types, pc, s, snap, ran, swapHaps, nSwaps, rateSwaps, rateUnph,
                   outs, fc, fMiss, fUnph, hetR, missR, swaps>>

(* bwd_finish returns; fwd_begin. *)
FwdBegin ==
    /\ pc = "bwd" /\ bc < 0
    /\ pc' = "fwd" /\ fc' = 0 /\ fMiss' = 0 /\ fUnph' = 0
    /\ UNCHANGED <<regime, types0, types, s, snap, ran, swapHaps, nSwaps, rateSwaps, rateUnph,
                   outs, bc, bMiss, bUnph, hetW, missW, hetR, missR, swaps>>

(* fwd_step at cluster fc, with phase_het and swap_haps at an unphased het,
   then fwd_end_step with impute_alleles at a missing or masked cluster. *)
FwdStep ==
    /\ pc = "fwd" /\ fc < N
    /\ \E swap, pass, passMasked \in BOOLEAN :
         LET c == fc
             t == types[s]
             hets == UnphClusters(snap)
             swapEnd == IF fUnph + 1 < Len(hets) THEN hets[fUnph + 2] ELSE N
             isHet == t[c] = "UNPH"
             markHet == isHet /\ regime # "burnin" /\ Ratio(pass)
             t1 == IF markHet THEN [t EXCEPT ![c] = "PHASED"] ELSE t
             isMM == MissOrMasked(t1[c])
             markMasked == t1[c] = "MASKED" /\ Ratio(passMasked)
             t2 == IF markMasked THEN [t1 EXCEPT ![c] = "PHASED"] ELSE t1
         IN  /\ isHet \/ (swap = FALSE /\ pass = FALSE)
             /\ isMM \/ passMasked = FALSE
             /\ types' = [types EXCEPT ![s] = t2]
             /\ IF isHet
                  THEN /\ hetR' = [hetR EXCEPT ![c] = fUnph]
                       /\ swapHaps' = swap
                       /\ nSwaps' = IF swap # swapHaps THEN nSwaps + 1 ELSE nSwaps
                       /\ outs' = [outs EXCEPT ![s] = Append(@, swap)]
                       /\ fUnph' = fUnph + 1
                       /\ swaps' = [k \in Cl |-> IF swap /\ c <= k /\ k < swapEnd
                                                 THEN swaps[k] + 1 ELSE swaps[k]]
                  ELSE UNCHANGED <<hetR, swapHaps, nSwaps, outs, fUnph, swaps>>
             /\ IF isMM
                  THEN missR' = [missR EXCEPT ![c] = fMiss] /\ fMiss' = fMiss + 1
                  ELSE UNCHANGED <<missR, fMiss>>
    /\ fc' = fc + 1
    /\ UNCHANGED <<regime, types0, pc, s, snap, ran, rateSwaps, rateUnph, bc, bMiss, bUnph, hetW, missW>>

(* fwd_finish returns; end_sample adds the sample's counts to the rate. *)
EndSample ==
    /\ pc = "fwd" /\ fc = N
    /\ rateSwaps' = rateSwaps + nSwaps
    /\ rateUnph' = rateUnph + Len(UnphClusters(snap))
    /\ pc' = "begin" /\ s' = s + 1
    /\ UNCHANGED <<regime, types0, types, snap, ran, swapHaps, nSwaps, outs, bc, bMiss, bUnph,
                   fc, fMiss, fUnph, hetW, missW, hetR, missR, swaps>>

Finish ==
    /\ pc = "begin" /\ s > S
    /\ pc' = "done"
    /\ UNCHANGED <<regime, types0, types, s, snap, ran, swapHaps, nSwaps, rateSwaps, rateUnph, outs,
                   bc, bMiss, bUnph, fc, fMiss, fUnph, hetW, missW, hetR, missR, swaps>>

Stop == pc = "done" /\ UNCHANGED vars

Next == Begin \/ BwdStep \/ FwdBegin \/ FwdStep \/ EndSample \/ Finish \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Begin \/ BwdStep \/ FwdBegin \/ FwdStep \/ EndSample \/ Finish)

-----------------------------------------------------------------------------
(* The forward cursors count the unphased hets and the missing or masked
   clusters before c, by the types at begin_sample. *)
FwdLockstep ==
    pc = "fwd" => fUnph = CountUnph(snap, 0 .. fc - 1) /\ fMiss = CountMM(snap, 0 .. fc - 1)

(* The backward cursors count down to the ones before the clusters still to
   store. *)
BwdLockstep ==
    pc = "bwd" => /\ bUnph = CountUnph(snap, 0 .. bc + 1) - 1
                  /\ bMiss = CountMM(snap, 0 .. bc) - 1

(* The forward pass reads each cluster's backward values from the index the
   backward pass stored them under, its rank among clusters of its kind. *)
IndicesMatch ==
    \A k \in Cl :
        /\ hetW[k] # None => hetW[k] = Rank(snap, k)
        /\ missW[k] # None => missW[k] = MissRank(snap, k)
        /\ hetR[k] # None => hetR[k] = hetW[k] /\ hetR[k] = Rank(snap, k)
        /\ missR[k] # None => missR[k] = missW[k] /\ missR[k] = MissRank(snap, k)

(* The pass changes no cluster it has not reached. *)
AheadUnchanged == pc = "fwd" => \A k \in fc .. N - 1 : types[s][k] = snap[k]

(* n_swaps counts the changes of swap_haps, which starts each sample false. *)
SwapsCountFlips ==
    pc = "fwd" => nSwaps = Flips(outs[s]) /\ swapHaps = LastOr(outs[s])

(* The swap rate sums every phased sample's flips and unphased hets. *)
RateSums ==
    pc = "done" => /\ rateSwaps = SumRan(LAMBDA x : Flips(outs[x]))
                   /\ rateUnph = SumRan(LAMBDA x : CountUnph(types0[x], Cl))

(* Each swap covers exactly the clusters from its het to the next unphased
   het, so each cluster is swapped as the het before it decided. *)
SwapRange == pc = "fwd" => \A k \in Cl : swaps[k] = SwapRef(snap, outs[s], k)

(* Nothing is marked phased in burn-in. *)
NoMarkInBurnin == regime = "burnin" => types = types0

(* At the last iteration every unphased and masked het ends phased. *)
NoneLeft == \A x \in Samples, k \in Cl : types[x][k] \notin {"UNPH", "MASKED"}
LastRunPhasesAll == <>(pc = "done" /\ (regime = "one" => NoneLeft))
=============================================================================
