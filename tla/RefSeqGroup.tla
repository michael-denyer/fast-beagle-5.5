----------------------------- MODULE RefSeqGroup -----------------------------
(***************************************************************************)
(* The sequence-coding buffer of src/vcf/ref_it.c (fill, apply_seq_coding, *)
(* flush_compressed, ref_it_next and the fill in ref_it_open) with         *)
(* seq_coder3_add and seq_coder3_flush of src/bref/seq_coder3.c. fill      *)
(* takes records from block_reader_next, which returns NULL after the last *)
(* and on every later call. A record that passes the marker filter goes to *)
(* low_freq, or to the coder with a NULL placeholder in low_freq. A new    *)
(* chromosome, a full coder (add returns false, then the add is retried    *)
(* once), and the end of the input flush the coder: its records fill the   *)
(* placeholders in order and low_freq moves to ready.                      *)
(*                                                                         *)
(* A record is a function from the H haplotypes to alleles, written with   *)
(* the alleles numbered in order of first appearance; the coder depends    *)
(* only on which haplotypes share an allele. Init chooses each record's    *)
(* chromosome, whether the filter drops it, and whether it is sequence     *)
(* coded. apply_seq_coding's bound on n_alleles is kept, so a coded record *)
(* has at most MaxCoded alleles; its bound on the major allele's count is  *)
(* a free choice, which gives the coder more inputs than the code does.    *)
(* The major allele is any allele of greatest count. The coder's sequence  *)
(* indices are modeled up to renaming: hapSeq names each sequence by its   *)
(* smallest haplotype. Each step is one record read, one flush decision,   *)
(* one add, or one call.                                                   *)
(***************************************************************************)
EXTENDS Naturals, Sequences, FiniteSets

CONSTANTS N, H, MaxNSeq, C

Haps == 1..H

(* it->max_seq_coded_alleles in ref_it_open. *)
MaxCoded == IF MaxNSeq < 255 THEN MaxNSeq ELSE 255

Range(f) == {f[x] : x \in DOMAIN f}
Count(v, a) == Cardinality({h \in Haps : v[h] = a})

(* Allele vectors with alleles numbered in order of first appearance. *)
Vectors == {v \in [Haps -> 0..(MaxCoded - 1)] :
                \A h \in Haps : v[h] > 0 => \E f \in 1..(h - 1) : v[f] = v[h] - 1}
MostFrequent(v) == {a \in Range(v) : \A b \in Range(v) : Count(v, b) <= Count(v, a)}

Zero == [h \in Haps |-> 0]
Records ==
    {[kind |-> "excl", chrom |-> 1, v |-> Zero, maj |-> 0]}
    \cup {[kind |-> "low", chrom |-> c, v |-> Zero, maj |-> 0] : c \in 1..C}
    \cup UNION {{[kind |-> "coded", chrom |-> c, v |-> v, maj |-> m] : m \in MostFrequent(v)} :
                    c \in 1..C, v \in Vectors}

Fresh == [h \in Haps |-> 1]

Item(id, g) == [k |-> "rec", id |-> id, grp |-> g]
End == [k |-> "end", id |-> 0, grp |-> {}]
None == [k |-> "none", id |-> 0, grp |-> {}]
Garbage == [k |-> "garbage", id |-> 0, grp |-> {}]

VARIABLES
    Recs,      \* the records block_reader_next returns, in file order
    i,         \* the record block_reader_next returns next; N + 1 for NULL
    cur,       \* the record fill is placing
    lastChrom, \* it->last_chrom; 0 for -1, none yet
    lowFreq,   \* it->low_freq: record numbers, 0 for a NULL placeholder
    group,     \* the coder's records, c->recs
    hapSeq,    \* c->hap_to_seq, each sequence named by its smallest haplotype
    nseq,      \* c->n_seq_map
    ready,     \* it->ready
    pending,   \* the record ref_it_next returns once fill returns
    pc,        \* "fill", "chrom", "code", "retry", "ret", "idle" or "assert"
    out,       \* what ref_it_next returned, in order
    expected   \* Expected for these records, fixed by Init

vars == <<Recs, i, cur, lastChrom, lowFreq, group, hapSeq, nseq, ready, pending, pc, out, expected>>

(* The number of distinct allele strings the haplotypes carry over the
   records of g. *)
NSeq(g) == Cardinality({[k \in 1..Len(g) |-> Recs[g[k]].v[h]] : h \in Haps})

Close(g) == IF g = <<>> THEN <<>> ELSE <<g>>

(* The groups, in order, when a coded record joins the open group unless
   the chromosome changed or the group would then carry more than MaxNSeq
   distinct allele strings. *)
RECURSIVE Greedy(_, _, _)
Greedy(j, g, c) ==
    IF j > N THEN Close(g)
    ELSE LET r == Recs[j] IN
         IF r.kind = "excl" THEN Greedy(j + 1, g, c)
         ELSE IF c # 0 /\ r.chrom # c THEN Close(g) \o Greedy(j, <<>>, r.chrom)
         ELSE IF r.kind = "low" THEN Greedy(j + 1, g, c)
         ELSE IF NSeq(Append(g, j)) <= MaxNSeq THEN Greedy(j + 1, Append(g, j), r.chrom)
         ELSE Close(g) \o Greedy(j + 1, <<j>>, r.chrom)

IsKept(j) == Recs[j].kind # "excl"

(* Every kept record in file order, each coded record with its group, and
   then NULL. *)
Expected ==
    LET kept == SelectSeq([j \in 1..N |-> j], IsKept)
        groups == Greedy(1, <<>>, 0)
        GroupOf(j) == Range(groups[CHOOSE k \in 1..Len(groups) : j \in Range(groups[k])])
    IN  [m \in 1..Len(kept) |->
            IF Recs[kept[m]].kind = "coded" THEN Item(kept[m], GroupOf(kept[m]))
            ELSE Item(kept[m], {})] \o <<End>>

Seqs == Range(hapSeq)

(* set_allele_map for record r: the new sequences its non-major alleles
   push, after the first per sequence, and add_major_allele's push for a
   sequence that holds both the major allele and another. *)
NewNSeq(r) ==
    LET NonMajor(s) == {r.v[h] : h \in {g \in Haps : hapSeq[g] = s /\ r.v[g] # r.maj}}
        HasMajor(s) == \E h \in Haps : hapSeq[h] = s /\ r.v[h] = r.maj
    IN  nseq
        + Cardinality(UNION {{<<s, a>> : a \in NonMajor(s)} : s \in Seqs})
        - Cardinality({s \in Seqs : NonMajor(s) # {}})
        + Cardinality({s \in Seqs : HasMajor(s) /\ NonMajor(s) # {}})

(* seq_coder3_add's new hap_to_seq: haplotypes stay together while they
   share the old sequence and the record's allele. *)
Refine(v) ==
    [h \in Haps |-> CHOOSE g \in Haps :
        /\ hapSeq[g] = hapSeq[h] /\ v[g] = v[h]
        /\ \A f \in Haps : (hapSeq[f] = hapSeq[h] /\ v[f] = v[h]) => g <= f]

(* seq_coder3_add succeeds for record cur and appends its placeholder. A
   false return restores n_seq_map to seq_to_cnt.n, one per sequence. *)
AddOk == NewNSeq(Recs[cur]) <= MaxNSeq
Add ==
    /\ group' = Append(group, cur)
    /\ hapSeq' = Refine(Recs[cur].v)
    /\ nseq' = NewNSeq(Recs[cur])
    /\ lowFreq' = Append(lowFreq, 0)

(* flush_compressed: seq_coder3_flush hands out the coder's records, and
   initializes the coder unless it held none; low_freq moves to ready with
   the k-th placeholder replaced by the k-th coded record. *)
Flushed(lf, g) ==
    [p \in 1..Len(lf) |->
        IF lf[p] # 0 THEN Item(lf[p], {})
        ELSE LET k == Cardinality({q \in 1..p : lf[q] = 0})
             IN  IF k <= Len(g) THEN Item(g[k], Range(g)) ELSE Garbage]
Flush ==
    /\ ready' = ready \o Flushed(lowFreq, group)
    /\ lowFreq' = <<>>
    /\ group' = <<>>
    /\ hapSeq' = IF group = <<>> THEN hapSeq ELSE Fresh
    /\ nseq' = IF group = <<>> THEN nseq ELSE 1

TypeOK ==
    /\ Recs \in [1..N -> Records]
    /\ i \in 1..(N + 1)
    /\ cur \in 0..N
    /\ lastChrom \in 0..C
    /\ lowFreq \in Seq(0..N)
    /\ group \in Seq(1..N)
    /\ hapSeq \in [Haps -> Haps]
    /\ nseq \in Nat
    /\ pc \in {"fill", "chrom", "code", "retry", "ret", "idle", "assert"}

Init ==
    /\ Recs \in [1..N -> Records]
    /\ i = 1
    /\ cur = 0
    /\ lastChrom = 0
    /\ lowFreq = <<>>
    /\ group = <<>>
    /\ hapSeq = Fresh
    /\ nseq = 1
    /\ ready = <<>>
    /\ pending = None
    /\ pc = "fill"
    /\ out = <<>>
    /\ expected = Expected

Done == Len(out) = Len(expected)

(* fill's loop with ready empty: block_reader_next gives a record, which
   the filter may drop, or NULL, which flushes and ends fill. *)
Read ==
    /\ pc = "fill" /\ ready = <<>>
    /\ IF i > N
         THEN /\ Flush /\ pc' = "ret"
              /\ UNCHANGED <<i, cur>>
       ELSE IF Recs[i].kind = "excl"
         THEN /\ i' = i + 1
              /\ UNCHANGED <<cur, lowFreq, group, hapSeq, nseq, ready, pc>>
       ELSE /\ i' = i + 1 /\ cur' = i /\ pc' = "chrom"
            /\ UNCHANGED <<lowFreq, group, hapSeq, nseq, ready>>
    /\ UNCHANGED <<Recs, lastChrom, pending, out, expected>>

(* A record on a new chromosome flushes the coder. *)
Chrom ==
    /\ pc = "chrom"
    /\ IF lastChrom # 0 /\ Recs[cur].chrom # lastChrom
         THEN Flush
         ELSE UNCHANGED <<lowFreq, group, hapSeq, nseq, ready>>
    /\ lastChrom' = Recs[cur].chrom
    /\ pc' = "code"
    /\ UNCHANGED <<Recs, i, cur, pending, out, expected>>

(* apply_seq_coding false: the record goes to low_freq. True: it goes to
   the coder, and a full coder is flushed before the retry. *)
Code ==
    /\ pc = "code"
    /\ IF Recs[cur].kind = "low"
         THEN /\ lowFreq' = Append(lowFreq, cur) /\ pc' = "fill"
              /\ UNCHANGED <<group, hapSeq, nseq, ready>>
       ELSE IF AddOk
         THEN /\ Add /\ pc' = "fill"
              /\ UNCHANGED ready
       ELSE /\ Flush /\ pc' = "retry"
    /\ UNCHANGED <<Recs, i, cur, lastChrom, pending, out, expected>>

(* The add after the flush; a second false is util_exit
   "java.lang.AssertionError: SeqCoder3.add". *)
Retry ==
    /\ pc = "retry"
    /\ IF AddOk
         THEN /\ Add /\ pc' = "fill"
         ELSE /\ pc' = "assert"
              /\ UNCHANGED <<lowFreq, group, hapSeq, nseq>>
    /\ UNCHANGED <<Recs, i, cur, lastChrom, ready, pending, out, expected>>

(* fill returns, and ref_it_next returns the record it took. *)
Return ==
    /\ pc = "ret" \/ (pc = "fill" /\ ready # <<>>)
    /\ pc' = "idle"
    /\ out' = IF pending = None THEN out ELSE Append(out, pending)
    /\ pending' = None
    /\ UNCHANGED <<Recs, i, cur, lastChrom, lowFreq, group, hapSeq, nseq, ready, expected>>

(* ref_it_next: NULL when ready is empty, else its first record, calling
   fill when that empties ready. *)
Take ==
    /\ pc = "idle" /\ ~Done
    /\ IF ready = <<>>
         THEN /\ out' = Append(out, End)
              /\ UNCHANGED <<ready, pending, pc>>
       ELSE IF Len(ready) = 1
         THEN /\ pending' = Head(ready) /\ ready' = <<>> /\ pc' = "fill"
              /\ UNCHANGED out
       ELSE /\ out' = Append(out, Head(ready)) /\ ready' = Tail(ready)
            /\ UNCHANGED <<pending, pc>>
    /\ UNCHANGED <<Recs, i, cur, lastChrom, lowFreq, group, hapSeq, nseq, expected>>

Stop == Done /\ pc = "idle" /\ UNCHANGED vars

Next == Read \/ Chrom \/ Code \/ Retry \/ Return \/ Take \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Read \/ Chrom \/ Code \/ Retry \/ Return \/ Take)

-----------------------------------------------------------------------------
(* The kept records come out in file order, each coded record with the
   group the reference puts it in, and then NULL. *)
OutPrefix ==
    /\ Len(out) <= Len(expected)
    /\ \A j \in 1..Len(out) : out[j] = expected[j]

(* Each placeholder in low_freq has one coded record in the coder. *)
PlaceholdersFilled == Cardinality({p \in 1..Len(lowFreq) : lowFreq[p] = 0}) = Len(group)

(* n_seq_map counts the distinct allele strings of the coder's records,
   one per sequence a haplotype carries, and stays within max_nseq. *)
CoderSeqs ==
    /\ nseq = NSeq(group)
    /\ nseq = Cardinality(Seqs)
    /\ nseq <= MaxNSeq

(* An empty coder is in its initial state, which seq_coder3_flush relies
   on when it returns early. *)
EmptyCoderFresh == group = <<>> => hapSeq = Fresh /\ nseq = 1

(* A group never spans two chromosomes. *)
GroupOneChrom == \A j, k \in Range(group) : Recs[j].chrom = Recs[k].chrom

(* The retried add never fails. *)
NoAssert == pc # "assert"

(* The caller sees every record and then NULL. *)
Completes == <>(out = expected)
=============================================================================
