------------------------------- MODULE VcfIndex -------------------------------
(***************************************************************************)
(* The .tbi bookkeeping of src/main/vcf_index.c: vcf_index_new,            *)
(* vcf_index_add (with parse) and vcf_index_write, as window_writer.c      *)
(* calls them. window_writer_open writes the header lines and then calls   *)
(* vcf_index_new with the offset after them, so header lines never reach   *)
(* vcf_index_add. put (window_writer.c:421) calls vcf_index_add for each   *)
(* record line, with the offset after its newline, and window_writer_close *)
(* calls vcf_index_write. A record line reaches the '#' branch only when   *)
(* its CHROM starts with '#'. Offsets are counted in lines: 0 is the end   *)
(* of the header lines and line i ends at i. Init chooses each line's      *)
(* kind. A line that does not parse is not something put can write, since  *)
(* POS is an integer; Kinds can include "bad" to see what the C does with  *)
(* one. The reference is htslib's tbx_index, which ends the header at the  *)
(* last '#' line before the first other line, skips later '#' lines, and   *)
(* fails on a line it cannot parse.                                        *)
(***************************************************************************)
EXTENDS Naturals, Sequences

CONSTANTS N, Kinds

ASSUME N \in Nat /\ Kinds \subseteq {"hash", "rec", "bad"}

VARIABLES
    Kind,     \* per line: "hash" (CHROM starts with '#'), "rec" or "bad"
    next,     \* the next line put passes to vcf_index_add
    hdrEnd,   \* x->u_header_end
    recs,     \* x->recs, as line numbers
    pc        \* "adding" or "written" (vcf_index_write ran)

vars == <<Kind, next, hdrEnd, recs, pc>>

TypeOK ==
    /\ Kind \in [1..N -> Kinds]
    /\ next \in 1..N + 1
    /\ hdrEnd \in 0..N
    /\ recs \in Seq(1..N)
    /\ pc \in {"adding", "written"}

Init ==
    /\ Kind \in [1..N -> Kinds]
    /\ next = 1
    /\ hdrEnd = 0
    /\ recs = <<>>
    /\ pc = "adding"

(* vcf_index_add for line next. *)
Add ==
    /\ pc = "adding" /\ next <= N
    /\ CASE Kind[next] = "hash" ->
              /\ hdrEnd' = IF recs = <<>> THEN next ELSE hdrEnd
              /\ UNCHANGED recs
         [] Kind[next] = "bad" -> UNCHANGED <<hdrEnd, recs>>
         [] OTHER -> recs' = Append(recs, next) /\ UNCHANGED hdrEnd
    /\ next' = next + 1
    /\ UNCHANGED <<Kind, pc>>

(* window_writer_close: vcf_index_write builds the index from hdrEnd and
   recs. *)
Write ==
    /\ pc = "adding" /\ next = N + 1
    /\ pc' = "written"
    /\ UNCHANGED <<Kind, next, hdrEnd, recs>>

Stop == next = N + 1 /\ UNCHANGED vars

Next == Add \/ Write \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Add \/ Write)

-----------------------------------------------------------------------------
(* tbx_index over the lines put wrote. *)
TabixFails == \E i \in 1..N : Kind[i] = "bad"

FirstOther == IF \E i \in 1..N : Kind[i] # "hash"
              THEN CHOOSE i \in 1..N : Kind[i] # "hash" /\ \A j \in 1..i - 1 : Kind[j] = "hash"
              ELSE N + 1

TabixHeaderEnd == FirstOther - 1

RECURSIVE Others(_)
Others(i) == IF i = 0 THEN <<>>
             ELSE IF Kind[i] = "hash" THEN Others(i - 1)
             ELSE Append(Others(i - 1), i)

(* The records indexed so far are tabix's records among the lines added so
   far, in file order. *)
IndexedLikeTabix == ~TabixFails => recs = Others(next - 1)

(* The written index starts where tabix's does. *)
HeaderEndLikeTabix == pc = "written" /\ ~TabixFails => hdrEnd = TabixHeaderEnd

(* window_writer_close always writes the .tbi. *)
Completes == <>(pc = "written")
=============================================================================
