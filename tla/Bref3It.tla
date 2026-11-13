------------------------------- MODULE Bref3It -------------------------------
(***************************************************************************)
(* The record buffer of src/bref/bref3_it.c: fill, the fill in             *)
(* bref3_it_open, and bref3_it_next. A bref3 file holds NB data blocks and *)
(* then an end-of-data block whose record count is 0. fill resets the      *)
(* buffer and reads blocks until one has a record that passes the marker   *)
(* filter or the end-of-data block sets at_end. bref3_it_next returns the  *)
(* next buffered record and calls fill when that empties the buffer, so it *)
(* reads one block ahead. Init chooses each block's records and their      *)
(* filter verdicts (no records stands for a negative count, which          *)
(* read_block also reads as an empty block), and where the file is cut     *)
(* short: a read at or past the cut fails with util_exit, as fread does at *)
(* the end of the file. The caller calls bref3_it_next twice after the     *)
(* first NULL. Each step is one call or one block read.                    *)
(***************************************************************************)
EXTENDS Naturals, Sequences

CONSTANTS NB, R

VARIABLES
    Blocks,    \* per data block, the filter verdict of each record
    CutAt,     \* the first block the file lacks or cuts short; NB + 2 for an intact file
    blk,       \* the block whose header read_int reads next; NB + 1 is the end-of-data block
    buf,       \* it->buf[0..buf_n), as records
    head,      \* it->buf_head
    atEnd,     \* it->at_end
    pending,   \* the record bref3_it_next returns once fill returns
    pc,        \* "fill", "idle" or "exit" (util_exit)
    out        \* what bref3_it_next returned, in order

vars == <<Blocks, CutAt, blk, buf, head, atEnd, pending, pc, out>>

Rec(b, j) == [k |-> "rec", b |-> b, j |-> j]
Null == [k |-> "null", b |-> 0, j |-> 0]
Exit == [k |-> "exit", b |-> 0, j |-> 0]
None == [k |-> "none", b |-> 0, j |-> 0]

Passes(r) == Blocks[r.b][r.j]

(* The records of block b that pass the filter, in file order. *)
Accepted(b) == SelectSeq([j \in 1..Len(Blocks[b]) |-> Rec(b, j)], Passes)

RECURSIVE AcceptedBefore(_)
AcceptedBefore(b) == IF b <= 1 THEN <<>> ELSE AcceptedBefore(b - 1) \o Accepted(b - 1)

DropLast(s) == SubSeq(s, 1, Len(s) - 1)

(* An intact file gives every accepted record in file order and then NULL
   on every call. A cut file exits when the reader reaches the cut, one
   block after the last record it could return, so that record is lost. *)
Expected ==
    IF CutAt <= NB + 1
    THEN DropLast(AcceptedBefore(CutAt)) \o <<Exit>>
    ELSE AcceptedBefore(NB + 1) \o <<Null, Null, Null>>

TypeOK ==
    /\ Blocks \in [1..NB -> Seq(BOOLEAN)]
    /\ CutAt \in 1..(NB + 2)
    /\ blk \in 1..(NB + 3)
    /\ buf \in Seq([k : {"rec"}, b : 1..NB, j : 1..R])
    /\ head \in 0..Len(buf)
    /\ atEnd \in BOOLEAN
    /\ pc \in {"fill", "idle", "exit"}

Init ==
    /\ Blocks \in [1..NB -> UNION {[1..n -> BOOLEAN] : n \in 0..R}]
    /\ CutAt \in 1..(NB + 2)
    /\ blk = 1
    /\ buf = <<>>
    /\ head = 0
    /\ atEnd = FALSE
    /\ pending = None
    /\ pc = "fill"
    /\ out = <<>>

Done == pc = "exit" \/ Len(out) = Len(Expected)

(* One pass of fill's loop: read_int reads a block's record count. The
   end-of-data block sets at_end; a data block's accepted records are
   added to the buffer; a read at the cut exits. *)
ReadBlock ==
    /\ pc = "fill" /\ buf = <<>> /\ ~atEnd
    /\ blk' = blk + 1
    /\ IF blk >= CutAt
         THEN /\ pc' = "exit" /\ out' = Append(out, Exit)
              /\ UNCHANGED <<buf, atEnd>>
       ELSE IF blk = NB + 1
         THEN /\ atEnd' = TRUE
              /\ UNCHANGED <<buf, pc, out>>
       ELSE /\ buf' = buf \o Accepted(blk)
            /\ UNCHANGED <<atEnd, pc, out>>
    /\ UNCHANGED <<Blocks, CutAt, head, pending>>

(* fill's loop ends, and bref3_it_next returns the record it took. *)
FillDone ==
    /\ pc = "fill" /\ (buf # <<>> \/ atEnd)
    /\ pc' = "idle"
    /\ out' = IF pending = None THEN out ELSE Append(out, pending)
    /\ pending' = None
    /\ UNCHANGED <<Blocks, CutAt, blk, buf, head, atEnd>>

(* bref3_it_next: NULL when the buffer is empty, else the next record,
   calling fill, which first resets the buffer, when that was the last. *)
Take ==
    /\ pc = "idle" /\ ~Done
    /\ IF head = Len(buf)
         THEN /\ out' = Append(out, Null)
              /\ UNCHANGED <<buf, head, pending, pc>>
       ELSE IF head + 1 = Len(buf)
         THEN /\ pending' = buf[head + 1]
              /\ buf' = <<>> /\ head' = 0 /\ pc' = "fill"
              /\ UNCHANGED out
       ELSE /\ out' = Append(out, buf[head + 1])
            /\ head' = head + 1
            /\ UNCHANGED <<buf, pending, pc>>
    /\ UNCHANGED <<Blocks, CutAt, blk, atEnd>>

Stop == Done /\ pc \in {"idle", "exit"} /\ UNCHANGED vars

Next == ReadBlock \/ FillDone \/ Take \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(ReadBlock \/ FillDone \/ Take)

-----------------------------------------------------------------------------
(* Records come out in file order, a NULL is followed only by NULLs, and a
   cut file exits where the reference says. *)
OutPrefix ==
    /\ Len(out) <= Len(Expected)
    /\ \A j \in 1..Len(out) : out[j] = Expected[j]

(* Nothing is read after the end-of-data block. *)
NoReadPastEnd == blk <= NB + 2

(* The caller sees every record and then NULL, or the exit. *)
Completes == <>(out = Expected)
=============================================================================
