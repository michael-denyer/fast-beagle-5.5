-------------------------------- MODULE VcfIt --------------------------------
(***************************************************************************)
(* The batch iterator of src/vcf/vcf_it.c: read_line, read_lines,          *)
(* fill_batch and vcf_it_next. The iterator holds the next unread line     *)
(* (line, have_line), moves up to Batch lines into the batch, and reading  *)
(* the line after each one can fail. read_lines runs under util_try, so a  *)
(* failed read ends the batch and is kept in read_error. vcf_it_next       *)
(* returns the records of the batch in order. It raises a line's parse     *)
(* error when that record is reached, and read_error when the batch's last *)
(* record is reached, before the marker filter. Lines 1..N are the data    *)
(* lines after the header. Init chooses each line's fate (Kind) and the    *)
(* line whose following read fails (FailAt, 0 for none). Reading line 1    *)
(* happens in vcf_it_open, outside util_try, and is not modeled. Each step *)
(* is one line read or one record taken.                                   *)
(***************************************************************************)
EXTENDS Naturals, Sequences

CONSTANTS N, Batch

ASSUME N \in Nat /\ Batch \in Nat \ {0}

VARIABLES
    Kind,      \* per line: "ok", "excl" (fails the marker filter) or "bad"
    FailAt,    \* reading the line after this one fails; 0 for none
    pos,       \* the file line in it->line
    have,      \* it->have_line
    batch,     \* it->lines[0..n), as file line numbers
    next,      \* it->next
    readErr,   \* it->read_error != NULL
    pc,        \* "idle", "fill" (inside read_lines) or "take"
    out        \* what vcf_it_next returned or raised, in order

vars == <<Kind, FailAt, pos, have, batch, next, readErr, pc, out>>

Rec(i) == [k |-> "rec", line |-> i]
ParseErr(i) == [k |-> "parse", line |-> i]
ReadErr(i) == [k |-> "read", line |-> i]
End == [k |-> "end", line |-> 0]
Raised(o) == o.k \in {"parse", "read"}

(* The same records and errors when every line is read, parsed and filtered
   on demand, as the comment on vcf_it says. *)
RECURSIVE OnDemand(_)
OnDemand(i) ==
    IF i > N THEN <<End>>
    ELSE IF Kind[i] = "bad" THEN <<ParseErr(i)>>
    ELSE IF FailAt = i THEN <<ReadErr(i)>>
    ELSE IF Kind[i] = "excl" THEN OnDemand(i + 1)
    ELSE <<Rec(i)>> \o OnDemand(i + 1)

Expected == OnDemand(1)

TypeOK ==
    /\ Kind \in [1..N -> {"ok", "excl", "bad"}]
    /\ FailAt \in 0..N
    /\ pos \in 1..N + 1
    /\ have \in BOOLEAN
    /\ batch \in Seq(1..N) /\ Len(batch) <= Batch
    /\ next \in 0..Batch
    /\ readErr \in BOOLEAN
    /\ pc \in {"idle", "fill", "take"}

Init ==
    /\ Kind \in [1..N -> {"ok", "excl", "bad"}]
    /\ FailAt \in 0..N
    /\ pos = 1
    /\ have = (N >= 1)
    /\ batch = <<>>
    /\ next = 0
    /\ readErr = FALSE
    /\ pc = "idle"
    /\ out = <<>>

Done == out # <<>> /\ (out[Len(out)] = End \/ Raised(out[Len(out)]))

(* vcf_it_next is called. *)
Call ==
    /\ pc = "idle" /\ ~Done
    /\ pc' = "take"
    /\ UNCHANGED <<Kind, FailAt, pos, have, batch, next, readErr, out>>

(* The batch is used up: the end of the file, or fill_batch. *)
Refill ==
    /\ pc = "take" /\ next = Len(batch)
    /\ IF have
         THEN /\ pc' = "fill" /\ batch' = <<>> /\ next' = 0 /\ readErr' = FALSE
              /\ UNCHANGED out
         ELSE /\ out' = Append(out, End) /\ pc' = "idle"
              /\ UNCHANGED <<batch, next, readErr>>
    /\ UNCHANGED <<Kind, FailAt, pos, have>>

(* One pass of the read_lines loop: line moves into the batch and read_line
   reads the one after it, which fails after line FailAt. The failure
   longjmps out of read_line, leaving have_line as it was. *)
ReadLine ==
    /\ pc = "fill" /\ Len(batch) < Batch /\ have
    /\ batch' = Append(batch, pos)
    /\ IF pos = FailAt
         THEN /\ readErr' = TRUE /\ pc' = "take" /\ UNCHANGED <<pos, have>>
         ELSE /\ pos' = pos + 1 /\ have' = (pos + 1 <= N) /\ UNCHANGED <<readErr, pc>>
    /\ UNCHANGED <<Kind, FailAt, next, out>>

(* read_lines returns: the batch is full or the file ended. *)
FillEnd ==
    /\ pc = "fill" /\ (Len(batch) = Batch \/ ~have)
    /\ pc' = "take"
    /\ UNCHANGED <<Kind, FailAt, pos, have, batch, next, readErr, out>>

(* vcf_it_next takes record batch[next]. *)
Take ==
    /\ pc = "take" /\ next < Len(batch)
    /\ LET i == batch[next + 1]
           last == next + 1 = Len(batch)
       IN  /\ next' = next + 1
           /\ CASE Kind[i] = "bad" -> out' = Append(out, ParseErr(i)) /\ pc' = "idle"
                [] last /\ readErr -> out' = Append(out, ReadErr(i)) /\ pc' = "idle"
                [] Kind[i] = "excl" -> UNCHANGED <<out, pc>>
                [] OTHER -> out' = Append(out, Rec(i)) /\ pc' = "idle"
    /\ UNCHANGED <<Kind, FailAt, pos, have, batch, readErr>>

Stop == Done /\ pc = "idle" /\ UNCHANGED vars

Next == Call \/ Refill \/ ReadLine \/ FillEnd \/ Take \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Call \/ Refill \/ ReadLine \/ FillEnd \/ Take)

-----------------------------------------------------------------------------
(* Every record and error comes in the order, and with the content, of
   reading on demand, and nothing follows an error or the end. *)
OnDemandPrefix ==
    /\ Len(out) <= Len(Expected)
    /\ \A j \in 1..Len(out) : out[j] = Expected[j]

(* A read error is never lost: the batch it ends is not left without it. *)
ReadErrKept == readErr => FailAt \in {batch[j] : j \in 1..Len(batch)}

(* The caller sees every record and then the end or the first error. *)
Completes == <>(out = Expected)
=============================================================================
