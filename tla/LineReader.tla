------------------------------ MODULE LineReader ------------------------------
(***************************************************************************)
(* The line splitting of src/blbutil/line_reader.c: fill,                  *)
(* line_reader_next_raw and line_reader_next, whose UTF-8 sanitising is    *)
(* not modeled. fill hands out the input one chunk at a time: an fread of  *)
(* a plain file, a gzip inflate, or one bgzf block. htslib never returns   *)
(* an empty block before the end (bgzf_read_block skips them), so every    *)
(* chunk before the end is non-empty. BgzipLines is r->bgzip_lines: TRUE   *)
(* splits as BGZipIt does, FALSE as BufferedReader.readLine does. Init     *)
(* chooses the input, every string over {a, CR, LF} of at most K bytes,    *)
(* and its split into chunks. The caller calls line_reader_next until it   *)
(* has seen the end twice. Each step is one pass of the loop in            *)
(* line_reader_next_raw, with the skip_lf test as a step of its own.       *)
(***************************************************************************)
EXTENDS Naturals, Sequences

CONSTANTS K, BgzipLines

CR == "CR"
LF == "LF"
Alpha == {"a", CR, LF}

VARIABLES
    bytes,     \* the decompressed input
    cuts,      \* a chunk ends after each of these bytes, and after the last
    gp,        \* r->buf + r->pos, as an index into bytes
    chunkEnd,  \* the last byte of the chunk in r->buf; r->pos == r->len is gp > chunkEnd
    eof,       \* r->eof
    skipLf,    \* r->skip_lf
    eofReads,  \* reads of the source that returned no bytes
    line,      \* the kstring being built
    any,       \* any in line_reader_next_raw
    pc,        \* "idle" or "loop" (inside line_reader_next_raw)
    out,       \* what line_reader_next returned, in order
    expected   \* Expected for this input, fixed by Init

vars == <<bytes, cuts, gp, chunkEnd, eof, skipLf, eofReads, line, any, pc, out, expected>>

Line(s) == [k |-> "line", s |-> s]
End == [k |-> "end", s |-> <<>>]

Drop(s, n) == SubSeq(s, n + 1, Len(s))

(* The first index in p..e whose byte is in S, or e + 1 for none. *)
FirstIn(p, e, S) ==
    IF \E q \in p..e : bytes[q] \in S
    THEN CHOOSE q \in p..e : bytes[q] \in S /\ \A r \in p..(q - 1) : bytes[r] \notin S
    ELSE e + 1

(* BufferedReader.readLine: a line ends at LF, CR or CR LF, and text after
   the last terminator is a line. *)
RECURSIVE ReaderSplit(_)
ReaderSplit(s) ==
    IF s = <<>> THEN <<>>
    ELSE LET t == CHOOSE q \in 1..(Len(s) + 1) :
                      /\ (q <= Len(s) => s[q] \in {CR, LF})
                      /\ \A r \in 1..(q - 1) : s[r] \notin {CR, LF}
         IN  IF t > Len(s) THEN <<s>>
             ELSE <<SubSeq(s, 1, t - 1)>> \o
                  ReaderSplit(Drop(s, IF s[t] = CR /\ t < Len(s) /\ s[t + 1] = LF THEN t + 1 ELSE t))

(* BGZipIt: a line ends at LF and loses one trailing CR, and text after the
   last LF is dropped. *)
StripCr(s) == IF s # <<>> /\ s[Len(s)] = CR THEN SubSeq(s, 1, Len(s) - 1) ELSE s

RECURSIVE BgzipSplit(_)
BgzipSplit(s) ==
    IF \A q \in 1..Len(s) : s[q] # LF THEN <<>>
    ELSE LET t == CHOOSE q \in 1..Len(s) : s[q] = LF /\ \A r \in 1..(q - 1) : s[r] # LF
         IN  <<StripCr(SubSeq(s, 1, t - 1))>> \o BgzipSplit(Drop(s, t))

RefLines == IF BgzipLines THEN BgzipSplit(bytes) ELSE ReaderSplit(bytes)

Expected == [j \in 1..Len(RefLines) |-> Line(RefLines[j])] \o <<End, End>>

(* The last byte of the chunk fill returns when the next byte is p. *)
ChunkEnd(p) ==
    IF \E c \in cuts : c >= p
    THEN CHOOSE c \in cuts : c >= p /\ \A d \in cuts : d >= p => c <= d
    ELSE Len(bytes)

TypeOK ==
    /\ bytes \in Seq(Alpha) /\ Len(bytes) <= K
    /\ cuts \subseteq 1..K
    /\ gp \in 1..(Len(bytes) + 1)
    /\ chunkEnd \in 0..Len(bytes)
    /\ eof \in BOOLEAN
    /\ skipLf \in BOOLEAN
    /\ eofReads \in Nat
    /\ line \in Seq(Alpha)
    /\ any \in BOOLEAN
    /\ pc \in {"idle", "loop"}

Init ==
    /\ bytes \in UNION {[1..n -> Alpha] : n \in 0..K}
    /\ cuts \in SUBSET (1..(Len(bytes) - 1))
    /\ gp = 1
    /\ chunkEnd = 0
    /\ eof = FALSE
    /\ skipLf = FALSE
    /\ eofReads = 0
    /\ line = <<>>
    /\ any = FALSE
    /\ pc = "idle"
    /\ out = <<>>
    /\ expected = Expected

Done == Len(out) = Len(expected)

(* line_reader_next is called. *)
Call ==
    /\ pc = "idle" /\ ~Done
    /\ pc' = "loop" /\ line' = <<>> /\ any' = FALSE
    /\ UNCHANGED <<bytes, cuts, gp, chunkEnd, eof, skipLf, eofReads, out, expected>>

(* r->pos == r->len: fill gives the next chunk, or the source has ended and
   the call returns the text it holds (BufferedReader) or false. Once eof
   is set, fill returns false without reading. *)
Fill ==
    /\ pc = "loop" /\ gp > chunkEnd
    /\ IF ~eof /\ gp <= Len(bytes)
         THEN /\ chunkEnd' = ChunkEnd(gp)
              /\ UNCHANGED <<eof, eofReads, pc, out>>
         ELSE /\ eof' = TRUE
              /\ eofReads' = IF eof THEN eofReads ELSE eofReads + 1
              /\ out' = Append(out, IF any /\ ~BgzipLines THEN Line(line) ELSE End)
              /\ pc' = "idle"
              /\ UNCHANGED chunkEnd
    /\ UNCHANGED <<bytes, cuts, gp, skipLf, line, any, expected>>

(* The last line ended at a CR: a LF that follows it, even as the first
   byte of a new chunk, is the rest of that terminator. *)
SkipLf ==
    /\ pc = "loop" /\ gp <= chunkEnd /\ skipLf
    /\ skipLf' = FALSE
    /\ gp' = IF bytes[gp] = LF THEN gp + 1 ELSE gp
    /\ UNCHANGED <<bytes, cuts, chunkEnd, eof, eofReads, line, any, pc, out, expected>>

(* The rest of the chunk up to the first LF, or for BufferedReader the
   first CR before it, joins the line. A terminator ends the call. *)
Scan ==
    /\ pc = "loop" /\ gp <= chunkEnd /\ ~skipLf
    /\ LET lf == FirstIn(gp, chunkEnd, {LF})
           t == IF BgzipLines THEN lf ELSE FirstIn(gp, lf - 1, {CR})
           l1 == line \o SubSeq(bytes, gp, t - 1)
       IN  /\ any' = TRUE
           /\ IF t > chunkEnd
                THEN /\ line' = l1 /\ gp' = chunkEnd + 1
                     /\ UNCHANGED <<skipLf, pc, out>>
                ELSE /\ skipLf' = (~BgzipLines /\ bytes[t] = CR)
                     /\ gp' = t + 1
                     /\ line' = IF BgzipLines THEN StripCr(l1) ELSE l1
                     /\ out' = Append(out, Line(line'))
                     /\ pc' = "idle"
    /\ UNCHANGED <<bytes, cuts, chunkEnd, eof, eofReads, expected>>

Stop == Done /\ pc = "idle" /\ UNCHANGED vars

Next == Call \/ Fill \/ SkipLf \/ Scan \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Call \/ Fill \/ SkipLf \/ Scan)

-----------------------------------------------------------------------------
(* Every result is the reference splitter's next line, and after the last
   line every call returns false, for every split into chunks. *)
OutPrefix ==
    /\ Len(out) <= Len(expected)
    /\ \A j \in 1..Len(out) : out[j] = expected[j]

(* The source is read at most once after it has ended. *)
EofReadOnce == eofReads <= 1

(* The caller sees every line and then false. *)
Completes == <>(out = expected)
=============================================================================
