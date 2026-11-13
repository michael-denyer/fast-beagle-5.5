------------------------------ MODULE BgenWriter ------------------------------
(***************************************************************************)
(* The BGEN output lifecycle: bgen_writer_open, bgen_writer_check_chrom    *)
(* (from window_writer_begin_window), bgen_writer_put, bgen_writer_close   *)
(* and write_sample_file in src/bgen/bgen_writer.c, in the order           *)
(* window_writer_close (src/main/window_writer.c) runs them: the VCF is    *)
(* closed, the .tbi written, then the BGEN members closed. The members are *)
(* created, completed and cleaned up by the actions of                     *)
(* tla/BgenCleanup.tla, which this spec instantiates rather than restates. *)
(* Init chooses the mode, each record's chromosome and kind, where windows *)
(* start, and which windows begin with an overlap record from the          *)
(* chromosome before (the comment at make_window in                        *)
(* src/vcf/sliding_window.c), whose chromosome window_writer_begin_window  *)
(* then checks. Chromosomes are 1..NChroms in file order; Autosomes are    *)
(* those plink2 accepts. A util_exit can come at any step (Fail); the      *)
(* checks and the zero-variant test exit by themselves. As in              *)
(* tla/BgenCleanup.tla, the main thread can keep going while exit() runs   *)
(* on another thread.                                                      *)
(***************************************************************************)
EXTENDS Integers, FiniteSets

CONSTANTS N, NChroms, Autosomes

ASSUME N \in Nat /\ NChroms \in Nat \ {0} /\ Autosomes \subseteq 1..NChroms

Members == {"bgen", "info", "sample"}

VARIABLES
    Mode,      \* "plink2" or "phased"
    Chrom,     \* per record: its chromosome, in file order
    Kind,      \* per record: what bgen_writer_encode made it
    WinStart,  \* the records after 1 that start a window
    Lag,       \* windows whose first target record is from the chromosome before
    pc,        \* the step the main thread takes next
    i,         \* the next record put passes to bgen_writer_put
    cur,       \* bw->chrom_index; 0 for -1
    checked,   \* chromosomes bgen_writer_check_chrom accepted
    nVar,      \* bw->n_variants
    written,   \* records whose bytes reached the .bgen
    hdrCount,  \* the variant count in the .bgen header; -1 until close writes it
    closedOK,  \* members fclose closed without error
    proc, created, completed, cleaned, removed

writer == <<Mode, Chrom, Kind, WinStart, Lag, pc, i, cur, checked, nVar, written, hdrCount, closedOK>>
cleanup == <<proc, created, completed, cleaned, removed>>
vars == <<writer, cleanup>>

C == INSTANCE BgenCleanup

Kinds == {"write", "filtered", "multiallelic"}

Sorted(f) == \A r \in 1..N - 1 : f[r] <= f[r + 1]

TypeOK ==
    /\ C!TypeOK
    /\ Mode \in {"plink2", "phased"}
    /\ pc \in {"open", "info", "begin", "put", "vcf", "tbi", "zero", "bgen", "infoclose",
               "sample", "sampleclose", "complete", "done"}
    /\ i \in 1..N + 1
    /\ cur \in 0..NChroms
    /\ checked \subseteq 1..NChroms
    /\ nVar \in 0..N
    /\ written \subseteq 1..N
    /\ hdrCount \in -1..N
    /\ closedOK \subseteq Members

Init ==
    /\ Mode \in {"plink2", "phased"}
    /\ Chrom \in {f \in [1..N -> 1..NChroms] : Sorted(f)}
    /\ Kind \in [1..N -> Kinds]
    /\ WinStart \in SUBSET (2..N)
    /\ Lag \in SUBSET {r \in WinStart : Chrom[r - 1] # Chrom[r]}
    /\ pc = "open"
    /\ i = 1
    /\ cur = 0
    /\ checked = {}
    /\ nVar = 0
    /\ written = {}
    /\ hdrCount = -1
    /\ closedOK = {}
    /\ C!Init

(* The process has not finished exit(). *)
Running == C!Live

(* A util_exit from the step the main thread is on. *)
Halt == C!Exit /\ UNCHANGED writer

(* pc moves on; nothing else changes. *)
Goto(p) == pc' = p /\ UNCHANGED <<Mode, Chrom, Kind, WinStart, Lag, i, cur, checked, nVar, written, hdrCount, closedOK>>

(* bgen_writer_open creates the .bgen, then the .info. *)
OpenBgen == Running /\ pc = "open" /\ C!Open("bgen") /\ Goto("info")

OpenInfo == Running /\ pc = "info" /\ C!Open("info") /\ Goto(IF N = 0 THEN "vcf" ELSE "begin")

(* window_writer_begin_window: bgen_writer_check_chrom on the chromosome of
   the window's first target record. *)
Begin ==
    /\ Running /\ pc = "begin"
    /\ LET c == IF i \in Lag THEN Chrom[i - 1] ELSE Chrom[i]
       IN  IF Mode = "plink2" /\ c \notin Autosomes
           THEN Halt
           ELSE /\ checked' = IF Mode = "plink2" THEN checked \cup {c} ELSE checked
                /\ pc' = "put"
                /\ UNCHANGED <<Mode, Chrom, Kind, WinStart, Lag, i, cur, nVar, written, hdrCount, closedOK>>
                /\ UNCHANGED cleanup

(* bgen_writer_put for record i: the chromosome check on a change of
   chromosome, for every kind, then the write. *)
Put ==
    /\ Running /\ pc = "put" /\ i <= N
    /\ LET c == Chrom[i]
           check == Mode = "plink2" /\ c # cur
       IN  IF check /\ c \notin Autosomes
           THEN Halt
           ELSE /\ cur' = IF check THEN c ELSE cur
                /\ checked' = IF check THEN checked \cup {c} ELSE checked
                /\ written' = IF Kind[i] = "write" THEN written \cup {i} ELSE written
                /\ nVar' = IF Kind[i] = "write" THEN nVar + 1 ELSE nVar
                /\ i' = i + 1
                /\ pc' = IF i = N THEN "vcf" ELSE IF i + 1 \in WinStart THEN "begin" ELSE "put"
                /\ UNCHANGED <<Mode, Chrom, Kind, WinStart, Lag, hdrCount, closedOK>>
                /\ UNCHANGED cleanup

(* window_writer_close: bgzf_close of the VCF, then vcf_index_write. *)
CloseVcf == Running /\ pc = "vcf" /\ Goto("tbi") /\ UNCHANGED cleanup

WriteTbi == Running /\ pc = "tbi" /\ Goto("zero") /\ UNCHANGED cleanup

(* bgen_writer_close: the zero-variant exit, then the variant count into the
   header. *)
Zero ==
    /\ Running /\ pc = "zero"
    /\ IF Mode = "plink2" /\ nVar = 0
       THEN Halt
       ELSE /\ hdrCount' = nVar
            /\ pc' = "bgen"
            /\ UNCHANGED <<Mode, Chrom, Kind, WinStart, Lag, i, cur, checked, nVar, written, closedOK>>
            /\ UNCHANGED cleanup

(* fclose of a member; a failing fclose is a Fail step. *)
Close(m, p) ==
    /\ closedOK' = closedOK \cup {m}
    /\ pc' = p
    /\ UNCHANGED <<Mode, Chrom, Kind, WinStart, Lag, i, cur, checked, nVar, written, hdrCount>>
    /\ UNCHANGED cleanup

CloseBgen == Running /\ pc = "bgen" /\ Close("bgen", "infoclose")

CloseInfo == Running /\ pc = "infoclose" /\ Close("info", "sample")

(* write_sample_file creates the .sample only now. *)
OpenSample == Running /\ pc = "sample" /\ C!Open("sample") /\ Goto("sampleclose")

CloseSample == Running /\ pc = "sampleclose" /\ Close("sample", "complete")

Complete == Running /\ pc = "complete" /\ C!Complete /\ Goto("done")

(* main returns 0. *)
Return == pc = "done" /\ C!Return /\ UNCHANGED writer

(* Any other util_exit: a failed write, open or fclose, an error on another
   thread, or one after window_writer_close returns (run_stats_close). *)
Fail == Running /\ Halt

Steps == OpenBgen \/ OpenInfo \/ Begin \/ Put \/ CloseVcf \/ WriteTbi \/ Zero \/ CloseBgen
         \/ CloseInfo \/ OpenSample \/ CloseSample \/ Complete \/ Return

Next == Steps \/ Fail \/ (C!Cleanup /\ UNCHANGED writer) \/ (C!End /\ UNCHANGED writer)
        \/ (C!Ended /\ UNCHANGED writer)

Spec == Init /\ [][Next]_vars /\ WF_vars(Steps) /\ WF_vars(C!Cleanup /\ UNCHANGED writer)
        /\ WF_vars(C!End /\ UNCHANGED writer)

-----------------------------------------------------------------------------
(* plink2 imports autosomes only, so no other chromosome reaches the .bgen. *)
WrittenAutosomes == Mode = "plink2" => \A r \in written : Chrom[r] \in Autosomes

(* In plink2 mode a record's chromosome was checked before its bytes were
   written. *)
CheckedBeforeWrite == Mode = "plink2" => \A r \in written : Chrom[r] \in checked

(* The header's variant count is the number of records in the file. *)
HeaderCountsWritten == hdrCount >= 0 => hdrCount = Cardinality(written)

(* bgen_files_complete runs only after the .bgen, .info and .sample members
   were all closed without error. *)
CompleteAfterClose == completed => closedOK = Members

(* plink2 mode with no record written exits before completing. *)
Plink2EmptyExits == Mode = "plink2" /\ written = {} => ~completed

(* A run that exits 0 leaves all three members, closed and counted. *)
SuccessKeepsAll ==
    proc = "exit0" => created \ removed = Members /\ closedOK = Members /\ hdrCount = Cardinality(written)

(* The run ends, with status 0 or 1. *)
Terminates == <>(proc \in {"exit0", "exit1"})
=============================================================================
