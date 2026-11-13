--------------------------- MODULE PhasingSchedule ---------------------------
(***************************************************************************)
(* The iteration schedule of one window: phase_window and                  *)
(* phase_stage1_variants in src/main/main.c, lr_threshold,                 *)
(* phase_data_init, phase_data_increment_it and                            *)
(* phase_data_advance_to_first_phasing_it in src/phase/phase_data.c, and   *)
(* the EM gating and PBWT direction of phase_ls_run_stage1 in              *)
(* src/phase/phase_ls.c. par.c allows burnin >= 1 and iterations >= 1, as  *)
(* Java's Par does. The threshold is abstracted to the branch of           *)
(* lr_threshold that computed it: "inf" (it < burnin), "one" (the last     *)
(* phasing iteration), "div0" (iterations = 1 and it > burnin, where the   *)
(* exponent is -1 / 0), "past" (it beyond the last phasing iteration, a    *)
(* negative exponent) and "finite" (the geometric fall). Init chooses      *)
(* whether the input is already phased and whether every marker is a       *)
(* stage-1 marker; each run's swap rate is chosen when it ends: "low"      *)
(* (at most 1%), "high", or "nan" (0 / 0 when no sample had an unphased    *)
(* het, for which <= is false), or "unread" once the run reached burnin,   *)
(* where && never reaches the comparison. Each step is one run of          *)
(* phase_ls_run_stage1 with the increment after it, the swap-rate test     *)
(* that follows, or one path choice in phase_window.                       *)
(***************************************************************************)
EXTENDS Integers, Sequences

CONSTANTS Burnin, Iterations

ASSUME Burnin \in Nat \ {0} /\ Iterations \in Nat \ {0}

VARIABLES
    phased,     \* targ_is_phased(w)
    allStage1,  \* fpd->n_stage1 == fpd->n_markers
    it,         \* pd->it
    lr,         \* the branch of lr_threshold that set pd->lr_threshold
    pc,         \* "start", "loop", "check", "stage2", "write" or "done"
    runs,       \* each phase_ls_run_stage1: its it, threshold, EM step, direction
    rates,      \* each run's swap rate
    traced,     \* the threshold phase_data_trace_iteration prints after each run
    stage2It,   \* pd->it when phase_ls_run_stage2 starts; -1 if it does not
    path        \* "none", "phased", "stage1" or "stage2"

vars == <<phased, allStage1, it, lr, pc, runs, rates, traced, stage2It, path>>

NIts == Burnin + Iterations
Regimes == {"inf", "one", "div0", "past", "finite"}
Rates == {"low", "high", "nan"}

(* lr_threshold(p, it), by the branch it returns from. *)
Regime(i) ==
    LET nItsM1 == Iterations - 1
    IN  IF i < Burnin THEN "inf"
        ELSE IF i = nItsM1 + Burnin THEN "one"
        ELSE IF nItsM1 = 0 THEN "div0"
        ELSE IF i - Burnin > nItsM1 THEN "past"
        ELSE "finite"

(* swap <= max_burnin_swap_rate *)
AtMostOnePercent(r) == r = "low"

(* while (pd->it < n_its) *)
LoopCond == it < NIts

Its == [i \in 1 .. Len(runs) |-> runs[i].it]

(* The runs a window must make when burn-in ends after k runs: burn-in
   iterations 0 .. k-1, then every phasing iteration once, in order. *)
Schedule(k) == [i \in 1 .. k |-> i - 1] \o [i \in 1 .. Iterations |-> Burnin + i - 1]

IsPrefix(s, t) == Len(s) <= Len(t) /\ \A i \in 1 .. Len(s) : s[i] = t[i]

TypeOK ==
    /\ phased \in BOOLEAN /\ allStage1 \in BOOLEAN
    /\ it \in 0 .. NIts
    /\ lr \in Regimes
    /\ pc \in {"start", "loop", "check", "stage2", "write", "done"}
    /\ runs \in Seq([it : 0 .. NIts, lr : Regimes, em : {"init", "update", "none"}, bwd : BOOLEAN])
    /\ rates \in Seq(Rates \cup {"unread"})
    /\ traced \in Seq(Regimes)
    /\ stage2It \in -1 .. NIts
    /\ path \in {"none", "phased", "stage1", "stage2"}

Init ==
    /\ phased \in BOOLEAN /\ allStage1 \in BOOLEAN
    /\ it = 0 /\ lr = Regime(0)
    /\ pc = "start"
    /\ runs = <<>> /\ rates = <<>> /\ traced = <<>>
    /\ stage2It = -1 /\ path = "none"

(* phase_window: the input is already phased, or stage 1 runs. *)
Begin ==
    /\ pc = "start"
    /\ IF phased THEN pc' = "write" /\ path' = "phased" ELSE pc' = "loop" /\ UNCHANGED path
    /\ UNCHANGED <<phased, allStage1, it, lr, runs, rates, traced, stage2It>>

(* One pass of the while loop: phase_ls_run_stage1 at it, with EM
   initialised at it 0 and updated during burn-in, and the backward PBWT on
   even iterations; then phase_data_increment_it. *)
Run ==
    /\ pc = "loop" /\ LoopCond
    /\ runs' = Append(runs, [it |-> it, lr |-> lr,
                             em |-> IF it = 0 THEN "init" ELSE IF it < Burnin THEN "update" ELSE "none",
                             bwd |-> it % 2 = 0])
    /\ it' = it + 1 /\ lr' = Regime(it + 1)
    /\ pc' = "check"
    /\ UNCHANGED <<phased, allStage1, rates, traced, stage2It, path>>

(* The swap rate of the run, the trace, and the early end of burn-in. *)
Check ==
    /\ pc = "check"
    /\ \E r \in IF it < Burnin THEN Rates ELSE {"unread"} :
         /\ rates' = Append(rates, r)
         /\ traced' = Append(traced, lr)
         /\ IF it < Burnin /\ AtMostOnePercent(r)
              THEN it' = Burnin /\ lr' = Regime(Burnin)
              ELSE UNCHANGED <<it, lr>>
    /\ pc' = "loop"
    /\ UNCHANGED <<phased, allStage1, runs, stage2It, path>>

(* The loop ends; phase_window picks the stage-1 or the stage-2 haplotypes. *)
LoopEnd ==
    /\ pc = "loop" /\ ~LoopCond
    /\ IF allStage1 THEN pc' = "write" /\ path' = "stage1" ELSE pc' = "stage2" /\ path' = "stage2"
    /\ UNCHANGED <<phased, allStage1, it, lr, runs, rates, traced, stage2It>>

(* phase_ls_run_stage2, whose pbwt_ibs_data_init reads pd->it. *)
Stage2 ==
    /\ pc = "stage2"
    /\ stage2It' = it /\ pc' = "write"
    /\ UNCHANGED <<phased, allStage1, it, lr, runs, rates, traced, path>>

(* The window is written. *)
Write ==
    /\ pc = "write" /\ pc' = "done"
    /\ UNCHANGED <<phased, allStage1, it, lr, runs, rates, traced, stage2It, path>>

Stop == pc = "done" /\ UNCHANGED vars

Next == Begin \/ Run \/ Check \/ LoopEnd \/ Stage2 \/ Write \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Begin \/ Run \/ Check \/ LoopEnd \/ Stage2 \/ Write)

-----------------------------------------------------------------------------
(* At most burnin + iterations runs. *)
RunBound == Len(runs) <= NIts

(* The runs so far are a prefix of a schedule, and a finished stage 1 ran a
   whole one: the burn-in runs come first, and every phasing iteration runs
   once, in order. *)
FollowsSchedule ==
    /\ \E k \in 1 .. Burnin : IsPrefix(Its, Schedule(k))
    /\ (pc \in {"stage2", "write", "done"} /\ ~phased) => \E k \in 1 .. Burnin : Its = Schedule(k)

(* Burn-in leaves early only after a rate of at most 1%, and always does. *)
EarlyExitOnLowRate ==
    \A i \in 1 .. Len(runs) - 1 :
        (runs[i + 1].it # runs[i].it + 1 <=> (rates[i] = "low" /\ runs[i].it + 1 < Burnin))

(* Every run uses a threshold its iteration should have: infinite in burn-in,
   1.0 on the last phasing iteration, and the geometric fall between; never
   one computed by division by zero or past the end. *)
RunThresholds ==
    \A i \in 1 .. Len(runs) :
        runs[i].lr = IF runs[i].it < Burnin THEN "inf"
                     ELSE IF runs[i].it = NIts - 1 THEN "one"
                     ELSE "finite"

(* The last run of stage 1 uses the 1.0 threshold. *)
LastRunOne == (pc = "done" /\ ~phased) => runs[Len(runs)].lr = "one"

(* EM initialises the parameters in the first run and updates them only in
   the other burn-in runs. *)
EmInBurnin ==
    \A i \in 1 .. Len(runs) :
        /\ (runs[i].em = "init") <=> (i = 1)
        /\ (runs[i].em = "update") <=> (i > 1 /\ runs[i].it < Burnin)

(* Stage 2 starts after the last increment, where pbwt_ibs_data_init takes
   the stage-2 candidate count, and only for a window with markers outside
   stage 1. *)
Stage2AfterStage1 ==
    /\ stage2It \in {-1, NIts}
    /\ path = "stage2" => ~phased /\ ~allStage1
    /\ pc = "done" => (stage2It = NIts <=> ~phased /\ ~allStage1)

(* A phased input runs nothing. *)
PhasedRunsNothing == phased => runs = <<>>

Terminates == <>(pc = "done")
=============================================================================
