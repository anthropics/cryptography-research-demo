# Key-recovery attack on 13-round LEA-128

This is research code implementing a conditional differential-linear
key-recovery attack on LEA-128 reduced to 13 rounds (the full cipher has 24).
It recovers the 128-bit master key from chosen-plaintext queries to an
encryption-only oracle under a single fixed key. It is written to be read and
re-run, not to be production software.

## Target and attack model

- Cipher: LEA-128 as specified at WISA 2013 (standard key schedule and round
  function), reduced to 13 rounds. No other modification. The included
  `lea.h` reproduces the published LEA-128 test vector at the full 24 rounds.
- Oracle: chosen-plaintext encryption only. One fixed unknown key. No
  decryption queries, no related keys, no side channels.
- Cost at default settings: about 2^28.6 chosen-plaintext blocks (the attack
  prints its exact oracle usage at the end of every completed attempt,
  successful or not), plus a CPU search
  whose dominant part is a beam search over about 2^27 ciphertext pairs and a
  final brute force of 2^32 key-schedule completions per beam candidate.

## How the attack works

1. The plaintext difference (0x80000000, 0x80400000, 0x80400010, 0x80400014),
   together with six key-dependent bit conditions on the plaintexts, forces
   the round-1 output difference to be 0x80000000 in every word with
   probability 1. The six conditions are controlled by six chosen plaintext
   bits (the "combo", 0..63); which combo is correct depends on the key.
2. An all-MSB difference is free through modular addition, so that difference
   propagates two more rounds with probability 1.
3. The resulting 13-round differential-linear distinguisher has biases up to
   about 4e-2 on bits of the round-11 difference.
4. Bit 0 of a modular sum is linear, so bit j of several deep-state difference
   words can be computed from the ciphertext given only bits 0..j-1 of four
   32-bit words derived from the last round keys (A0 = rk[12][0],
   T12 = rk[12][1]^rk[12][2], T13 = rk[12][1]^rk[12][3],
   U = rk[12][1]^rk[11][0]). This triangular structure enables bit-by-bit
   recovery. The attack scores four such observables with complementary
   key-dependent dead bits, which keeps signal at essentially every bit
   position across keys.

The pipeline in `attack.c`:

- Phase 0: fetch two known plaintext/ciphertext blocks (used only to verify
  candidate keys; they are obtained from the same oracle).
- Phase 1: for each of the 64 combos, gather N1 pairs and run a small beam
  over the low 7 key bits; the correct combo separates clearly at the default
  N1 = 2^20.
- Phase 2: gather N2 pairs (default 2^27) for the best combo and run a
  bit-by-bit beam search (default width 64) over the low MAXBIT (default 24)
  bits of each of the four derived words. Each candidate is scored
  statelessly: the four-word subtraction chain is recomputed per ciphertext
  pair over the stored (rotated) ciphertext words, and the score correlates
  one delta bit per level across all pairs, accumulated cumulatively across
  levels.
- Phase 3: for each beam candidate in order, brute-force the remaining
  8 high bits of each word (2^32 combinations), expanding each guess through
  the inverted key schedule and testing a full 13-round encryption against a
  known block. A hit yields the master key directly.
- The recovered key is re-verified against a second known block, and in
  self-test mode additionally compared against the simulated oracle's key.

A failed run says so explicitly and exits nonzero: if no beam candidate
survives phase 3, the key was not recovered at those settings. `run.sh`
then escalates, at most twice, with FRESH data (a new seed) and a much wider
beam each time: width 64 at the first attempt, then 512 (with doubled
phase-1 data), then 1024 (with doubled phase-1 and phase-2 data). Width
costs time, not memory, so the rungs differ
mainly in how long they run. About 44% of random-key draws need at
least one escalation: some keys have runs of consecutive key-bit positions where
the distinguisher's signal is at the noise floor (the dead zones are a
property of the key, not of the data draw), and the beam holds the true
prefix through such a zone only if it is wide enough. Each attempt prints
its settings and memory estimate; if all rungs fail, run.sh reports every
setting tried and exits nonzero.

## Key-to-key difficulty variation

Not all keys are equally hard. The distinguisher's per-bit signal depends on
the key itself: on just over half of keys the signal stays strong enough at
every bit position for the default-width beam to hold the true prefix, but
some keys have runs of consecutive positions where the signal sits at the
noise floor.
These dead zones are a property of the key, not of the data draw - they
reproduce across independent draws. Through a dead zone the beam carries the
true prefix on inherited score alone, so recovery comes down to whether the
true prefix's rank among all candidates stays inside the beam width until
signal returns. A key whose dead zone spans many consecutive bit positions
can exceed any practical width (the candidate pool grows 16x per
zero-signal level); against such keys this package's honest output is its
explicit ATTACK FAILED report.

Two worked examples, taken from the keys that failed the full ladder during
earlier random-key validation batteries (distinct from the 50-key battery
reported under Resource requirements). Both are measured with the TRUTH
diagnostic (print-only; see Resource requirements) across widths and
independent data draws; every measurement below used N1=2^20 and N2=2^28
pairs, with SEED values 1, 2, 3 as the independent draws, e.g.:

```
K=e9c625310539c262559c9012dfca82a7
TRUTH=$K ./attack 1048576 268435456 1024 32 24 1 $K
```

- Key e9c625310539c262559c9012dfca82a7: the true prefix is dropped at bit 12
  with rank near 8,200 - 8,193 and 8,196 on two independent draws at width
  1024, and 8,193 again at width 2048 on a controlled same-data comparison.
  The rank is essentially width-independent: widening the beam added only
  candidates that score below the true prefix (when the pool doubled on the
  same data, the set ranked above it was unchanged), so every fixed width
  below roughly 8,200 fails identically, and the measured rank says exactly
  what width would have been needed on those draws.
- Key b58c1fc2d49053533463a37dc3b21afa: worse, and noisier across draws. On
  three independent draws it lost the true prefix at bit 10 with rank 2,113
  of 32,768 (width 2048); at bit 11 with rank 29,961 of 32,768 (width 2048);
  and at bit 11 with rank 27,936 of 65,536 (width 4096, after surviving bit
  10 at rank 1,763). At the bit-11 level the true prefix's rank is deep in
  the candidate pool (the 91st and 43rd percentile on the two draws that
  reached it) - consistent with a uniformly random position, i.e.
  statistically indistinguishable from noise - so no practical fixed width
  rescues this key.

The width-versus-runtime trade, for anyone who wants to push further: memory
is width-independent (32 bytes per ciphertext pair, regardless of width),
and phase-2 runtime scales roughly linearly with width. Our width-2048 runs
on the two keys above took about 14 hours each; a width-32,768 attempt on a
b58c1fc2-class key is 16x that width, so by the same linear scaling it
extrapolates to roughly 8-10 days on the same machine - and still carries
no guarantee, since the measured bit-11 ranks (27,936 and 29,961) sit close
to that width on the two draws that reached the level. One practical note
for anyone trying it: phase 3 tries at most PH3MAX beam candidates (an
environment variable, default 1024), so a much wider beam should raise
PH3MAX to match - otherwise a true prefix held deep in the final beam is
never tried. Keys of this class need a qualitatively different statistic,
not more of this one.

## Building and running

```
make                       # cc -O3, needs only libc, pthreads, libm
make check                 # verifies lea.h against the published LEA-128
                           # test vector (prints KAT PASS)
./run.sh random            # self-test: attack a random key via the built-in
                           # local oracle (no network, no dependencies)
```

Memory up front: the default settings need about 4 GiB of RAM, and a hard
key that drives run.sh through its full escalation ladder needs about
8 GiB at the final rung (details under Resource requirements). On machines
without enough RAM (and without swap) the driver aborts with a clear message
instead of thrashing; with swap enabled, some thrashing can precede the abort.

To run the attack over HTTP against the included oracle server (requires
Python 3 with numpy — `pip install -r requirements.txt`):

```
python3 oracle_server.py &         # prints the secret key, serves on :8666
./run.sh                           # attacks it; compare the printed KEY
```

To point it at another oracle implementing the same API
(POST /encrypt with body {"pt":"<hex>"} returning {"ct":"<hex>"},
up to 1 MiB per request), set ORACLE_HOST (an IPv4 literal) and
ORACLE_PORT. The client expects Content-Length on responses.

Direct invocation (all parameters on the command line):

```
./attack [N1 [N2 [BEAM [THREADS [MAXBIT [SEED [SIMKEY]]]]]]]
```

Defaults: N1=2^20, N2=2^27, BEAM=64, THREADS=8, MAXBIT=24. THREADS accepts
1..1024; run.sh overrides the default with the machine's core count (capped
at 1024). With a SIMKEY
argument (32 hex chars) the oracle is simulated in-process with that key;
without it the HTTP oracle is used. All chosen plaintexts are generated from
SEED by a fixed PRNG, so runs are reproducible given (key, settings, SEED),
independent of the thread count: results are identical at any THREADS value,
which only changes wall-clock time (verified: the output is byte-identical,
apart from the wall-clock lines, across a wide range of THREADS values).

Exit codes (attack, and run.sh after its ladder): 0 key recovered and
verified; 1 honest failure - the key was not found at the attempted
settings; 2 bad parameters; 3 resource allocation failure (memory or
threads); 4 oracle communication failure; 128+N killed by signal N
(137 means the kernel's out-of-memory killer).

An experimental adaptive-width beam ships in the code but is OFF by default
(environment ADAPTIVE=1 enables it; DEADCAP, FLATX, GAPTHR, FLATTHR tune
it). Every number in this README was measured with it off - the validated
configuration is the classic fixed-width beam.

## Resource requirements

- Machine: any modern multi-core machine. The numbers below are from runs
  at 32 threads; the result never depends on the thread count (see the
  determinism note under Building and running), only the runtime does.
- Memory: about 4 GiB at default settings, about 8 GiB at the escalation
  ladder's final rung. The attack prints its estimate at phase 2 start and
  aborts cleanly if the machine cannot fit it.
- Runtime: a typical key recovers in about 5 minutes. Keys that need one
  escalation usually finish in about half an hour; the hardest runs,
  through the full three-rung ladder, take roughly 2.5-4 hours. Phase 2
  dominates and scales roughly as N2 x BEAM / threads. Build with the
  Makefile's -O3: the scoring loop relies on compiler autovectorization,
  and weaker optimization costs several times more in phase 2.
- Validation: 50 fresh random keys (generated up front into a key
  manifest) attacked end to end through run.sh (self-test mode):
  47 recovered (94%); about 44% of the keys needed at least one
  escalation. `validate.sh` reproduces this battery.
- Self-test note: `validate.sh`'s quick check uses a fixed known-good key
  on purpose. A fresh RANDOM key fails the FIRST rung on about 44%
  of draws (the dead-zone class above) - that is the expected behavior the
  ladder exists for, not a broken build - so random-key testing must allow
  the full ladder to run.
- Diagnostics: in simulation mode, setting TRUTH=<the same 32 hex chars as
  SIMKEY> makes the beam print the true candidate's rank at every bit level
  (print-only; never affects the search). This shows exactly where and how
  a failing key loses the true prefix.

## Files

- `attack.c` - the whole attack (phases 0-3), pthreads
- `lea.h` - compact LEA-128 implementation written from the specification (variable rounds)
- `oracle.h` - HTTP oracle client and the built-in simulation mode
- `oracle_server.py` - local oracle server (numpy)
- `Makefile` - build (plain cc -O3; see Resource requirements on flags)
- `run.sh` - one-command driver with the bounded escalation ladder
- `validate.sh` - multi-key validation battery
- `kat.c` - known-answer test for lea.h (`make check`)
