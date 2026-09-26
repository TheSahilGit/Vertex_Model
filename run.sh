#!/usr/bin/env bash
# Launch the already-built ./bin/vertexmain in the background via nohup,
# so it keeps running after this shell/session exits -- same pattern the
# shellscripts/tcsh/ cluster scripts already use
# (`nohup ./bin/vertexmain > nohup.out &`), just as a standalone
# convenience wrapper. Does NOT build -- run `make`/`make fbounds`
# yourself first (see Makefile/README.md).
#
# Usage: ./run.sh
#
# Must be run with bin/vertexmain already built, and mesh/v_in.dat etc.
# and para_Simulation.dat already set up (see run_build_initial.sh /
# README.md) -- vertexmain opens its input/output files via relative
# paths from the CURRENT directory, which is why this script cd's to its
# own location (the repo root) first, regardless of where it's invoked
# from.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# Matched by exact process name (`pgrep -x`, i.e. /proc/<pid>/comm) first,
# THEN narrowed down by comparing /proc/<pid>/cwd to this repo root.
# Deliberately NOT `pgrep -f` (command-line substring match): that would
# false-positive on any unrelated process whose command line merely
# mentions the text "bin/vertexmain" (a comment, an echo, this very
# script being displayed) -- confirmed by testing. `-x vertexmain` only
# matches an actual running vertexmain binary, and cwd then distinguishes
# "running from HERE" from an unrelated instance in a different
# clone/copy of this repo.
REPO_ROOT="$(pwd)"
for pid in $(pgrep -x vertexmain 2>/dev/null || true); do
  if [ "$(readlink -f "/proc/$pid/cwd" 2>/dev/null)" = "$REPO_ROOT" ]; then
    echo "A vertexmain process from this repo is already running (PID $pid)" >&2
    echo "-- refusing to start a second one (it would clobber data/)." >&2
    exit 1
  fi
done

if [ ! -x ./bin/vertexmain ]; then
  echo "bin/vertexmain not found -- build it first (make or make fbounds)." >&2
  exit 1
fi

LOGFILE=nohup.out
nohup ./bin/vertexmain > "$LOGFILE" 2>&1 &
PID=$!
echo "vertexmain started in background (PID $PID), logging to $LOGFILE"
