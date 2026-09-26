#!/usr/bin/env bash
# Run the already-built ./bin/generate_initial_mesh, overwriting
# mesh/v_in.dat, mesh/inn_in.dat, mesh/num_in.dat, mesh/para_MeshDims.dat
# from para_MeshGen.dat (edit that first if you want different mesh
# parameters -- see PARAMETERS.md). Does NOT build -- run
# `make initial_config` yourself first (see Makefile/README.md).
#
# Runs in the FOREGROUND, not backgrounded via nohup like run.sh: this is
# a one-shot, seconds-scale preprocessing step (unlike vertexmain, which
# can run for hours/days), and run.sh depends on these mesh/ files being
# completely written before it starts -- backgrounding this one would risk
# launching vertexmain against a half-written mesh. Ask if you'd rather
# have this backgrounded too.
#
# Usage: ./run_build_initial.sh
#
# Must be run from anywhere -- it cd's to its own location (the repo
# root) first, since generate_initial_mesh opens para_MeshGen.dat and
# writes mesh/... via relative paths from the CURRENT directory.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# See run.sh's matching comment -- `pgrep -x vertexmain` (exact process
# name) rather than `-f` (command-line substring, which false-positives
# on any unrelated process whose command line merely mentions the text
# "bin/vertexmain" -- confirmed by testing), then cwd distinguishes
# "running from HERE" from an unrelated instance in a different
# clone/copy of this repo.
REPO_ROOT="$(pwd)"
for pid in $(pgrep -x vertexmain 2>/dev/null || true); do
  if [ "$(readlink -f "/proc/$pid/cwd" 2>/dev/null)" = "$REPO_ROOT" ]; then
    echo "A vertexmain process from this repo is currently running (PID $pid)" >&2
    echo "-- refusing to regenerate the mesh out from under it. Stop it" >&2
    echo "first if you really want a fresh mesh, then re-run this script." >&2
    exit 1
  fi
done

if [ ! -x ./bin/generate_initial_mesh ]; then
  echo "bin/generate_initial_mesh not found -- build it first (make initial_config)." >&2
  exit 1
fi

./bin/generate_initial_mesh
