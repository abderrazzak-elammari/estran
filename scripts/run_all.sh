#!/bin/sh
# Every case behind the README, then every figure. About an hour on 2 cores,
# most of it the fine island mesh.
#   NP=8 ./scripts/run_all.sh
set -e
cd "$(dirname "$0")/.."
NP=${NP:-4}
MPIRUN=${MPIRUN:-mpirun}

make -s estran meshes
mkdir -p results
run() { echo "== $1"; $MPIRUN -np "$NP" ./estran "$1" | tail -1; }

run cases/lake.nml
run cases/bump.nml
for k in 1 2 3 4; do
  sed -e "s|dambreak_3.msh|dambreak_$k.msh|" -e "s|results/dambreak'|results/dambreak_$k'|" \
    cases/dambreak.nml > results/dambreak_$k.nml
  run results/dambreak_$k.nml
done
run cases/dambreak_dry.nml
run cases/thacker.nml
run cases/island_a.nml
run cases/island_c.nml

# run-up sensitivity for case A: half the cell size on the island, friction
sed -e "s|meshes/island.msh|meshes/island_fine.msh|" -e "s|results/island_a'|results/island_a_fine'|" \
    -e "s|dt_out    = 0.5|dt_out    = 20.0|" cases/island_a.nml > results/island_a_fine.nml
run results/island_a_fine.nml
sed -e "s|results/island_a'|results/island_a_manning'|" -e "s|dt_out    = 0.5|dt_out    = 20.0|" \
    -e "s|^/$|  manning = 0.012\n/|" cases/island_a.nml > results/island_a_manning.nml
run results/island_a_manning.nml

for s in plot_bump plot_dambreak plot_thacker plot_island plot_partition; do
  python3 scripts/$s.py
done
