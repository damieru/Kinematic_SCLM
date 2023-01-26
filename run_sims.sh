# ==== Calculates the velocity Fields ===
#cd ../../Velocity\ Fields
#./comp_run_parallel

# ====== Runs Kinematic Bulk Model ======
#cd ../../bulk_models/kinematic\ bulk\ ice/
#./compile

# == Runs Superdroplet Kinematic model ==
#cd ../../sd_models/kinematic\ SD\ ice/
./compile.sh
export OMP_NUM_THREADS=20;
#./ksdh5 < ice_D0.0_i0.010.in
#./ksdh5 < ice_D0.0_i0.050.in
#./ksdh5 < ice_D0.0_i0.100.in
./ksdh5 < Input/input.in

