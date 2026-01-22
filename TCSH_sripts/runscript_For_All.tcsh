#!/usr/bin/tcsh

# -----------------------------
# Base directory
# -----------------------------
set BASEDIR = "Vertex_Model"

# -----------------------------
# Types
# -----------------------------
set types = ( \
    "Run_Apolar_cell_motility" \
    "Run_Polar_cell_motility" \
    "Run_Fluctuating_contractility" \
    "Run_Mechano_chemical_regulation" \
)

# -----------------------------
# Parameter values
# -----------------------------

## Keep the order intact. 
#  values1 --> Apolar_cell_motility, etas
#  values2 --> Polar_cell_motility, vo
#  values3 --> Run_Fluctuating_contractility, active_contr_strength
#  values4 --> Run_Mechano_chemical_regulation, coupling_noise_strength
##

# 1st set

#set values1 = ("0.02d0" "0.03d0" "0.04d0")
#set values2 = ("1.0d-1" "2.0d-1" "2.5d-1")
#set values3 = ("6.0d-1" "7.0d-1" "8.0d-1")
#set values4 = ("1.0d0" "1.1d0" "1.2d0")

## 2nd set

#set values1 = ("0.045" "0.05")
#set values2 = ("3.0d-1" "3.5d-1")


# -----------------------------
# Loop over types
# -----------------------------
@ itype = 1
foreach type ($types)

    echo "======================================"
    echo "Processing type: $type"
    echo "======================================"


     # Reset variables to avoid carry-over
      unset vals
      unset param_name

      # --------------------------------------------------
      # Select parameter set for this type (SAFE)
      # --------------------------------------------------
      if ($itype == 1 && $?values1) then
          set vals = ($values1)
          set param_name = "Apolar_cell_motility"

      else if ($itype == 2 && $?values2) then
          set vals = ($values2)
          set param_name = "vo"

      else if ($itype == 3 && $?values3) then
          set vals = ($values3)
          set param_name = "active_contr_strength"

      else if ($itype == 4 && $?values4) then
          set vals = ($values4)
          set param_name = "coupling_noise_strength"
      endif

      # --------------------------------------------------
      # Skip if no values defined
      # --------------------------------------------------
      if (! $?vals || $#vals == 0) then
          echo "No values defined for $type — skipping"
          @ itype++
          continue
      endif



    # -----------------------------
    # Loop over values
    # -----------------------------
    foreach val ($vals)

        set newDir = "${type}/${param_name}_${val}"
        echo "Creating run: $newDir"

        mkdir -p $type
        rsync -av --exclude 'nrun2_*' --exclude '.git*' $BASEDIR/ $newDir

        cd $newDir

        ## Modify parameter

        sed -i "/${param_name}/s/[0-9.eEd+-]\+/${val}/" para1_in.dat

        if ($type == "Run_Polar_cell_motility") then
            sed -i "/if_ABP/s/.false./.true./" para1_in.dat
        
        else if ($type == "Run_Fluctuating_contractility") then
            sed -i "/if_active_contractility/s/.false./.true./" para1_in.dat
        
        else if ($type == "Run_Mechano_chemical_regulation") then
            sed -i "/if_RhoROCK/s/.false./.true./" para1_in.dat
            sed -i "/if_RK4/s/.false./.true./" para1_in.dat
            sed -i "/if_coupling_noise/s/.false./.true./" para1_in.dat
        endif


        ## Run simulation

        sh clear.sh
        sh compile.sh
        nohup ./vertexmain.exe > nohup.out &

        cd ../../
        echo "Done: $newDir"

    end

    @ itype++
end

echo "All simulations created successfully."

