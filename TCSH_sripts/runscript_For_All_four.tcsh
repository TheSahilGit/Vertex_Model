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

set values1   = ("0.04d0")
set values2   = ("0.4d0")

# Two-parameter sweep
set values3   = ("0.14d0")
set values3_3 = ("0.8d0")
set beta_val = "0.08d0"

set values4   = ("1.0d0")


# -----------------------------
# Loop over types
# -----------------------------
@ itype = 1
foreach type ($types)

    echo "======================================"
    echo "Processing type: $type"
    echo "======================================"

    unset vals
    unset param_name

    if ($itype == 1 && $?values1) then
        set vals = ($values1)
        set param_name = "Apolar_motility_strength"

    else if ($itype == 2 && $?values2) then
        set vals = ($values2)
        set param_name = "polar_motility_strength"

    else if ($itype == 3 && $?values3) then
        set vals = ($values3)
        set vals2 = ($values3_3)
        set param_name  = "active_contr_strength"
        set param_name2 = "tau_contr"

    else if ($itype == 4 && $?values4) then
        set vals = ($values4)
        set param_name = "coupling_noise_strength"
    endif


    if (! $?vals || $#vals == 0) then
        echo "No values defined for $type — skipping"
        @ itype++
        continue
    endif


    # =====================================================
    # TWO PARAMETER SWEEP
    # =====================================================

    if ($type == "Run_Fluctuating_contractility") then

        foreach val ($vals)
            foreach tau ($vals2)

                set newDir = "${type}/${param_name}_${val}_tau_${tau}"
                echo "Creating run: $newDir"

                if (-d $newDir) then
                    echo "Already exists — skipping"
                    continue
                endif

                mkdir -p $type
                rsync -avz --exclude 'nrun2_*' --exclude '.git*' $BASEDIR/ $newDir

                cd $newDir

                # Update parameters
                sed -i "/${param_name}/s/[0-9.eEd+-]\+/${val}/" para1_in.dat
                sed -i "/${param_name2}/s/[0-9.eEd+-]\+/${tau}/" para1_in.dat
                sed -i "/beta/s/[0-9.eEd+-]\+/${beta_val}/" para1_in.dat

                sed -i "/if_active_contractility/s/.false./.true./" para1_in.dat


                # -----------------------------
                # RUN SIMULATION
                # -----------------------------

                sh clear.sh
                sh compile.sh
                
                echo "Launching simulation..."
                nohup ./vertexmain.exe > nohup.out &

                cd ../../
                echo "Done + Launched: $newDir"

            end
        end


    # =====================================================
    # SINGLE PARAMETER CASES
    # =====================================================

    else

        foreach val ($vals)

            set newDir = "${type}/${param_name}_${val}"
            echo "Creating run: $newDir"

            if (-d $newDir) then
                echo "Already exists — skipping"
                continue
            endif

            mkdir -p $type
            rsync -avz --exclude 'nrun2_*' --exclude '.git*' $BASEDIR/ $newDir

            cd $newDir

            sed -i "/${param_name}/s/[0-9.eEd+-]\+/${val}/" para1_in.dat

            if ($type == "Run_Apolar_cell_motility") then
                sed -i "/Apolar_cell_motility/s/.false./.true./" para1_in.dat

            else if ($type == "Run_Polar_cell_motility") then
                sed -i "/if_polar_motility/s/.false./.true./" para1_in.dat

            else if ($type == "Run_Mechano_chemical_regulation") then
                sed -i "/if_RhoROCK/s/.false./.true./" para1_in.dat
                sed -i "/if_RK4/s/.false./.true./" para1_in.dat
                sed -i "/if_coupling_noise/s/.false./.true./" para1_in.dat
            endif


            # -----------------------------
            # RUN SIMULATION
            # -----------------------------

            sh clear.sh
            sh compile.sh

            echo "Launching simulation..."
            nohup ./vertexmain.exe > nohup.out &

            cd ../../
            echo "Done + Launched: $newDir"

        end
    endif

    @ itype++
end

echo "======================================"
echo "All simulations created AND launched."
echo "======================================"

