#!/usr/bin/tcsh

# -----------------------------
# Base directory
# -----------------------------
set BASEDIR = "Vertex_Model"

# -----------------------------
# Polar parameter values
# -----------------------------
## vo values

#set values2   = ("1.0d-1" "2.0d-1" "2.5d-1" "3.0d-1")

## Dr values
#set values2_2 = ("1.5d0" "2.0d0" "3.0d0" "4.0d0")



set values2   = ("2.5d-1")
set values2_2 = ("1.0d-2" "5.0d-2" "1.0d-1" "2.0d-1" "2.5d-1")



# -----------------------------
# Target type
# -----------------------------
set type = "Run_Polar_cell_motility"

echo "======================================"
echo "Processing Polar cell motility runs"
echo "======================================"

# -----------------------------
# Loop over vo and Dr
# -----------------------------
foreach vo ($values2)
    foreach Dr ($values2_2)

        set newDir = "${type}/vo_${vo}_Dr_${Dr}"
        echo "Creating run: $newDir"

        mkdir -p $type
        rsync -av --exclude 'nrun2_*' --exclude '.git*' $BASEDIR/ $newDir

        cd $newDir

        # -----------------------------
        # Modify parameters
        # -----------------------------
        sed -i "/vo/s/[0-9.eEd+-]\+/${vo}/" para1_in.dat
        sed -i "/Dr/s/[0-9.eEd+-]\+/${Dr}/" para1_in.dat

        # -----------------------------
        # Polar-specific flags (UNCHANGED)
        # -----------------------------
        sed -i "/if_ABP/s/.false./.true./" para1_in.dat

        # -----------------------------
        # Run simulation
        # -----------------------------

        sh clear.sh
        sh compile.sh
        nohup ./vertexmain.exe > nohup.out &

        cd ../../
        echo "Done: $newDir"

    end
end

echo "======================================"
echo "All Polar (vo, Dr) simulations created"
echo "======================================"
