#!/usr/bin/tcsh

# ==========================================================
# Choose analysis (JUST CHANGE THIS LINE)
# ==========================================================

set ANALYSIS_NAME = "Analysis_Circularity"
#set ANALYSIS_NAME = "Analysis_MSD_cellID"
#set ANALYSIS_NAME = "Analysis_Qt"

set ANALYSIS_SRC = "Vertex_Model/Analysis_Codes_matlab/${ANALYSIS_NAME}.m"

# ==========================================================
# Run types (top-level directories)
# ==========================================================

set types = ( \
    "Run_Apolar_cell_motility" \
    "Run_Polar_cell_motility" \
    "Run_Fluctuating_contractility" \
    "Run_Mechano_chemical_regulation" \
)

# ===  For only Polar Case ====
#set types = ( \
#    "Run_Polar_cell_motility" \
#)

echo "======================================"
echo "Starting MATLAB analysis (no display)"
echo "Analysis selected: $ANALYSIS_NAME"
echo "======================================"

# ==========================================================
# Loop over all run directories
# ==========================================================
foreach type ($types)

    if (! -d $type) then
        echo "Skipping missing directory: $type"
        continue
    endif

    foreach rundir (`ls -d ${type}/*`)

        if (! -d $rundir) then
            continue
        endif

        echo "Processing analysis in: $rundir"

        # Ensure analysis directory exists
        mkdir -p $rundir/Analysis_Codes_matlab

        # Copy MATLAB analysis code
        cp $ANALYSIS_SRC $rundir/Analysis_Codes_matlab/

        # Run MATLAB headless
        cd $rundir/Analysis_Codes_matlab

        #pwd
        if (-e circularity.dat) then
          rm circularity.dat
        endif
        if (-e msd.dat) then
          rm msd.dat
        endif
        if (-e Qt.dat) then
          rm Qt.dat
        endif

        matlab -nodisplay -nosplash -nodesktop << EOF > matlab.out
try
    $ANALYSIS_NAME
    exit
catch ME
    disp(getReport(ME))
    exit(1)
end
EOF

        # Rename output if present (generic + safe)
        set runname = `basename $rundir`

        # echo $runname

        rm ${ANALYSIS_NAME}_${runname}.dat

        if (-e circularity.dat) then
            mv circularity.dat ${ANALYSIS_NAME}_${runname}.dat
        endif

        if (-e msd.dat) then
            mv msd.dat ${ANALYSIS_NAME}_${runname}.dat
        endif

        if (-e Qt.dat) then
            mv Qt.dat ${ANALYSIS_NAME}_${runname}.dat
        endif

        cd ../../../
        echo "Done: $rundir"

    end
end

echo "======================================"
echo "ALL ANALYSES COMPLETED"
echo "======================================"

