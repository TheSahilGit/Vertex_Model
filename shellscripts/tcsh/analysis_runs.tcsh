#!/usr/bin/tcsh

# ==========================================================
# Choose analysis (JUST CHANGE THIS LINE)
# ==========================================================

set ANALYSIS_NAME = "Analysis_Circularity"
#set ANALYSIS_NAME = "Analysis_MSD_cellID"
#set ANALYSIS_NAME = "Analysis_Qt"

set ANALYSIS_SRC = "Vertex_Model/analysis/${ANALYSIS_NAME}.m"

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

        # Ensure analysis/analysisdata directories exist (log.txt:
        # Analysis_Circularity.m/Analysis_MSD_cellID.m now write their
        # output to ../analysisdata/ instead of alongside the scripts)
        mkdir -p $rundir/analysis
        mkdir -p $rundir/analysisdata

        # Copy MATLAB analysis code
        cp $ANALYSIS_SRC $rundir/analysis/

        # Run MATLAB headless
        cd $rundir/analysis

        #pwd
        if (-e ../analysisdata/circularity.dat) then
          rm ../analysisdata/circularity.dat
        endif
        if (-e ../analysisdata/msd.dat) then
          rm ../analysisdata/msd.dat
        endif
        if (-e ../analysisdata/Qt.dat) then
          rm ../analysisdata/Qt.dat
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

        rm ../analysisdata/${ANALYSIS_NAME}_${runname}.dat

        if (-e ../analysisdata/circularity.dat) then
            mv ../analysisdata/circularity.dat ../analysisdata/${ANALYSIS_NAME}_${runname}.dat
        endif

        if (-e ../analysisdata/msd.dat) then
            mv ../analysisdata/msd.dat ../analysisdata/${ANALYSIS_NAME}_${runname}.dat
        endif

        if (-e ../analysisdata/Qt.dat) then
            mv ../analysisdata/Qt.dat ../analysisdata/${ANALYSIS_NAME}_${runname}.dat
        endif

        cd ../../../
        echo "Done: $rundir"

    end
end

echo "======================================"
echo "ALL ANALYSES COMPLETED"
echo "======================================"

