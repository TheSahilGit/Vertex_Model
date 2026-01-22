#!/usr/bin/tcsh


# Another example to change para1_in.dat and submit batchjob. 
# this is from ___.___.__.40 machine, so change the path



# Set the main directory (adjust if necessary)
set MY_DIR=/mnt/hdd2/Sahil/VertexModel/VertexModel_NewCode_V4_applyPerturb_64x64

# Define the parameter values
set etas_arr = ("0.0d0" "0.01d0" "0.02d0" "0.04d0")
set sin_arr = (3 5 9 17 33)

# Loop over all combinations
foreach etas_max ($etas_arr)
    foreach sin_pert ($sin_arr)
        # Construct a new directory name which includes the parameter values
        set newDir = "Run_etas${etas_max}_sinPerturb${sin_pert}"
        echo "Creating new case in directory: $newDir"
        
        # Copy the BaseCase directory to the new directory
        cp -r BaseCase $newDir
        
        cd $newDir
        sed -i "/etas_max/s/0.00d0/${etas_max}/g" para1_in.dat
        sed -i "/sin_perturb_waveNumber/s/2/${sin_pert}/g" para1_in.dat

        cat para1_in.dat

        sh clear.sh
        sh compile.sh
        nohup ./vertexmain.exe  > nohup.out &

        cd ../

        
        echo "Modified para1_in.dat in $newDir"
    end
end

echo "All cases created."
