#!/usr/bin/tcsh


# Another example to change para1_in.dat and submit batchjob. 
# this is from Dobby, so change the path



# Set the main directory (adjust if necessary)
set MY_DIR=/media/data/sahil/VertexModel/

# Define the parameter values
set etas_arr = ("0.0d0" "0.005d0" "0.01d0" "0.02d0" "0.04d0")
set oscl_freq_arr = ("1.0d-4" "1.0d-3" "1.0d-2" "1.0d-1" "1.0d0" "1.0d1" "1.0d2")

# Loop over all combinations
foreach etas_max ($etas_arr)
    foreach oscl_freq ($oscl_freq_arr)
        # Construct a new directory name which includes the parameter values
        set newDir = "Run_etas${etas_max}_oscl_freq${oscl_freq}"
        echo "Creating new case in directory: $newDir"
        
        # Copy the BaseCase directory to the new directory
        cp -r Vertex_Model $newDir
        
        cd $newDir
        sed -i "/etas_max/s/0.00d0/${etas_max}/g" para1_in.dat
        sed -i "/Oscl_freq_wo/s/1.0d-4/${oscl_freq}/g" para1_in.dat

        cat para1_in.dat

        sh clear.sh
        sh compile.sh
        nohup ./vertexmain.exe  > nohup.out &

        cd ../

        
        echo "Modified para1_in.dat in $newDir"
    end
end

echo "All cases created."
