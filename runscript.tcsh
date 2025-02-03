#!/usr/bin/tcsh

# Set the working directory
set MY_DIR=/media/data/sahil/Work/VertexModel_NewCode_V4_Yield/

cd BaseCode

# Get the value from para1_in.dat
set sudden_shearStrength0 = `cat para1_in.dat | grep sudden_shearStrength | awk '{print $1}'`
set sudden_shearStrength_in = `echo $sudden_shearStrength0 | sed 's/d/E/'`
set sudden_shearStrength_dec = `printf "%.15f\n" "$sudden_shearStrength_in"`

echo "Decimal format: $sudden_shearStrength_dec"

cd ../ 

# Define the array with Fortran-style values
set etas_max = (0 1.0d-3 2.0d-3 5d-3 1d-2 2d-2 4d-2)

# Loop through the array
foreach value ($etas_max)
    # Convert 'd' to 'E' for Unix compatibility
    set etas_max_unix = `echo $value | sed 's/d/E/'`
    
    # Convert to decimal format for calculations if needed
    set etas_max_dec = `echo $etas_max_unix | awk '{printf "%.15f\n", $1}'`
    
    echo "Converted etas_max_dec: $etas_max_dec"
    
    # Use sudden_shearStrength_dec in the loop
    set sudden_shearStrength_use = $sudden_shearStrength_dec  # This is the correct way to assign
    
    # Loop to perform repeated actions
    set jj = 1 
    while ($jj <= 7)
        # Use the converted values in the cp command
        cp -r BaseCode/ Mot_${etas_max_dec}_Shear_${sudden_shearStrength_use}
    
        # Go into the directory
        cd Mot_${etas_max_dec}_Shear_${sudden_shearStrength_use}
    
        # Update the values in para1_in.dat
        sed -i "/etas_max_/s/0.0d0/$etas_max_dec/g" para1_in.dat
        sed -i "/sudden_shearStrength/s/1.0d-7/$sudden_shearStrength_use/g" para1_in.dat
    
        echo "Updated para1_in.dat with new values"
        cat para1_in.dat

        sh clear.sh
        sh compile.sh
        nohup ./vertexmain.exe > nohup.out &
    
        # Multiply by 10 for next iteration
        set sudden_shearStrength_use = `echo "$sudden_shearStrength_use * 10.0" | bc -l` 
    
        # Increment jj
        @ jj++
    
        # Return to the previous directory
        cd ../ 
    end 
end

