#!/usr/bin/env zsh
#
#### Job name 
#BSUB -J "prescribedInflow[1-4]"

#BSUB -o /home/ns179556/DG_moments_dealii/Tests/build/prescribedInflow.out


#BSUB -N -u neerajsarna587@gmail.com
#BSUB -W 0:30

#BSUB -M 1024

echo LSB_JOBINDEX: $LSB_JOBINDEX
cd /home/ns179556/DG_moments_dealii/Tests/build/

case "$LSB_JOBINDEX" in
        1)
                ./test.out 4 -p ../../Input_Files/2D_Systems/input_G20.in -p ../../Input_Files/Input_Free_Flow_Channel/input_free_flow_channel_prescribed_Vn_0p001.in
        ;;
        2)
                ./test.out 4 -p ../../Input_Files/2D_Systems/input_G20.in -p ../../Input_Files/Input_Free_Flow_Channel/input_free_flow_channel_prescribed_Vn_0p01.in
        ;;

        3)
                ./test.out 4 -p ../../Input_Files/2D_Systems/input_G20.in -p ../../Input_Files/Input_Free_Flow_Channel/input_free_flow_channel_prescribed_Vn_0p1.in
        ;;

        4)
                ./test.out 4 -p ../../Input_Files/2D_Systems/input_G20.in -p ../../Input_Files/Input_Free_Flow_Channel/input_free_flow_channel_prescribed_Vn_1p0.in
        ;;


esac

