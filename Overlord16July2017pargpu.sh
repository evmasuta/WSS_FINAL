# Kills all zombie instances of MATLAB.  Remove if you dare.  But seriously, don't, this stupid error sucks.
killall MATLAB
# options: (executable inputPath outputPath exptName timesteps parallelitySwitch parallelityCutoff onlyVel}
./shear16July2017pargpu.sh 2016Oct24_Ex1951_4dflow_Ser12 2016Oct24_Ex1951_4dflow_Ser12_out EMM037A_final 20 1 0 1 no &&
./shear16July2017pargpu.sh 2016Oct27_Ex1985_4dflow_Ser10 2016Oct27_Ex1985_4dflow_Ser10_out EMM037B_final 20 1 0 1 no &&
./shear16July2017pargpu.sh 2015Sep23_Ex648_4dflow_Ser12 2015Sep23_Ex648_4dflow_Ser12_out EMM037D_final 20 1 0 1 no &&
./shear16July2017pargpu.sh 2016Nov15_Ex14970_4dflow_Ser8 2016Nov15_Ex14970_4dflow_Ser8_out EMM037E_final 20 1 0 1 no &&
killall MATLAB