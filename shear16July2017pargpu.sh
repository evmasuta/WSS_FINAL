# N.B. The following syntax must be used for this shell script:
# ./gaussGradOverlord_xxxxxx {input PATH w/o EXTENSION} {output PATH} {Experiment Title} {dimT} {Gauss Dim} {run package.flow?} 

# clean any crashed files
rm -r "$1"
rm -r "$2"
rm -r "$2vel"
# phase 0: unzip raw data files
unzip "$1.zip"

cd $1/mag
numFiles=$(ls -1 | wc -l)
cd ../..

rm -r "$2"
mkdir "$2"
mkdir "$2vel"
mkdir "$2/mag"
mkdir "$2/vel_0_corr"
mkdir "$2/vel_1_corr"
mkdir "$2/vel_2_corr"
mkdir "$2vel/mag"
mkdir "$2vel/vel_0_corr"
mkdir "$2vel/vel_1_corr"
mkdir "$2vel/vel_2_corr"
# https://www.mathworks.com/matlabcentral/newsreader/view_thread/90363
# https://www.mathworks.com/matlabcentral/answers/39067-run-matlab-function-with-arguments-on-linux-terminal
TEST="trueShearv17_RECONparGPU('$1','$2','$3','$4','$numFiles','$5','$6','$7');quit;"
echo $TEST
/usr/local/MATLAB/R2016b/bin/matlab -r $TEST

cp -r "$1/mag" "$2vel"
cp -r "$1/vel_0_corr" "$2vel"
cp -r "$1/vel_1_corr" "$2vel"
cp -r "$1/vel_2_corr" "$2vel"

zip -r $2.zip $2
zip -r $2vel.zip $2vel

rm -r "$1"
rm -r "$2"
rm -r "$2vel"

/home/flow-user/scripts/curl-hsiaolab $2.zip


rm -r $2.zip
rm -r $2vel.zip
