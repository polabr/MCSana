#!bash/bin

source /grid/fermiapp/products/uboone/setup_uboone.sh

setup root v6_04_02 -q e7:prof


cd /uboone/app/users/pabraten/larlite

source config/setup.sh
#source config/setup_reco.sh

source /uboone/app/users/cadams/pystack/setup.sh

#export OPENCV_INCDIR=/uboone/app/users/vgenty/opencv/include
#export OPENCV_LIBDIR=/uboone/app/users/vgenty/opencv/lib