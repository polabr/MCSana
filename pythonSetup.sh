# From David C's guide: https://github.com/davidc1/Documentation/wiki/Installing-python-packages-on-the-gpvms
source /grid/fermiapp/products/uboone/setup_uboone.sh
setup uboonecode v06_26_01_03 -qe10:prof
export PYTHONPATH=/uboone/app/users/$USER/python_install/lib/python2.7/site-packages:$PYTHONPATH