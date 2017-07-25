#!/bin/bash

virtualenv --system-site-packages venv
source venv/bin/activate
pip install /home/bast/pnc/numerov

source ~/.localrc

cooley --general=general.yml \
       --reduced-mass=force.yml \
       --force-field=force.yml \
       --exp-values=property.yml > numerov-output.yml

python pnc.py numerov-output.yml > pnc-output.yml

python plot.py --force-field=force.yml \
               --exp-values=property.yml \
               --numerov-output=numerov-output.yml \
               --img=plot.png

exit 0
