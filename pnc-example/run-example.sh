#!/bin/bash

python ../numerov.py --general=general.yml --force-field=force.yml --exp-values=property.yml > numerov-output.yml
python pnc.py numerov-output.yml > pnc-output.yml
python plot.py --force-field=force.yml --exp-values=property.yml --numerov-output=numerov-output.yml --img=plot.png

exit 0
