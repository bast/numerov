#!/bin/bash

cooley --general=general.yml \
       --reduced-mass=force.yml \
       --force-field=force.yml \
       --exp-values=property.yml > numerov-output.yml

python plot.py --numerov-output=numerov-output.yml \
               --img=plot.png

exit 0
