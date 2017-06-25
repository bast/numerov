#!/bin/bash

python ../numerov.py input.yml > output.yml
python pnc.py output.yml > pnc-output.yml
python plot.py output.yml plot.png

exit 0
