#!/bin/bash

python ../numerov.py input.yml > output.yml
python pnc.py output.yml > pnc-output.yml

exit 0
