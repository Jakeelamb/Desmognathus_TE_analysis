#!/usr/bin/env python3

import os

for file in os.listdir("/nfs/home/jlamb/Projects/RepeatMasker/Data/"):
    if file.endswith(".fna"):
        os.system("sbatch /nfs/home/jlamb/Projects/RepeatMasker/Scripts/repeat_mask.sh " + file)
