#!/usr/bin/env python

import os

thelist=list(range(1200,1215))
for j in thelist:
	print("\nhistograms_Run%05d.root\n" % (j))
	os.system("python3 reconstruction.py configFile.txt -r %d -j4" % (j))


