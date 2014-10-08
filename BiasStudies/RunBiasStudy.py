#!/usr/bin/env python
import getopt, ROOT, math, sys, os, subprocess, string, re
from rootutils import *
from ROOT import *
from array import array 

signal_fraction       = array('d')
signal_fraction_error = array('d')

bias_1 = []

for mass in xrange(300, 500, 100):
	
	h_bias_1 = TH1F("bias_1","bias_1", 100, -1, 1)
	h_bias_1.Sumw2()
	bias_1.append(h_bias_1)

	for i in xrange(1000):
		cmd = "./BiasStudyRooFit.py " + str(mass)
		print "Running: " + cmd
		proc        = subprocess.Popen( cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
		output      = proc.communicate()[0]
		outputlines = output.split("\n")

		#print output
		for line in outputlines:
			if re.search("signal_fraction", line):
				print line.split()
				signal_fraction.append(float(line.split()[1]))
				signal_fraction_error.append(float(line.split()[2]))
				print signal_fraction[len(signal_fraction)-1]
				print signal_fraction_error[len(signal_fraction_error)-1]
