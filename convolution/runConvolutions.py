import os
from itertools import product
signaltypes = ['RSGravitonToGG', 'RSGravitonToQQbar', 'QstarToJJ']
for mass, signaltype in product(xrange(300, 6600, 100), signaltypes):
	command = 'python convolution_v3.py '+str(mass)+' 1e6 '+ signaltype
	print command
	os.system(command)
