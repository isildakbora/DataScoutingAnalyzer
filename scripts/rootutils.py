from ROOT import TFile, TH1F, TMath, TF1, TCanvas, TLegend, gROOT, gPad, TObject, TPad, TGraphAsymmErrors, TGraph, TGraphErrors
from array import array
from numpy import divide, sqrt
import math, sys

def setTDRStyle():
	gROOT.ProcessLine('.X ~/setTDRStyle.C')

def MakeGraph(h, LUMI):
	vx, vexl, vexh, vy, veyl, veyh = array('d'), array('d'), array('d'), array('d'), array('d'), array('d')
	a = 0.3173/2.
	for i in range(h.GetNbinsX()):
		dx = h.GetBinWidth(i+1)
		x  = h.GetBinCenter(i+1)
		y  = h.GetBinContent(i+1)
		norm = 1./(dx*LUMI)
		print x, norm*y
		if y > 0.:
			vx.append(x)
			vy.append(norm*y)
			vexl.append(dx/2.)
			vexh.append(dx/2.)
			if y > 30.:
				el = norm*sqrt(y)
				eh = norm*sqrt(y)
			else:
				yl = 0.5*TMath.ChisquareQuantile(a,2.*y)
				yh = 0.5*TMath.ChisquareQuantile(1-a,2.*(y+1))
				el = norm*(y-yl)
				eh = norm*(yh-y)
			veyl.append(el)
			veyh.append(eh)
	g = TGraphAsymmErrors(len(vx), vx, vy, vexl, vexh, veyl, veyh)
	return g

#TGraph to python array
def TGraph2Array(graph, xy):
	data = array('d')
	for i in range(0, graph.GetN()):
		if xy == 'x':
			data.append(graph.GetX()[i])
		elif xy == 'y':
			data.append(graph.GetY()[i])
		else:
			print 'Please indicate x or y!'
	return data

def progressbar(progress):
	sys.stdout.write('\r['+int(progress)*'|'+'%'+str(2*progress)+(50-int(progress))*' '+']')
	sys.stdout.flush()

def keepGUIalive():
	rep = ''
	while not rep in ['q','Q']:
		rep = raw_input('enter "q" to quit: ')
		if 1 < len(rep):
			rep = rep[0]


def CrystalBall(x, par):
	xcur = x[0]
	alpha = par[0]
	n = par[1]
	mu = par[2]
	sigma = par[3]
	N = par[4]
	exp = TF1("exp","exp(x)",1e-20,1e20)

	if alpha < 0.:
		A = pow((n/(-1*alpha)),n)*exp.Eval((-1)*alpha*alpha/2)
		B = n/(-1*alpha) + alpha
	else:
		A = pow((n/alpha),n)*exp.Eval((-1)*alpha*alpha/2)
		B = n/alpha - alpha

	if (xcur-mu)/sigma > (-1)*alpha:
		f = N*exp.Eval((-1)*(xcur-mu)*(xcur-mu)/ (2*sigma*sigma))
	else:
		f = N*A*pow((B-(xcur-mu)/sigma),(-1*n))
	return f
