import glob
import os
import numpy as n
import sys
import time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import Distance
from astropy import units as u
from astropy import wcs
from astropy.wcs import WCS

def getCompleteness(brick,raR,decR):
	""" get the completeness values for a list of ra and dec in a brick """
	# make random points in delimit ra dec min max
	dpl_g="../decalsDepthMasks/decals-"+brick+"-depth-g.fits.gz"
	dpl_r="../decalsDepthMasks/decals-"+brick+"-depth-r.fits.gz"
	dpl_z="../decalsDepthMasks/decals-"+brick+"-depth-z.fits.gz"
	# print dpl_g
	w = WCS(dpl_g)
	imageG=fits.open(dpl_g)
	x, y = w.all_world2pix(raR,decR, 0)
	value_g=n.array([imageG[0].data[int(x[ii])][int(y[ii])] for ii in range(len(x))])
	imageG.close()
	# print "1"
	w = WCS(dpl_r)
	imageG=fits.open(dpl_r)
	x, y = w.all_world2pix(raR,decR, 0)
	value_r=n.array([imageG[0].data[int(x[ii])][int(y[ii])] for ii in range(len(x))])
	imageG.close()
	# print "2"
	w = WCS(dpl_z)
	imageG=fits.open(dpl_z)
	x, y = w.all_world2pix(raR,decR, 0)
	value_z=n.array([imageG[0].data[int(x[ii])][int(y[ii])] for ii in range(len(x))])
	imageG.close()
	# print "comp ok"
	return value_g, value_r ,value_z

def selectELG_280(brickPath, brick, gLimitC, rLimitC, zLimitC):
	""" selects the ELGs from a brick to a completeness limit in grz with the new 280 algorithm (still under tests). """
	hdu=fits.open(brickPath)
	dat=hdu[1].data
	hdu.close()
	# print "brick opened"
	# column re-definition
	g = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[1] / dat['decam_mw_transmission'].T[1])
	r = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[2] / dat['decam_mw_transmission'].T[2])
	z = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[4] / dat['decam_mw_transmission'].T[4])
	radius=n.max(n.transpose([dat['shapeExp_r'],dat['shapeDev_r']]),axis=1)
	gAp3=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[3][1] / dat['decam_mw_transmission'].T[1] )
	#W1 = 22.5 - 2.5 * n.log10(dat['wise_flux'].T[0] / dat['wise_mw_transmission'].T[0])
	# bad detection exclusion
	noJunk=(gAp3>21.5)&(radius<2)&(dat['brick_primary'])&(dat['decam_anymask'].T[1]==0)& (dat['decam_anymask'].T[2]==0)&(dat['decam_anymask'].T[4]==0)
	# color selection
	color240=(21.5< g)&(g < 22.8)&(gAp3<22.8+0.5)&( 0.2< g-r)&(g-r< 0.7)&( 0.25 < r-z)&(r-z < 1.4)&(  r-z > 0.45*(g-r) + 0.4)&( r-z < 0.8*(g-r) + 1)
	color280=(21.5< g)&(g < 22.8)&(gAp3<22.8+0.5)&( 0.3< g-r)&(g-r< 0.1+0.67*(r-z))&(g-r< 2.0-(r-z))&(0.6<r-z )#&( 0.5< r-W1)
	color=(color280)|(color240)
	selection=(noJunk)&(color)
	# completeness cuts
	value_g=n.zeros_like(dat['dec'])
	value_r=n.zeros_like(dat['dec'])
	value_z=n.zeros_like(dat['dec'])
	# print len(value_g),len(selection),"selected",len(selection.nonzero()[0])
	if len(selection.nonzero()[0]) > 0 :
		value_g[selection], value_r[selection] ,value_z[selection]=getCompleteness( brick, dat['ra'][selection], dat['dec'][selection])
		# selection
		observedCompleteGRZ=(value_g>gLimitC)&(value_r>rLimitC)&(value_z>zLimitC)
		sel=(selection)&(observedCompleteGRZ)
		# print value_g[selection], value_r[selection], value_z[selection]
		print "n selected",len(dat['brickname'][sel])
		if len(sel.nonzero()[0])>0:
			new_columns = dat.columns
			hdu2 = fits.BinTableHDU.from_columns(new_columns)
			hdu2.data = hdu2.data[sel]
			return hdu2
		else :
			return -1
	else :
		return -1

def selectELG_240(brickPath, brick, gLimitC, rLimitC, zLimitC):
	""" selects the ELGs from a brick to a completeness limit in grz. """
	hdu=fits.open(brickPath)
	dat=hdu[1].data
	hdu.close()
	# print "brick opened"
	# column re-definition
	g = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[1] / dat['decam_mw_transmission'].T[1])
	r = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[2] / dat['decam_mw_transmission'].T[2])
	z = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[4] / dat['decam_mw_transmission'].T[4])
	radius=n.max(n.transpose([dat['shapeExp_r'],dat['shapeDev_r']]),axis=1)
	gAp3=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[3][1] / dat['decam_mw_transmission'].T[1] )
	# bad detection exclusion
	noJunk=(gAp3>21.5)&(radius<2)&(dat['brick_primary'])&(dat['decam_anymask'].T[1]==0)& (dat['decam_anymask'].T[2]==0)&(dat['decam_anymask'].T[4]==0)
	# color selection
	color=(21.5< g)&(g < 22.8)&(gAp3<22.8+0.5)&( 0.2< g-r)&(g-r< 0.7)&( 0.25 < r-z)&(r-z < 1.4)&(  r-z > 0.45*(g-r) + 0.4)&( r-z < 0.8*(g-r) + 1)
	selection=(noJunk)&(color)
	# completeness cuts
	value_g=n.zeros_like(dat['dec'])
	value_r=n.zeros_like(dat['dec'])
	value_z=n.zeros_like(dat['dec'])
	# print len(value_g),len(selection),"selected",len(selection.nonzero()[0])
	if len(selection.nonzero()[0]) > 0 :
		value_g[selection], value_r[selection] ,value_z[selection]=getCompleteness( brick, dat['ra'][selection], dat['dec'][selection])
		# selection
		observedCompleteGRZ=(value_g>gLimitC)&(value_r>rLimitC)&(value_z>zLimitC)
		sel=(selection)&(observedCompleteGRZ)
		# print value_g[selection], value_r[selection], value_z[selection]
		print "n selected",len(dat['brickname'][sel])
		if len(sel.nonzero()[0])>0:
			new_columns = dat.columns
			hdu2 = fits.BinTableHDU.from_columns(new_columns)
			hdu2.data = hdu2.data[sel]
			return hdu2
		else :
			return -1
	else :
		return -1

def selectLRG_DESI(brickPath, brick, gLimitC, rLimitC, zLimitC):
	""" selects the LRG for DESI from a brick to a completeness limit in grzW1. """
	hdu=fits.open(brickPath)
	dat=hdu[1].data
	hdu.close()
	# print "brick opened"
	# column re-definition
	g = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[1] / dat['decam_mw_transmission'].T[1])
	r = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[2] / dat['decam_mw_transmission'].T[2])
	z = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[4] / dat['decam_mw_transmission'].T[4])
	w1 = 22.5 - 2.5 * n.log10(dat['wise_flux'].T[0] / dat['wise_mw_transmission'].T[0])
	radius=n.max(n.transpose([dat['shapeExp_r'],dat['shapeDev_r']]),axis=1)
	gAp3=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[3][1] / dat['decam_mw_transmission'].T[1] )
	# bad detection exclusion
	noJunk=(gAp3>21.5)&(radius<2)&(dat['brick_primary'])&(dat['decam_anymask'].T[1]==0)& (dat['decam_anymask'].T[2]==0)&(dat['decam_anymask'].T[4]==0)
	# color selection
	color=(z<19.2)&(0.9<r-z)&(r-z<1.7)&(0.5<z-w1)&(z-w1<1.3)&(z-w1>r-z-0.55) &(z-w1<r-z-0.1)
	selection=(noJunk)&(color)
	# completeness cuts
	value_g=n.zeros_like(dat['dec'])
	value_r=n.zeros_like(dat['dec'])
	value_z=n.zeros_like(dat['dec'])
	# print len(value_g),len(selection),"selected",len(selection.nonzero()[0])
	if len(selection.nonzero()[0]) > 0 :
		value_g[selection], value_r[selection] ,value_z[selection]=getCompleteness( brick, dat['ra'][selection], dat['dec'][selection])
		# selection
		observedCompleteGRZ=(value_g>gLimitC)&(value_r>rLimitC)&(value_z>zLimitC)
		sel=(selection)&(observedCompleteGRZ)
		# print value_g[selection], value_r[selection], value_z[selection]
		print "n selected",len(dat['brickname'][sel])
		if len(sel.nonzero()[0])>0:
			new_columns = dat.columns
			hdu2 = fits.BinTableHDU.from_columns(new_columns)
			hdu2.data = hdu2.data[sel]
			return hdu2
		else :
			return -1
	else :
		return -1


def selectRandoms(brick,gLimitC, rLimitC, zLimitC):
	""" selects the randoms from a brick to a completeness limit in grz. """
	raRA , decRA,value_gA, value_rA ,value_zA= n.loadtxt( "../decalsRandoms/random-" +brick +".dat.gz", unpack=True)
	observedCompleteGRZ=(value_gA>gLimitC)&(value_rA>rLimitC)&(value_zA>zLimitC)
	print "randoms selected",len(raRA), len(decRA[observedCompleteGRZ])
	return raRA[observedCompleteGRZ], decRA[observedCompleteGRZ]
