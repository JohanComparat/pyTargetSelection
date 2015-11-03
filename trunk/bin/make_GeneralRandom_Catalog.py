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
import glob

brickSumm="../decals/dr1/decals-bricks-in-dr1.fits"
cat=fits.open(brickSumm)
dat=cat[1].data

nR=3000
#sel=((dat['ra']<100))&
sel=((dat['has_g'])&(dat['has_r'])&(dat['has_z'])) # (dat['ra']>300)|
bricks,ramins,ramaxs,decmins,decmaxs=dat['brickname'],dat['ra1'],dat['ra2'],dat['dec1'],dat['dec2']

def makeRDCatalog(nR,brick,ramin,ramax,decmin,decmax):
	# make random points in delimit ra dec min max
	raR=n.random.uniform(low=ramin, high=ramax, size=nR)
	decR=n.random.uniform(low=decmin, high=decmax, size=nR)
	dpl_g="../decalsDepthMasks/decals-"+brick+"-depth-g.fits.gz"
	dpl_r="../decalsDepthMasks/decals-"+brick+"-depth-r.fits.gz"
	dpl_z="../decalsDepthMasks/decals-"+brick+"-depth-z.fits.gz"
	w = WCS(dpl_g)
	imageG=fits.open(dpl_g)
	x, y = w.all_world2pix(raR,decR, 0)
	value_g=n.array([imageG[0].data[int(x[ii])][int(y[ii])] for ii in range(len(x))])
	imageG.close()
	w = WCS(dpl_r)
	imageG=fits.open(dpl_r)
	x, y = w.all_world2pix(raR,decR, 0)
	value_r=n.array([imageG[0].data[int(x[ii])][int(y[ii])] for ii in range(len(x))])
	imageG.close()
	w = WCS(dpl_z)
	imageG=fits.open(dpl_z)
	x, y = w.all_world2pix(raR,decR, 0)
	value_z=n.array([imageG[0].data[int(x[ii])][int(y[ii])] for ii in range(len(x))])
	imageG.close()
	n.savetxt("../decalsRandoms/random-"+brick+".dat.gz", n.transpose([ raR , decR,value_g, value_r ,value_z]) , fmt= '%3.7f %3.7f %10.4f %10.4f %10.4f')

indexes=range(len(bricks[sel]))
n.random.shuffle(indexes)

for ii in indexes:
	print bricks[sel][ii]
	el=bricks[sel][ii]
	exist_g= glob.glob("../decalsDepthMasks/decals-"+ el+ "-depth-g.fits.gz")
	exist_r= glob.glob("../decalsDepthMasks/decals-"+ el+ "-depth-r.fits.gz")
	exist_z= glob.glob("../decalsDepthMasks/decals-"+ el+ "-depth-z.fits.gz")
	existRD=glob.glob("../decalsRandoms/random-"+el+".dat.gz")
	if len(existRD)==0 and len(exist_g)==1 and len(exist_r)==1 and len(exist_z)==1 :
		makeRDCatalog(nR,bricks[sel][ii],ramins[sel][ii], ramaxs[sel][ii],decmins[sel][ii], decmaxs[sel][ii])


