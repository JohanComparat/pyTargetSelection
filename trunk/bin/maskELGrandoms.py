import numpy as n
import glob
import pymangle
m1=pymangle.Mangle("../starMasks/tycho2mask-0Vmag10.pol")
m2=pymangle.Mangle("../starMasks/tycho2mask-10Vmag11.pol")

def makeClusteringSample(tracer):

	ll=glob.glob(tracer+"_*_radec.gz")

	ras,decs=[],[]
	for el in ll:
		dat=n.loadtxt(el,unpack=True)
		print el, len(dat[0])
		bo_mask=m1.contains(dat[0],dat[1])
		bs_mask=m2.contains(dat[0],dat[1])
		outsideMask=(bo_mask==0)&(bs_mask==0)
		ra=dat[0][(dat[0]>=0)&(outsideMask)]
		dec=dat[1][(dat[0]>=0)&(outsideMask)]
		ras=n.hstack((ras,ra))
		decs=n.hstack((decs,dec))


	nR=len(ras)
	rd=n.random.rand(nR)
	threshold=0.1
	ok=(rd<threshold)
	n.savetxt(tracer+".random.masked", n.transpose([ ras[ok], decs[ok] , n.ones_like(ras[ok])]), fmt='%s')


tracer ="../decalsSelections/randoms/randoms_ELG_280"
makeClusteringSample(tracer)
tracer ="../decalsSelections/randoms/randoms_ELG_240"
makeClusteringSample(tracer)

