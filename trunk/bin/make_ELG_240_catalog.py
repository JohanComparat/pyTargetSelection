from lib_select_ELG import *

suffix="ELG_240" # suffix to be written in all output files

brickSumm="../decals/dr1/decals-bricks-in-dr1.fits"
cat=fits.open(brickSumm)
dat=cat[1].data

ids=n.arange(len(dat['ra']))
RRAS=n.arange(0.,360.1,2.)
for jj in range(len(RRAS)-1):
	suffixRD="_"+str(RRAS[jj])+"_ra_"+str(RRAS[jj+1])+"_"
	print suffixRD
	print "--------------------------------------------------"
	sel=(dat['ra']>=RRAS[jj])& (dat['ra']<RRAS[jj+1])&((dat['has_g'])& (dat['has_r'])&(dat['has_z'])) # (dat['ra']>300)|
	raRs,decRs=[],[]
	if len(sel.nonzero()[0])>0:
		for ii in ids[sel]:
			el=dat['brickname'][ii]
			print "========"
			print "brick no",el
			path_to_brick = "../decals/dr1/tractor/"+el[:3]+"/tractor-"+el+".fits"
			exist = glob.glob(path_to_brick)
			exist_g = glob.glob( "../decalsDepthMasks/decals-" + el + "-depth-g.fits.gz" )
			exist_r = glob.glob( "../decalsDepthMasks/decals-" + el + "-depth-r.fits.gz" )
			exist_z = glob.glob( "../decalsDepthMasks/decals-" + el + "-depth-z.fits.gz" )
			exist_RD = glob.glob( "../decalsRandoms/random-" + el + ".dat.gz" )
			brickOutName = "../decalsSelections/bricks/" + path_to_brick.split('/')[-1][:-5]+ "-"+suffix+".fits"
			exist_selection = glob.glob( brickOutName )
			# print exist
			if len(exist)==1 and len(exist_g)==1 and len(exist_r)==1 and len(exist_z)==1 and len(exist_RD)==1 and len(exist_selection)==0:
				zLimitVAL=12.755
				rLimitVAL=30.05661087
				gLimitVAL=62.79716079 # 5 sigma z,r,g lim 22.6, 22.6, 23
				hduBis = selectELG_240(path_to_brick,el,gLimitVAL,rLimitVAL,zLimitVAL)
				if hduBis == -1 :
					continue
				else:
					hduOut = hduBis
					hduOut.writeto( brickOutName )
					raRI, decRI = selectRandoms(el, gLimitVAL, rLimitVAL, zLimitVAL)
					raRs = n.hstack(( raRs, raRI ))
					decRs = n.hstack(( decRs, decRI ))
			else:
				print "pass" , len(exist), len(exist_g), len(exist_r), len(exist_z), len(exist_RD), len(exist_selection)


		DAT =  n.column_stack((raRs,decRs))
		n.savetxt("../decalsSelections/randoms_"+suffix+suffixRD+"radec.gz",DAT,fmt= '%s')


# java commands to concatenate bricks into a complete catalog
print "now concatenates the catalog"
c1="ls ../decalsSelections/bricks/*"+suffix+".fits > flAll"
c2="""java -jar stilts.jar tcat ifmt=fits in=@flAll out=../decalsSelections/catalog_""" +suffix+ """.fits"""
os.system(c1)
os.system(c2)
c3="rm flAll"
os.system(c3)

