# Applies the bright star mask to Target catalog.

import numpy as n
import glob
import pymangle
import astropy.io.fits as fits
 
inputCatalog="../decalsSelections/catalog_ELG_240.fits"
outputCatalog="../decalsSelections/catalog_ELG_240-masked.fits"

# define the list of masks to be applied :

#mL=n.hstack(("../decalsMasks/bright_star_mask_pix.ply","../decalsMasks/bright_object_mask_rykoff_pix.ply",glob.glob("../decalsMasks/tycho2mask*.pol")))
mL=glob.glob("../starMasks/tycho*.pol")
 
# open the catalog
hd=fits.open(inputCatalog)
ras=hd[1].data['ra']
decs=hd[1].data['dec']
mask=n.zeros_like(ras)

for i in range(len(mL)):
    print "N=",i+1," stands for mask=",mL[i]
    m0=pymangle.Mangle(mL[i]) 
    bt_mask=m0.contains(ras,decs)
    mask[(bt_mask==1)]=n.ones_like(mask[(bt_mask==1)])*(i+1)


# adds a columns to the fits file "mask" where mask>=1 means it is inside one of the masks masked.
MaskCol = fits.Column(name="mask",format="I", array=mask )
new_columns = hd[1].data.columns + MaskCol

hdu = fits.BinTableHDU.from_columns(new_columns)
#hdu.writeto("catalogs/allELG240-masked.fits")
hdu.writeto(outputCatalog)
