# this scsript converts a masked ELG catalog to a fits file that can go through the tiling code 
import numpy as n
import astropy.io.fits as fits

inputCat="../decalsSelections/catalog_ELG_280-masked.fits"
outputCat="../decalsSelections/catalog_ELG_280-masked-readyForTiling.fits"

 
# open the catalog
hd=fits.open(inputCat)
dat=hd[1].data

magAp_1a_g=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[2][1] / dat['decam_mw_transmission'].T[1] )
magAp_1a_r=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[2][2] / dat['decam_mw_transmission'].T[2] )
magAp_1a_z=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[2][4] / dat['decam_mw_transmission'].T[4] )
magAp_1a_u = n.ones_like(magAp_1a_g) * 40.
magAp_1a_i = n.ones_like(magAp_1a_g) * 40.

c0 = fits.Column(name="ID",format='D', array=n.empty_like(magAp_1a_g) )

#c0 = fits.Column(name="FIBER2MAG_u",format='D', array=magAp_1a_u )
#c1 = fits.Column(name="FIBER2MAG_g",format='D', array=magAp_1a_g )
#c2 = fits.Column(name="FIBER2MAG_r",format='D', array=magAp_1a_r )
#c3 = fits.Column(name="FIBER2MAG_i",format='D', array=magAp_1a_i )
#c4 = fits.Column(name="FIBER2MAG_z",format='D', array=magAp_1a_z )
FIBERMAGMatrix = n.array([ magAp_1a_u, magAp_1a_g, magAp_1a_r, magAp_1a_i, magAp_1a_z ]).T
c1 = fits.Column(name="FIBER2MAG",format='5D', array=FIBERMAGMatrix )

magAp_1a_g=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[4][1] / dat['decam_mw_transmission'].T[1] )
magAp_1a_r=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[4][2] / dat['decam_mw_transmission'].T[2] )
magAp_1a_z=22.5 - 2.5 * n.log10(dat['decam_apflux'].T[4][4] / dat['decam_mw_transmission'].T[4] )
magAp_1a_u = n.ones_like(magAp_1a_g) * 40.
magAp_1a_i = n.ones_like(magAp_1a_g) * 40.

FIBERMAGMatrix = n.array([ magAp_1a_u, magAp_1a_g, magAp_1a_r, magAp_1a_i, magAp_1a_z ]).T
c2 = fits.Column(name="FIBERMAG",format='5D', array=FIBERMAGMatrix )

#c5 = fits.Column(name="FIBERMAG_u",format='D', array=magAp_1a_u )
#c6 = fits.Column(name="FIBERMAG_g",format='D', array=magAp_1a_g )
#c7 = fits.Column(name="FIBERMAG_r",format='D', array=magAp_1a_r )
#c8 = fits.Column(name="FIBERMAG_i",format='D', array=magAp_1a_i )
#c9 = fits.Column(name="FIBERMAG_z",format='D', array=magAp_1a_z )

def create_5_emptyColsInt(name):
	FIBERMAGMatrix = n.array([ n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g) ]).T
	c_1 = fits.Column(name=name,format='5I', array=FIBERMAGMatrix )
	return c_1

def create_5_emptyCols(name):
	FIBERMAGMatrix = n.array([ n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g), n.empty_like(magAp_1a_g) ]).T
	c_1 = fits.Column(name=name,format='5D', array=FIBERMAGMatrix )
	return c_1

c3 = create_5_emptyColsInt("CALIB_STATUS")

c4 = fits.Column(name="EPOCH",format='D', array=n.empty_like(magAp_1a_g) )
c5 = fits.Column(name="PMRA",format='D', array=n.empty_like(magAp_1a_g) )
c6 = fits.Column(name="PMDEC",format='D', array=n.empty_like(magAp_1a_g) )
c7 = fits.Column(name="EBOSS_TARGET_ID",format='D', array=n.empty_like(magAp_1a_g) )
c8 = fits.Column(name="THING_ID_TARGETING",format='I', array=n.empty_like(magAp_1a_g) )

c9 = fits.Column(name="OBJC_TYPE",format='I', array=n.empty_like(magAp_1a_g) )
c10 = fits.Column(name="OBJC_FLAGS",format='I', array=n.empty_like(magAp_1a_g) )
c11 = fits.Column(name="OBJC_FLAGS2",format='I', array=n.empty_like(magAp_1a_g) )
c12 = fits.Column(name="FLAGS",format='D', array=n.empty_like(magAp_1a_g) )
c13 = fits.Column(name="FLAGS2",format='D', array=n.empty_like(magAp_1a_g) )

c14 = create_5_emptyCols("PSF_FWHM")
c15 = create_5_emptyCols("PSFFLUX")
c16 = create_5_emptyCols("PSFFLUX_IVAR")

c17 = create_5_emptyCols("EXTINCTION")
c18 = create_5_emptyCols("FIBERFLUX")

c19 = create_5_emptyCols("FIBERFLUX_IVAR")
c20 = create_5_emptyCols("FIBER2FLUX")

c21 = create_5_emptyCols("FIBER2FLUX_IVAR")
c22 = create_5_emptyCols("MODELFLUX")

c23 = create_5_emptyCols("MODELFLUX_IVAR")
c24 = create_5_emptyCols("MODELMAG")

c25 = create_5_emptyCols("MODELMAG_IVAR")

c26 = fits.Column(name="RESOLVE_STATUS",format='I', array=n.empty_like(magAp_1a_g) )
c27 = fits.Column(name="W1_MAG",format='D', array=n.empty_like(magAp_1a_g) )
c28 = fits.Column(name="W1_MAG_ERR",format='D', array=n.empty_like(magAp_1a_g) )
c29 = fits.Column(name="W1_NANOMAGGIES",format='D', array=n.empty_like(magAp_1a_g) )
c30 = fits.Column(name="W1_NANOMAGGIES_IVAR",format='D', array=n.empty_like(magAp_1a_g) )

c31 = fits.Column(name="W2_NANOMAGGIES",format='D', array=n.empty_like(magAp_1a_g) )
c32 = fits.Column(name="W2_NANOMAGGIES_IVAR",format='D', array=n.empty_like(magAp_1a_g) )
c33 = fits.Column(name="HAS_WISE_PHOT",format='L', array=n.empty_like(magAp_1a_g) )
c34 = fits.Column(name="OBJID_TARGETING",format='D', array=n.empty_like(magAp_1a_g) )
sourcetype=n.array([ "ELG_DECALS_TEST1" for i in range(len(magAp_1a_g)) ])
c35 = fits.Column(name="sourcetype",format='A16', array=sourcetype)

new_columns = hd[1].data.columns + c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 + c11 + c12 + c13 + c14 + c15 + c16 + c17 + c18 + c19 + c20 + c21 + c22 + c23 + c24 + c25 + c26 + c27 + c28 + c29 + c30 + c31 + c32 + c33 + c34 + c35 #+ c36 + c37 + c38 + c39 + c40 + c41 + c42 + c43 + c44 + c45 + c46 + c47 + c48 + c49 + c50 + c51 + c52 + c53 + c54 + c55 + c56 + c57 + c58 + c59 + c60 + c61 + c62 + c63 + c64 + c65 + c66 + c67 + c68 + c69 + c70 + c71 + c72 + c73 + c74 + c75 + c76 + c77 + c78 + c79 + c80 + c81 + c82 + c83 + c84 + c85 + c86 + c87 + c88 + c89 + c90 + c91 + c92 + c93 + c94 + c95 

hdu = fits.BinTableHDU.from_columns(new_columns)
#hdu.writeto("catalogs/allELG240-masked.fits")
hdu.writeto(outputCat)


