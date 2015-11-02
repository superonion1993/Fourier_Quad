	real pi
	parameter (pi=3.1415926)

	integer NMAX
	parameter (NMAX=20000)

	integer npx,npy
	parameter (npx=3000)
	parameter (npy=5000)

	integer ns
	parameter (ns=48)
		
	integer len_g,len_s,nstar_max,ngal_max,npara
	parameter (len_g=40)
	parameter (len_s=15)
	parameter (ngal_max=len_g*len_g)
	parameter (nstar_max=len_s*len_s)
	parameter (npara=24)

	integer nschip,nC,np
	parameter (nschip=200)
	parameter (nC=36)
	parameter (np=5000)

	integer npcorr
	parameter (npcorr=nstar_max*nstar_max*nC)
	
	real SNR_thresh_star,saturation_thresh
	parameter (SNR_thresh_star=100)
	parameter (saturation_thresh=40000)

	real thresh_high,thresh_low
	parameter (thresh_high=4.)
	parameter (thresh_low=-4.)
	
	integer PSF_method
	parameter (PSF_method=2)

c PSF_method=0: Average PSF
c PSF_method=1: Kriging
c PSF_method=2: Polynomial Fitting
c PSF_method=3: Radial Weighting	
	
	integer norder
	parameter (norder=3)
	real rc,alpha
	parameter (rc=2000)
	parameter (alpha=-0.5)

	real PSFr_ratio
	parameter (PSFr_ratio=0.5)

	integer denoise_gal,norm_gal
	parameter (denoise_gal=0)
	parameter (norm_gal=0)

	integer margin
	parameter (margin=6)

c--------For Shear Field test----------------------------------
	integer nbin
	parameter (nbin=8)

c--------Merge shear catalog------------------------------------

	real dangle
	parameter (dangle=0.1)
	
	real w1_RA_min,w1_RA_max,w1_Dec_min,w1_Dec_max
	parameter (w1_RA_min=30.)
	parameter (w1_RA_max=39.)
	parameter (w1_Dec_min=-11.5)
	parameter (w1_Dec_max=-3.5)

	integer n_RA,n_Dec,nglimit,nexpo,ngfieldmax
	parameter (n_RA=(w1_RA_max-w1_RA_min)/dangle)
	parameter (n_Dec=(w1_Dec_max-w1_Dec_min)/dangle)
	parameter (nglimit=600)
	parameter (ngfieldmax=40000)
	parameter (nexpo=7)
c---------------------------------------------------------------

