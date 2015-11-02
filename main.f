	program main
        implicit none   
	include 'para.inc'

c The flow of Image Processing:

c I. THELI (data reduction), including background removal, astrometry correction, source Identification. It provides reduced single exposures, head files for astrometry, a weight/flag image for bad pixels, and source catalog for stars/galaxies;

c I-I. Photo-z measurement;

c II. Image (star/galaxy) Confirmation, including removal of bad images. It provides image stamps and catalogs for the following three kinds of objects: 1) bright stars for PSF reconstruction; 2) galaxies; 3) noises.

c III. PSF reconstruction/Shear measurement. It produces a shear catalog.

c IV. Shear statistics is calculated. 

        character fname(20)*6
	character exponame(-4:4,-4:3,nexpo)*6,fieldname(-4:4,-4:3)*6
	common /field_info_pass/ exponame,fieldname

	call assign_field_info()

	call field_processing(-4,0)
	call field_processing(-3,3)
	call field_processing(-3,2)
	call field_processing(-3,-4)
	call field_processing(-3,-3)
	call field_processing(-3,-2)
	call field_processing(-2,3)
	call field_processing(-4,-3)
	call field_processing(-4,1)
	call field_processing(-4,2)
	call field_processing(-4,3)
	call field_processing(1,0)
	call field_processing(1,-2)
	call field_processing(1,-3)
	call field_processing(1,-4)
	call field_processing(1,1)
	call field_processing(2,0)
	call field_processing(2,-1)
	call field_processing(2,-2)
	call field_processing(2,-3)


	fname(1)=fieldname(-3,2)
	fname(2)=fieldname(-3,-4)
	fname(3)=fieldname(-3,-2)
	fname(4)=fieldname(-2,3)
	fname(5)=fieldname(-4,0)
	fname(6)=fieldname(-4,-3)
	fname(7)=fieldname(-4,1)
	fname(8)=fieldname(-4,2)
	fname(9)=fieldname(-4,3)
	fname(10)=fieldname(1,0)
	fname(11)=fieldname(1,-3)
	fname(12)=fieldname(1,-4)
	fname(13)=fieldname(2,0)
	fname(14)=fieldname(2,-1)
	fname(15)=fieldname(2,-2)
	fname(16)=fieldname(2,-3)
	fname(17)=fieldname(-3,-3)
	fname(18)=fieldname(-3,3)
	fname(19)=fieldname(1,-2)
	fname(20)=fieldname(1,1)

c	call shear_FD_test_multiple_field(20,fname)
c	call shear_FD_test_multiple_field_CFHTlens(20,fname)
	
	stop
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine field_processing(ixf,iyf)
        implicit none   
	include 'para.inc'

	character exponame(-4:4,-4:3,nexpo)*6,fieldname(-4:4,-4:3)*6
	common /field_info_pass/ exponame,fieldname

	integer i,j,k,u,v,ixf,iyf
	real gcat(ngfieldmax,nexpo,npara)
	integer ngcat(n_RA,n_Dec),mngcat(n_RA,n_Dec,nglimit,2),ngtot
	common /gcat_pass/ gcat,mngcat,ngcat,ngtot	

	do i=1,n_RA
	  do j=1,n_Dec
	    ngcat(i,j)=0
	    do k=1,nglimit
	      mngcat(i,j,k,1)=0
	      mngcat(i,j,k,2)=0
	    enddo
	  enddo
	enddo
	ngtot=0
	do i=1,ngfieldmax
	  do j=1,nexpo
	    do k=1,npara
	      gcat(i,j,k)=0
	    enddo
	  enddo
	enddo

	do i=1,nexpo
	  call read_info(fieldname(ixf,iyf),exponame(ixf,iyf,i))
 	  call measure_shear(fieldname(ixf,iyf),exponame(ixf,iyf,i))
	  call merge_exposures(fieldname(ixf,iyf),exponame(ixf,iyf,i))
	enddo

	
	call merger_shear_catalog(fieldname(ixf,iyf))
	call shrink_catalog(fieldname(ixf,iyf))

c	call w_SNR_relation(fieldname(ixf,iyf))


        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine assign_field_info()
	implicit none
	include 'para.inc'

	character exponame(-4:4,-4:3,nexpo)*6,fieldname(-4:4,-4:3)*6
	common /field_info_pass/ exponame,fieldname

	fieldname(2,-3)='w1p2m3'
	exponame(2,-3,1)='875613'
	exponame(2,-3,2)='875614'
	exponame(2,-3,3)='875616'
	exponame(2,-3,4)='875617'
	exponame(2,-3,5)='875618'
	exponame(2,-3,6)='875619'
	exponame(2,-3,7)='875620'

	fieldname(2,-2)='w1p2m2'
	exponame(2,-2,1)='827318'
	exponame(2,-2,2)='827319'
	exponame(2,-2,3)='827320'
	exponame(2,-2,4)='827321'
	exponame(2,-2,5)='827322'
	exponame(2,-2,6)='827323'
	exponame(2,-2,7)='827324'

	fieldname(2,-1)='w1p2m1'
	exponame(2,-1,1)='826729'
	exponame(2,-1,2)='826730'
	exponame(2,-1,3)='826731'
	exponame(2,-1,4)='826732'
	exponame(2,-1,5)='826733'
	exponame(2,-1,6)='826734'
	exponame(2,-1,7)='826735'

	fieldname(2,0)='w1p2m0'
	exponame(2,0,1)='816763'
	exponame(2,0,2)='816764'
	exponame(2,0,3)='816765'
	exponame(2,0,4)='816766'
	exponame(2,0,5)='816767'
	exponame(2,0,6)='816768'
	exponame(2,0,7)='816769'

	fieldname(1,1)='w1p1p1'
	exponame(1,1,1)='761475'
	exponame(1,1,2)='761476'
	exponame(1,1,3)='761477'
	exponame(1,1,4)='761478'
	exponame(1,1,5)='761479'
	exponame(1,1,6)='761480'
	exponame(1,1,7)='761481'

	fieldname(1,-4)='w1p1m4'
	exponame(1,-4,1)='963756'
	exponame(1,-4,2)='963757'
	exponame(1,-4,3)='963758'
	exponame(1,-4,4)='963759'
	exponame(1,-4,5)='963760'
	exponame(1,-4,6)='963761'
	exponame(1,-4,7)='963762'

	fieldname(1,-3)='w1p1m3'
	exponame(1,-3,1)='875443'
	exponame(1,-3,2)='875444'
	exponame(1,-3,3)='875445'
	exponame(1,-3,4)='875446'
	exponame(1,-3,5)='875447'
	exponame(1,-3,6)='875486'
	exponame(1,-3,7)='875487'

	fieldname(1,-2)='w1p1m2'
	exponame(1,-2,1)='831153'
	exponame(1,-2,2)='831154'
	exponame(1,-2,3)='831155'
	exponame(1,-2,4)='831156'
	exponame(1,-2,5)='831157'
	exponame(1,-2,6)='831158'
	exponame(1,-2,7)='831159'

	fieldname(1,0)='w1p1m0'
	exponame(1,0,1)='819869'
	exponame(1,0,2)='819870'
	exponame(1,0,3)='819871'
	exponame(1,0,4)='819872'
	exponame(1,0,5)='819873'
	exponame(1,0,6)='819874'
	exponame(1,0,7)='819875'

	fieldname(-4,3)='w1m4p3'
	exponame(-4,3,1)='881261'
	exponame(-4,3,2)='881262'
	exponame(-4,3,3)='881263'
	exponame(-4,3,4)='881264'
	exponame(-4,3,5)='881265'
	exponame(-4,3,6)='881266'
	exponame(-4,3,7)='881267'

	fieldname(-4,2)='w1m4p2'
	exponame(-4,2,1)='881254'
	exponame(-4,2,2)='881255'
	exponame(-4,2,3)='881256'
	exponame(-4,2,4)='881257'
	exponame(-4,2,5)='881258'
	exponame(-4,2,6)='881259'
	exponame(-4,2,7)='881260'

	fieldname(-4,1)='w1m4p1'
	exponame(-4,1,1)='880977'
	exponame(-4,1,2)='880978'
	exponame(-4,1,3)='880979'
	exponame(-4,1,4)='880980'
	exponame(-4,1,5)='880981'
	exponame(-4,1,6)='880982'
	exponame(-4,1,7)='880983'

	fieldname(-4,-3)='w1m4m3'
	exponame(-4,-3,1)='880193'
	exponame(-4,-3,2)='880194'
	exponame(-4,-3,3)='880195'
	exponame(-4,-3,4)='880196'
	exponame(-4,-3,5)='880197'
	exponame(-4,-3,6)='880198'
	exponame(-4,-3,7)='880199'


	fieldname(-4,0)='w1m4m0'
	exponame(-4,0,1)='880434'
	exponame(-4,0,2)='880436'
	exponame(-4,0,3)='880437'
	exponame(-4,0,4)='880438'
	exponame(-4,0,5)='880439'
	exponame(-4,0,6)='880440'
	exponame(-4,0,7)='880441'

	fieldname(-3,3)='w1m3p3'
	exponame(-3,3,1)='880101'
	exponame(-3,3,2)='880102'
	exponame(-3,3,3)='880103'
	exponame(-3,3,4)='880104'
	exponame(-3,3,5)='880105'
	exponame(-3,3,6)='880106'
	exponame(-3,3,7)='880107'


	fieldname(-3,2)='w1m3p2'
	exponame(-3,2,1)='880094'
	exponame(-3,2,2)='880095'
	exponame(-3,2,3)='880096'
	exponame(-3,2,4)='880097'
	exponame(-3,2,5)='880098'
	exponame(-3,2,6)='880099'
	exponame(-3,2,7)='880100'

	fieldname(-3,-4)='w1m3m4'
	exponame(-3,-4,1)='964850'
	exponame(-3,-4,2)='964851'
	exponame(-3,-4,3)='964852'
	exponame(-3,-4,4)='964853'
	exponame(-3,-4,5)='964854'
	exponame(-3,-4,6)='964855'
	exponame(-3,-4,7)='964856'

	fieldname(-3,-3)='w1m3m3'
	exponame(-3,-3,1)='879636'
	exponame(-3,-3,2)='879637'
	exponame(-3,-3,3)='879638'
	exponame(-3,-3,4)='879639'
	exponame(-3,-3,5)='879640'
	exponame(-3,-3,6)='879641'
	exponame(-3,-3,7)='879642'
c	exponame(8)='879643'


	fieldname(-3,-2)='w1m3m2'
	exponame(-3,-2,1)='879644'
	exponame(-3,-2,2)='879645'
	exponame(-3,-2,3)='879646'
	exponame(-3,-2,4)='879647'
	exponame(-3,-2,5)='879648'
	exponame(-3,-2,6)='879649'
	exponame(-3,-2,7)='879650'

	fieldname(-2,3)='w1m2p3'
	exponame(-2,3,1)='820409'
	exponame(-2,3,2)='820410'
	exponame(-2,3,3)='820411'
	exponame(-2,3,4)='820412'
	exponame(-2,3,5)='820413'
	exponame(-2,3,6)='820414'
	exponame(-2,3,7)='820415'

	fieldname(0,0)='w1m0m0'	
	exponame(0,0,1)='819719'
	exponame(0,0,2)='819720'
	exponame(0,0,3)='819721'
	exponame(0,0,4)='819722'
	exponame(0,0,5)='819723'
	exponame(0,0,6)='819724'
	exponame(0,0,7)='819725'
c	exponame(8)='816895'
c	exponame(9)='816896'
c	exponame(10)='816897'
	

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
