	subroutine measure_shear(fieldname,exponame)
        implicit none   
	include 'para.inc'

	character exponame*6,fieldname*6

	call prepare_PSF(fieldname,exponame)
	call process_galaxy(fieldname,exponame)

        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prepare_PSF(fieldname,exponame)
        implicit none
	include 'para.inc'
	
	character exponame*6,s*2,filename*80,fieldname*6

	integer nstar
	real star_collect(nstar_max,ns,ns),star_para(nstar_max,npara)
	common /star_pass/ star_collect,star_para,nstar

	real star_power(nstar_max,ns,ns)

	integer nPSF(nC)
	real PSF(nC,nstar_max,ns,ns),PSF_para(nC,nstar_max,npara)
	common /PSF_pass/ PSF,PSF_para,nPSF

	real star(ns,ns),power(ns,ns)	
	integer iC,i,j,u,v,nn,nn1,nn2

30      format(I2.2)

        do iC=1,36
          call read_star_catalog(fieldname,exponame,iC)

	  write(*,*) 'Exposure: ',fieldname,exponame
          write(*,*) 'Preparing stars on Chip ',iC,' with '
     .,nstar,' stars'

	  do i=1,nstar

            do u=1,ns
              do v=1,ns
                star(u,v)=star_collect(i,u,v)
              enddo
            enddo	

            call get_power(ns,ns,star,power)

            call smooth_image55ln(ns,ns,power,3)

	    call regularize_power(ns,ns,power)

            do u=1,ns
              do v=1,ns
                PSF(iC,i,u,v)=power(u,v)
                star_power(i,u,v)=power(u,v)
              enddo
            enddo	
	    do j=1,npara
	      PSF_para(iC,i,j)=star_para(i,j)
	    enddo

	  enddo
	  nPSF(iC)=nstar

          write(s,30) iC

          filename='../'//fieldname//'/step2/star_power'//exponame//'_'
     .//s//'.fits'

    	  nn1=ns*len_s
    	  nn2=ns*(int(nstar/len_s)+1)

 	  call write_stamps(nstar_max,1,nstar,ns,ns
     .,star_power,nn1,nn2,filename)

	enddo
 

	if (PSF_method.eq.1) call prepare_star_kriging(exponame)

	do iC=1,36
	  if (PSF_method.eq.1) then
	    call kriging_setup(iC)  
	  elseif (PSF_method.eq.2) then 
	    call poly_PSF_setup(iC)
	  elseif (PSF_method.eq.0) then
	    call mean_PSF_setup(iC)
	  endif
	enddo

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_PSF(iC,x,y,psf)
	implicit none
	include 'para.inc'

	integer iC
	real x,y,psf(ns,ns)

        if (PSF_method.eq.0) then
          call average_PSF(iC,x,y,psf)
        elseif (PSF_method.eq.1) then
          call interpolate_kriging_PSF(iC,x,y,psf)
        elseif (PSF_method.eq.2) then 
          call interpolate_poly_PSF(iC,x,y,psf)
        elseif (PSF_method.eq.3) then
          call interpolate_radial_PSF(iC,x,y,psf)
        endif
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
        subroutine process_galaxy(fieldname,exponame)
        implicit none
	include 'para.inc'
	
	character exponame*6,s*2,filename*80,fieldname*6
	integer ngal
	real gal_collect(ngal_max,ns,ns),gal_para(ngal_max,npara)
 	real noise_collect(ngal_max,ns,ns)
	common /gal_pass/ gal_collect,noise_collect,gal_para,ngal

	real gal(ns,ns),noise(ns,ns),psf(ns,ns),gf1,gf2,temp	
	real gal_p(ns,ns),noise_p(ns,ns),x,y,g1,g2,de	
	integer iC,i,j,u,v

	
        do iC=1,36
	  call read_gal_catalog(fieldname,exponame,iC)

	  write(*,*) 'Exposure: ',fieldname,exponame
          write(*,*) 'Processing galaxies on Chip ',iC,' with '
     .,ngal,' galaxies'

	  do i=1,ngal
            do u=1,ns
              do v=1,ns
                gal(u,v)=gal_collect(i,u,v)
                noise(u,v)=noise_collect(i,u,v)
              enddo
            enddo	
            call get_power(ns,ns,gal,gal_p)
            call get_power(ns,ns,noise,noise_p)

	    if (denoise_gal.eq.1) then
              call smooth_image55ln(ns,ns,gal_p,3)
              call smooth_image55ln(ns,ns,noise_p,3)
	    endif

	    call trim_power(ns,ns,gal_p)
	    call trim_power(ns,ns,noise_p)

	    if (norm_gal.eq.1) then
	      call normalize_gal(ns,ns,gal_p,noise_p)
	    endif

	    x=gal_para(i,1)
	    y=gal_para(i,2)
	    call get_PSF(iC,x,y,psf)
	    call get_shear(gal_p,noise_p,psf,g1,g2,de)
	    gal_para(i,17)=g1
	    gal_para(i,18)=g2
	    gal_para(i,19)=de

	    temp=2./(gal_para(i,13)+gal_para(i,16))
	    gf1=-(gal_para(i,13)*temp-1.)
	    gf2=-0.5*(gal_para(i,14)+gal_para(i,15))*temp
	    gal_para(i,20)=gf1
	    gal_para(i,21)=gf2

	  enddo

          write(s,30) iC
          filename='../'//fieldname//'/step2/shear_info'//exponame//'_'
     .//s//'.dat'
          open(unit=10,file=filename)
          rewind 10

  	  write(10,*) '#1 gal No.'
	  write(10,*) '#2 gal x position (pixel)'
	  write(10,*) '#3 gal y position (pixel)'
  	  write(10,*) '#4 gal RA'
	  write(10,*) '#5 gal Dec'
 	  write(10,*) '#6 gal SNR'
	  write(10,*) '#7 gal area'
	  write(10,*) '#8 gal flux'
	  write(10,*) '#9 gal Z_B'
	  write(10,*) '#10 gal flag'
	  write(10,*) '#11 gal e1'
	  write(10,*) '#12 gal e2'
	  write(10,*) '#13 gal weight'
	  write(10,*) '#14 field_distortion_matrix (1,1)'
	  write(10,*) '#15 field_distortion_matrix (1,2)'
	  write(10,*) '#16 field_distortion_matrix (2,1)'
	  write(10,*) '#17 field_distortion_matrix (2,2)'
	  write(10,*) '#18 g1'
	  write(10,*) '#19 g2'
  	  write(10,*) '#20 de'
	  write(10,*) '#21 g1_field'
  	  write(10,*) '#22 g2_field'

	  do i=1,ngal
	    write(10,*) i,gal_para(i,1),gal_para(i,2),gal_para(i,3)
     .,gal_para(i,4),gal_para(i,5),gal_para(i,6),gal_para(i,7)
     .,gal_para(i,8),gal_para(i,9),gal_para(i,10),gal_para(i,11)
     .,gal_para(i,12),gal_para(i,13),gal_para(i,14),gal_para(i,15)
     .,gal_para(i,16),gal_para(i,17),gal_para(i,18),gal_para(i,19)
     .,gal_para(i,20),gal_para(i,21)

	  enddo
  	  close(10)

	enddo

30      format(I2.2)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE get_shear(gal,noise,psf,g1,g2,de)
        implicit none
	include 'para.inc'

        real g1,g2,de
        real gal(ns,ns),noise(ns,ns),psf(ns,ns),fac(3,ns,ns)
        integer i,j,n_2,cc
	real ks,peak,thresh,area,ks_2,kx,kx2,ky,ky2,k2,k,temp	
	real filter_deriv,filter
	
	peak=psf(1,1)
	do i=1,ns
	  do j=1,ns
	    if (psf(i,j).gt.peak) peak=psf(i,j)
	  enddo
	enddo

        thresh=exp(-1.)*peak

        area=0.
        do i=1,ns
          do j=1,ns
	    if (psf(i,j).ge.thresh) area=area+1.
	  enddo
	enddo
	
	ks=sqrt(area/pi)	 
	ks_2=(ks*PSFr_ratio)**(-2)
	
	thresh=peak*1e-4
	n_2=ns/2
	cc=1+n_2

	do i=1,ns
	  kx=i-cc
	  kx2=kx*kx
	  do j=1,ns
	    ky=j-cc
	    ky2=ky*ky
	    k2=kx2+ky2
	    k=sqrt(k2)
	    if (psf(i,j).gt.thresh) then
	      temp=filter(k2,ks_2)/psf(i,j)
	      fac(1,i,j)=-temp*0.5*(kx2-ky2)
	      fac(2,i,j)=-temp*kx*ky
	      fac(3,i,j)
     .=k2*(temp+filter_deriv(k2,ks_2)*k*0.25/psf(i,j))
	    else
	      fac(1,i,j)=0.
	      fac(2,i,j)=0.
	      fac(3,i,j)=0.
	    endif
	  enddo
	enddo

        g1=0.
        g2=0.
        de=0.

        do i=1,ns
          do j=1,ns
            g1=g1+(gal(i,j)-noise(i,j))*fac(1,i,j)
            g2=g2+(gal(i,j)-noise(i,j))*fac(2,i,j)
            de=de+(gal(i,j)-noise(i,j))*fac(3,i,j)
          enddo
        enddo

        return
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function filter(k2,ks_2)
	implicit none
	
	real filter,k2,ks_2

	filter=exp(-k2*ks_2)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function filter_deriv(k2,ks_2)
	implicit none
	
	real filter_deriv,ks_2,k2

	filter_deriv=-ks_2*2.*sqrt(k2)*exp(-k2*ks_2)	

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

