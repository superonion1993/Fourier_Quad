        subroutine gen_gal_catalog(fieldname,catname,exponame,iC)
        implicit none
        include 'para.inc'
	
        character catname*80,exponame*6,s*2,filename*80,fname*6
	character fieldname*6
	integer iC,ierror,flag,counts(0:6),nflag

        real ra,dec,FWHMr,Kron_r,x,y,dm(2,2),flagg
     .,flux_r,class_star,PSFe1,PSFe2,e1,e2,weight,SNRr,MASK,Z_B
     .,fitclass,mmm,c2,LP_Mi,MAG_i,star_flag	
	
        integer ngal,i,j,x1,x2,y1,y2,u,v,halfns,nbad,nn,nn1,nn2
        real gal_collect(ngal_max,ns,ns),gal_para(ngal_max,npara)
	real noise_collect(ngal_max,ns,ns),bad_collect(ngal_max,ns,ns)
        real gal(ns,ns),noise(ns,ns)     

	real peak,area,flux,SNR,FWHM
	common /stamp_pass/ peak,area,flux,SNR,FWHM


c flag=0: No problem detected;
c flag=1: stamp out of the image boundary;
c flag=2: stamp center is not bright enough or too bright;
c flag=3: source is too large or too small;
c flag=4: multiple defects within source region;
c flag=5: multiple sources within the stamp;
c flag=6: cannot find noise stamp.

        do i=0,6
          counts(i)=0
        enddo
	
        halfns=ns/2	  
 
        open(unit=10,file=catname,status='old',iostat=ierror)
        rewind 10

        if (ierror.ne.0) then
	  pause 'Catalog file error!!'
        endif
   
        read(10,*) 

        ngal=0
	nbad=0

        do while (ierror.ge.0)
c          read(10,*,iostat=ierror) fname,ra,dec,flag,FWHMr
c     .,Kron_r,flux_r,class_star,PSFe1,PSFe2,e1,e2,weight,SNRr,MASK,Z_B	

          read(10,*,iostat=ierror) ra,dec,flagg,flux_r
     .,e1,e2,weight,fitclass,SNRr,MASK,Z_B,mmm,c2,LP_Mi,star_flag,MAG_i

 	  if (ierror.ge.0) then
	    call coordinate_transfer(ra,dec,x,y,0)
	    call field_distortion(x,y,dm)
	    x1=int(x-halfns) 
  	    x2=x1+ns-1
	    y1=int(y-halfns)
	    y2=y1+ns-1

            call check_source(1,x1,x2,y1,y2,flag,gal)

            if (flag.eq.0) then
              call find_noise(x1,x2,y1,y2,noise,nflag)
              if (nflag.eq.0) then
	        ngal=ngal+1
                do u=1,ns
                  do v=1,ns
                    gal_collect(ngal,u,v)=gal(u,v)
                    noise_collect(ngal,u,v)=noise(u,v)
                  enddo
                enddo
 	        gal_para(ngal,1)=x 
	        gal_para(ngal,2)=y
	        gal_para(ngal,3)=ra
	        gal_para(ngal,4)=dec	  
  	        gal_para(ngal,5)=SNR	  
  	        gal_para(ngal,6)=area	  
  	        gal_para(ngal,7)=flux	  
  	        gal_para(ngal,8)=Z_B	  
	        gal_para(ngal,9)=flagg	
	        gal_para(ngal,10)=e1	
	        gal_para(ngal,11)=e2	
	        gal_para(ngal,12)=weight	
	        gal_para(ngal,13)=dm(1,1)	  
	        gal_para(ngal,14)=dm(1,2)  
	        gal_para(ngal,15)=dm(2,1)	  
	        gal_para(ngal,16)=dm(2,2)  
	      else
	        flag=6
              endif
	    endif

            if (flag.ne.1.and.flag.ne.2.and.flag.ne.0) then
              nbad=nbad+1	
              do u=1,ns
                do v=1,ns
                  bad_collect(nbad,u,v)=gal(u,v)
                enddo
              enddo
            endif
            counts(flag)=counts(flag)+1
	  endif	 
        enddo

	write(*,*) 'Fine gals:',counts(0)
	write(*,*) 'Not bright enough or too bright:',counts(2)
	write(*,*) 'Too large or too small:',counts(3)
	write(*,*) 'Multiple defects:',counts(4)
	write(*,*) 'Noise not found:',counts(6)

        close(10)

30      format(I2.2)
        write(s,30) iC

	nn1=ns*len_g
 	nn2=ns*(int(ngal/len_g)+1)

	filename='../'//fieldname//'/step1/gal_'//exponame//'_'
     .//s//'.fits'
	call write_stamps(ngal_max,1,ngal,ns,ns
     .,gal_collect,nn1,nn2,filename)

	filename='../'//fieldname//'/step1/noise'//exponame//'_'
     .//s//'.fits'
	call write_stamps(ngal_max,1,ngal,ns,ns
     .,noise_collect,nn1,nn2,filename)

 	nn2=ns*(int(nbad/len_g)+1)

	filename='../'//fieldname//'/step1/bad_gal_'//exponame//'_'
     .//s//'.fits'
	call write_stamps(ngal_max,1,nbad,ns,ns
     .,bad_collect,nn1,nn2,filename)

	
        filename='../'//fieldname//'/step1/gal_info'//exponame//'_'
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

	do i=1,ngal
	  write(10,*) i,gal_para(i,1),gal_para(i,2),gal_para(i,3)
     .,gal_para(i,4),gal_para(i,5),gal_para(i,6),gal_para(i,7)
     .,gal_para(i,8),gal_para(i,9),gal_para(i,10),gal_para(i,11)
     .,gal_para(i,12),gal_para(i,13),gal_para(i,14),gal_para(i,15)
     .,gal_para(i,16)
	enddo
	close(10)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine gen_star_catalog(fieldname,catname,exponame,iC)
        implicit none
	include 'para.inc'
	
        character catname*80,exponame*6,s*2,filename*80,fieldname*6
	real star_collect(nstar_max,ns,ns),star_para(nstar_max,npara)
        real ra,dec,x,y,a,b,c,star(ns,ns),mag,flux_r
	real bad_collect(nstar_max,ns,ns)
      	integer x1,x2,y1,y2,u,v,flag,nstar,i,j,ierror,counts(0:6)
	integer iC,nbad,nn,nn1,nn2
	real ra_d,dec_d,x_d,y_d

	real peak,area,flux,SNR,FWHM
	common /stamp_pass/ peak,area,flux,SNR,FWHM

	do i=0,6
	  counts(i)=0
	enddo

        open(unit=10,file=catname,status='old',iostat=ierror)
        rewind 10

        if (ierror.ne.0) then
	  pause 'Catalog file error!!'
        endif
    
	read(10,*,iostat=ierror)
	nstar=0
	nbad=0
        do while (ierror.ge.0)
c          read(10,*,iostat=ierror) ra,dec,x,y,a,b,c
          read(10,*,iostat=ierror) x,y,flux_r,mag
      	
	  if (ierror.ge.0) then
	    x1=int(x-ns/2)
	    x2=x1+ns-1
	    y1=int(y-ns/2)
	    y2=y1+ns-1

	    call check_source(2,x1,x2,y1,y2,flag,star)
	
	    counts(flag)=counts(flag)+1

	    if (flag.ne.1) then
  	      if (flag.eq.0) then
                nstar=nstar+1	
 	        star_para(nstar,1)=x
	        star_para(nstar,2)=y
	        x_d=x
	        y_d=y
	        call coordinate_transfer(ra_d,dec_d,x_d,y_d,1)
	        star_para(nstar,3)=ra_d
	        star_para(nstar,4)=dec_d	  
	        star_para(nstar,5)=SNR	  
	        star_para(nstar,6)=flag	  
                do u=1,ns
                  do v=1,ns
                    star_collect(nstar,u,v)=star(u,v)
                  enddo
                enddo
	      else
	        nbad=nbad+1
                do u=1,ns
                  do v=1,ns
                    bad_collect(nbad,u,v)=star(u,v)
                  enddo
                enddo
	      endif
	    endif
	  endif	
        enddo

	write(*,*) 'Fine stars:',counts(0)
	write(*,*) 'Not bright enough or too bright:',counts(2)
	write(*,*) 'Too large:',counts(3)
	write(*,*) 'Multiple defects:',counts(4)

        close(10)

30      format(I2.2)
        write(s,30) iC

	nn1=ns*len_s
 	nn2=ns*(int(nstar/len_s)+1)

	filename='../'//fieldname//'/step1/star_'//exponame//'_'
     .//s//'.fits'
	call write_stamps(nstar_max,1,nstar
     .,ns,ns,star_collect,nn1,nn2,filename)

 	nn2=ns*(int(nbad/len_s)+1)

	filename='../'//fieldname//'/step1/bad_star_'//exponame//'_'
     .//s//'.fits'
	call write_stamps(nstar_max,1,nbad
     .,ns,ns,bad_collect,nn1,nn2,filename)
	
        filename='../'//fieldname//'/step1/star_info'//exponame//'_'
     .//s//'.dat'
        open(unit=10,file=filename)
        rewind 10
	write(10,*) '#1 star No.'
	write(10,*) '#2 star x position (pixel)'
	write(10,*) '#3 star y position (pixel)'
	write(10,*) '#4 star RA'
	write(10,*) '#5 star Dec'
	write(10,*) '#6 star SNR'
	write(10,*) '#7 star flag'
	do i=1,nstar
	  write(10,*) i,star_para(i,1),star_para(i,2),star_para(i,3)
     .,star_para(i,4),star_para(i,5),star_para(i,6)
	enddo
	close(10)
	
        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine find_noise(x1,x2,y1,y2,stamp,flag)
	implicit none
	include 'para.inc'

	integer x1,x2,y1,y2,flag
	real stamp(ns,ns),mimim(ns,ns)

	integer nx,ny
	real array(npx,npy)
	common /image_pass/ array,nx,ny
	real mark(npx,npy)
	common /mark_pass/ mark

	integer i,j,ix1,iy1,ix2,iy2,nd,u,v,jx,jy,n1,n2,k
	real peak,pmax,temp
	
	temp=100000.
	pmax=temp

	do i=1,3
	  do j=1,3
	    if (i.ne.2.or.j.ne.2) then
	      ix1=x1+ns*(i-2)
	      iy1=y1+ns*(j-2)
	      ix2=ix1+ns-1
	      iy2=iy1+ns-1
              if (ix1.lt.1.or.ix2.gt.nx.or.iy1.lt.1.or.iy2.gt.ny) 
     . goto 10	     
	      nd=0
	      peak=array(ix1,iy1)
	      do u=ix1,ix2
	        do v=iy1,iy2
	          if (mark(u,v).eq.thresh_low) then
                    nd=nd+1
	          else
	            if (array(u,v).gt.peak) peak=array(u,v)
	          endif
	        enddo
	      enddo
	      if (nd.le.10) then
	        if (peak.lt.pmax) then
	          pmax=peak
	          jx=ix1
	          jy=iy1
	        endif
	      endif
	    endif
10	  enddo
	enddo

	if (pmax.lt.temp) then

	  do i=jx,jx+ns-1
	    u=i-jx+1
	    do j=jy,jy+ns-1
	      v=j-jy+1
              stamp(u,v)=array(i,j)
	      mimim(u,v)=mark(i,j)
	    enddo
	  enddo
	
	  call decorate_stamp(mimim,stamp)

	  flag=0
	else
	  flag=1
	endif
	  
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_stamp_info(mark,stamp)
	implicit none
	include 'para.inc'

	real mark(ns,ns),stamp(ns,ns),thresh

	real med,sig
	common /noise_para_pass/ med,sig

	real peak,area,flux,SNR,FWHM
	common /stamp_pass/ peak,area,flux,SNR,FWHM

	integer i,j,xp,yp

	peak=-100.
	do i=1,ns
	  do j=1,ns
	    if (stamp(i,j).gt.peak) then
	      peak=stamp(i,j)
	      xp=i
	      yp=j
	    endif
	  enddo
	enddo

	thresh=peak*0.5
	
	flux=0.
	area=0.
	do i=1,ns
	  do j=1,ns
	    if (stamp(i,j).ge.thresh) then
	      flux=flux+stamp(i,j)
	      area=area+1.
	    endif
	  enddo
	enddo
	
	SNR=flux/(sig*sqrt(area))
	FWHM=2.*sqrt(area/pi)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine count_peaks(mark,stamp,bump_ratio)
	implicit none
	include 'para.inc'

	real mark(ns,ns),stamp(ns,ns),map(ns,ns),bump_ratio
	real peak_2nd,peak	
	integer i,j

	do i=1,ns
	  do j=1,ns
            map(i,j)=stamp(i,j)
	  enddo
	enddo

	call smooth_image55(ns,ns,map)
	call smooth_image55(ns,ns,map)
c	call smooth_image55(ns,ns,map)

	peak_2nd=0.
	peak=0.

	do i=1,ns
	  do j=1,ns
	    if (mark(i,j).eq.2.5.and.map(i,j).gt.map(i-1,j)
     ..and.map(i,j).gt.map(i+1,j).and.map(i,j).gt.map(i-1,j+1)
     ..and.map(i,j).gt.map(i,j+1).and.map(i,j).gt.map(i+1,j+1)
     ..and.map(i,j).gt.map(i-1,j-1).and.map(i,j).gt.map(i,j-1)
     ..and.map(i,j).gt.map(i+1,j-1)) then
	      if (map(i,j).gt.peak) then
	        peak_2nd=peak
	        peak=map(i,j)
	      elseif (map(i,j).le.peak
     ..and.map(i,j).gt.peak_2nd) then
	        peak_2nd=map(i,j)
	      endif        
	    endif
	  enddo
	enddo

        bump_ratio=peak_2nd/peak

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine decorate_stamp(mark,stamp)
	implicit none
	include 'para.inc'

	real mark(ns,ns),stamp(ns,ns)
	real pix(np)	
	common /noise_sample_pass/ pix
	real ran1
	integer n1,n2,i,j,k,u,v

	real aa,bb,cc,arr(ns*ns,3)
	integer nr

	nr=0
	do i=1,ns
	  do j=1,ns
	    if (mark(i,j).ge.-2.and.mark(i,j).le.2) then
              nr=nr+1
              arr(nr,1)=i
              arr(nr,2)=j
              arr(nr,3)=stamp(i,j)	    
	    endif
	  enddo
	enddo
	
	call find_slope_2D(ns*ns,nr,arr,aa,bb,cc)

	n1=np/10
	n2=np*9/10

	do i=1,ns
	  do j=1,ns
	    if (mark(i,j).eq.thresh_low) then
	      k=int(n1+(n2-n1)*ran1())
	      stamp(i,j)=pix(k)
 	    elseif (mark(i,j).ge.3) then
	      do u=max(i-1,1),min(i+1,ns)
	        do v=max(j-1,1),min(j+1,ns)
	          k=int(n1+(n2-n1)*ran1())
	          stamp(u,v)=pix(k)
	        enddo
	      enddo
	    else
  	      stamp(i,j)=stamp(i,j)-aa-bb*i-cc*j
	    endif
	  enddo
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine check_source(stype,x1,x2,y1,y2,flag,stamp)
	implicit none
	include 'para.inc'

	integer x1,x2,y1,y2,flag,stype
	real stamp(ns,ns)

	integer nx,ny
	real array(npx,npy)
	common /image_pass/ array,nx,ny
	real mark(npx,npy)
	common /mark_pass/ mark

	real med,sig
	common /noise_para_pass/ med,sig

	real mimim(ns,ns),aa,bb,cc,temp,nt,bump_ratio
	integer i,j,k,xc,yc,offx,offy,changed,u,v,n1,n2,ncc,nd,sarea

	real peak,area,flux,SNR,FWHM
	common /stamp_pass/ peak,area,flux,SNR,FWHM

c flag=0: No problem detected;
c flag=1: stamp out of the image boundary;
c flag=2: stamp center is not bright enough or too bright;
c flag=3: source is too large or too small;
c flag=4: multiple defects within source region;
c flag=5: multiple sources within the stamp.

	flag=0

	do i=1,ns
	  do j=1,ns
	    stamp(i,j)=0.
	  enddo
	enddo

c/////////To check if the stamp is out of the image boundary\\\\\\\\\\\	
        if (x1.lt.1.or.x2.gt.nx.or.y1.lt.1.or.y2.gt.ny) then
	  flag=1
	  return
	endif


c///////// Assign the image to the stamp\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	do i=1,ns
	  do j=1,ns
	    stamp(i,j)=array(i+x1-1,j+y1-1)
	  enddo
	enddo

c/////////To check if the stamp center is bright enough\\\\\\\\\\\\\\\

	xc=(x1+x2)/2
	yc=(y1+y2)/2
	if (mark(xc,yc).lt.3) then
	  flag=2
	  return
	endif

c/////////To identify the source region\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	offx=x1-1
	offy=y1-1
	do i=1,ns
	  do j=1,ns
	    mimim(i,j)=mark(i+offx,j+offy)
	  enddo
	enddo	
	mimim(xc-offx,yc-offy)=2.5
	changed=1
	sarea=1
	do while (changed.eq.1)	 
	  changed=0
   	  do i=1,ns
	    do j=1,ns
	      if (mimim(i,j).eq.2.5) then
	        do u=max(i-2,1),min(i+2,ns)
	          do v=max(j-2,1),min(j+2,ns)
	            if (mimim(u,v).ge.3) then
	              mimim(u,v)=2.5
	              changed=1
	              sarea=sarea+1
	            endif
	          enddo
	        enddo
	      endif
	    enddo
	  enddo
	enddo

c/////////To check if the source region is too small or too large (touching the stamp boundary)\\\\\\\\\\

	if (sarea.lt.9) then
	  flag=3
	  return
	endif

	do i=1,ns
	  do j=1,margin
	    if (mimim(i,j).eq.2.5.or.mimim(i,ns-j+1).eq.2.5) then
	      flag=3
	      return
	    endif
	  enddo
	enddo

	do j=1,ns
	  do i=1,margin
 	    if (mimim(i,j).eq.2.5.or.mimim(ns-i+1,j).eq.2.5) then
	      flag=3
	      return
	    endif
	  enddo
	enddo

c/////////To check if the defects touch the source region\\\\\\\\\\\\\\\\\

        do i=1,ns
          do j=1,ns
            if (mimim(i,j).eq.thresh_low) then
	      ncc=0
	      nd=0
	      nt=0
	      temp=0
              do u=max(i-1,1),min(i+1,ns)
                do v=max(j-1,1),min(j+1,ns)
	          if (u.ne.i.or.v.ne.j) then
                    if (mimim(u,v).eq.2.5) ncc=ncc+1
                    if (mimim(u,v).eq.thresh_low) nd=nd+1
	            nt=nt+1
	            temp=temp+stamp(u,v)
	          endif
                enddo
              enddo
	      if (ncc.gt.0.and.nd.gt.0) then
	        flag=4
	        return
	      endif
	      if (nd.eq.0) then
	        stamp(i,j)=temp/nt
                if (ncc.gt.0) then
	          mimim(i,j)=2.5
	        else
	          mimim(i,j)=int((stamp(i,j)-med)/sig)
	          if (mimim(i,j).le.thresh_low) then
	            mimim(i,j)=thresh_low
	          elseif (mimim(i,j).le.0) then
	            mimim(i,j)=0
	          elseif (mimim(i,j).ge.thresh_high) then
	            mimim(i,j)=thresh_high
	          endif
	        endif	     
	      endif
            endif
          enddo
        enddo

c//////////To decorate the defects \\\\\\\\\\\\\\

	call decorate_stamp(mimim,stamp)

c//////////Get information (SNR, FWHM, ...) of the source\\\\\\\\\\\\\\

	call get_stamp_info(mimim,stamp)

	if (peak.gt.saturation_thresh) then
	  flag=2
	  return
	endif

	if (stype.eq.2) then
	  if (SNR.lt.SNR_thresh_star) then
            flag=3
            return
	  else	     
	    call count_peaks(mimim,stamp,bump_ratio)
	    if (bump_ratio.ge.0.02) then
              flag=5
              return
	    endif	
	  endif
	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

