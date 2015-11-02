        subroutine read_sextractor_cat()
        implicit none
        include 'para.inc'
	
        character filename*80
	integer ierror,i,j,nso	
	real x(NMAX,2),mag,rad,xx,yy,xc,yc

	nso=0
        filename='../W4m1p2/cat_asc/995283p.cat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10

        if (ierror.ne.0) then
	  pause 'sex cat reading error!!'
        endif

	do i=1,4
	  read(10,*)
	enddo

        do while (ierror.ge.0)
	  read(10,*,iostat=ierror) xx,yy,rad,mag
	  if (ierror.ge.0) then
            nso=nso+1
	    x(nso,1)=mag
	    x(nso,2)=rad
	  endif
	enddo
	close(10)

        filename='mag_r.fits'
	xc=14.5
	yc=1.9
	call map_points(NMAX,nso,x,xc,yc,500,500,filename)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine read_star_catalog(fieldname,exponame,iC)
        implicit none
        include 'para.inc'
	
        character exponame*6,s*2,filename*80,fieldname*6
	integer iC,ierror,aa,nstar,i,j,u,v,nn1,nn2,nn	
	real a(npara)
	real star_collect(nstar_max,ns,ns),star_para(nstar_max,npara)
	common /star_pass/ star_collect,star_para,nstar
 
30      format(I2.2)
        write(s,30) iC

	nstar=0
        filename='../'//fieldname//'/step1/star_info'//exponame//'_'
     .//s//'.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10

        if (ierror.ne.0) then
	  pause 'Star catalog reading error!!'
        endif

	do i=1,7
	  read(10,*)
	enddo

        do while (ierror.ge.0)
	  read(10,*,iostat=ierror) aa,a(1),a(2),a(3),a(4),a(5),a(6)
	  if (ierror.ge.0) then
            nstar=nstar+1
            do i=1,6
              star_para(nstar,i)=a(i)
            enddo
	  endif
	enddo
	close(10)
	
 	nn1=ns*len_s
 	nn2=ns*(int(nstar/len_s)+1)

	filename='../'//fieldname//'/step1/star_'//exponame//'_'
     .//s//'.fits'
	call read_stamps(nstar_max,1,nstar,ns,ns
     .,star_collect,nn1,nn2,filename)

	nn=nstar
	nstar=0
	do i=1,nn
	  if (star_para(i,6).eq.0) then
	    nstar=nstar+1
            do j=1,6
              star_para(nstar,j)=star_para(i,j)
            enddo
	    do u=1,ns
	      do v=1,ns
	        star_collect(nstar,u,v)=star_collect(i,u,v) 
	      enddo
	    enddo
	  endif
	enddo

c	call show_star_posi()

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine read_gal_catalog(fieldname,exponame,iC)
        implicit none
        include 'para.inc'
	
        character exponame*6,s*2,filename*80,fieldname*6
	integer iC,ierror,aa,i,j,u,v,nn,nn1,nn2	
	real a(npara)
	integer ngal
	real gal_collect(ngal_max,ns,ns),gal_para(ngal_max,npara)
 	real noise_collect(ngal_max,ns,ns)
	common /gal_pass/ gal_collect,noise_collect,gal_para,ngal


30      format(I2.2)
        write(s,30) iC

	ngal=0
        filename='../'//fieldname//'/step1/gal_info'//exponame//'_'
     .//s//'.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10

        if (ierror.ne.0) then
	  pause 'Gal catalog reading error!!'
        endif
	do i=1,17
	  read(10,*)
	enddo

        do while (ierror.ge.0)
	  read(10,*,iostat=ierror) aa,a(1),a(2),a(3),a(4),a(5),a(6)
     .,a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15),a(16)
	  if (ierror.ge.0) then
	    ngal=ngal+1
	    do i=1,npara
	      gal_para(ngal,i)=a(i)
	    enddo
	  endif
	enddo
	close(10)
	
 	nn1=ns*len_g
 	nn2=ns*(int(ngal/len_g)+1)

	filename='../'//fieldname//'/step1/gal_'//exponame//'_'
     .//s//'.fits'
	call read_stamps(ngal_max,1,ngal,ns,ns
     .,gal_collect,nn1,nn2,filename)

	filename='../'//fieldname//'/step1/noise'//exponame//'_'
     .//s//'.fits'
	call read_stamps(ngal_max,1,ngal,ns,ns
     .,noise_collect,nn1,nn2,filename)

	nn=ngal
	ngal=0
	do i=1,nn
	  if (gal_para(i,9).eq.0) then
	    ngal=ngal+1
            do j=1,npara
              gal_para(ngal,j)=gal_para(i,j)
            enddo
	    do u=1,ns
	      do v=1,ns
	        gal_collect(ngal,u,v)=gal_collect(i,u,v) 
	        noise_collect(ngal,u,v)=noise_collect(i,u,v) 
	      enddo
	    enddo
	  endif
	enddo

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine show_star_posi()
        implicit none
        include 'para.inc'
	
        character filename*80
	integer nstar,i,j,u,v,x,y,off	
	integer mag,nx,ny
	parameter (mag=4)
	parameter (nx=2400/mag)
	parameter (ny=4800/mag)
	real star_collect(nstar_max,ns,ns),star_para(nstar_max,npara)
	common /star_pass/ star_collect,star_para,nstar

	real map1(nx,ny),map2(nx,ny),map3(nx,ny),ave(ns,ns)
	real star(ns,ns),power(nstar_max,ns,ns),resi(ns,ns)
	real pp(ns,ns)

	do i=1,nx
	  do j=1,ny
	    map1(i,j)=0.
	    map2(i,j)=0.
	    map3(i,j)=0.
	  enddo
	enddo

	do i=1,ns
	  do j=1,ns
	    ave(i,j)=0.
	  enddo
	enddo

	off=ns/2 	
	do i=1,nstar
          x=int(star_para(i,1)/mag+0.5)
          y=int(star_para(i,2)/mag+0.5)
          do u=1,ns
            do v=1,ns
	      star(u,v)=star_collect(i,u,v) 
	      map1(x-off+u,y-off+v)=star_collect(i,u,v) 
	    enddo
	  enddo
          call get_power(ns,ns,star,pp)
          call smooth_image55ln(ns,ns,pp,3)
	  call regularize_power(ns,ns,pp)	  
          do u=1,ns
            do v=1,ns
	      power(i,u,v)=pp(u,v)
	      ave(u,v)=ave(u,v)+pp(u,v)
	      map2(x-off+u,y-off+v)=pp(u,v) 
	    enddo
	  enddo
	enddo
        do u=1,ns
          do v=1,ns
            ave(u,v)=ave(u,v)/nstar
          enddo
        enddo
	do i=1,nstar
          x=int(star_para(i,1)/mag+0.5)
          y=int(star_para(i,2)/mag+0.5)
          do u=1,ns
            do v=1,ns
	      map3(x-off+u,y-off+v)=power(i,u,v)-ave(u,v) 
	    enddo
	  enddo
	enddo

        filename='star_stamp_posi.fits'
	call writeimage(filename,nx,ny,nx,ny,map1)
        filename='star_power_posi.fits'
	call writeimage(filename,nx,ny,nx,ny,map2)
        filename='star_resi_posi.fits'
	call writeimage(filename,nx,ny,nx,ny,map3)

	pause 'star_posi generated!'

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
