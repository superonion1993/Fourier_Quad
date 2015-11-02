	subroutine read_info(fieldname,exponame)
        implicit none   
	include 'para.inc'

	character fieldname*6
	character exponame*6,weightname*80,gal_catname*80
	character filename*80,headname*80,catname*80
	integer nx,ny,i,j,iC

	real array(npx,npy)
	common /image_pass/ array,nx,ny
	real weight(npx,npy)
	common /weight_pass/ weight
	real mark(npx,npy)
	common /mark_pass/ mark

	integer nstep
	character s,ss*2

20     format(I1.1)
30     format(I2.2)

	gal_catname='../'//fieldname//'/cat/'//fieldname//'.dat'

	do iC=1,nC
	  if (iC.le.9) then
	    write(s,20) iC
	    filename='../'//fieldname//'/'//exponame//
     .'p_'//s//'OFCBC.fits'
	    weightname='../'//fieldname//'/'//exponame//
     .'p_'//s//'OFCBC.weight.fits'
	    headname='../'//fieldname//'/headers/'//exponame//
     .'p_'//s//'.head'
	    catname='../'//fieldname//'/cat/'//exponame//
     .'p_'//s//'OFCBC.star.asc'
	  else
	    write(ss,30) iC
	    filename='../'//fieldname//'/'//exponame//
     .'p_'//ss//'OFCBC.fits'
	    weightname='../'//fieldname//'/'//exponame//
     .'p_'//ss//'OFCBC.weight.fits'
	    headname='../'//fieldname//'/headers/'//exponame//
     .'p_'//ss//'.head'
	    catname='../'//fieldname//'/cat/'//exponame//
     .'p_'//ss//'OFCBC.star.asc'
	  endif

          call readimage(filename,nx,ny,npx,npy,array)
          call readimage(weightname,nx,ny,npx,npy,weight)
  	  call read_head(headname)

          call mark_image()
  	  call gen_star_catalog(fieldname,catname,exponame,iC)	 
	  call gen_gal_catalog(fieldname,gal_catname,exponame,iC)

	  write(*,*) 'Field: ',fieldname
	  write(*,*) 'Exposure: ',exponame
	  write(*,*) 'Chip No: ',iC

	enddo

        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine read_head(filename)
        implicit none   

	character filename*80
        real CRPIX(2),CD(2,2),CD_1(2,2),CRVAL(2)
	real PV(2,0:10)
	common /coor_pass/ PV,CD,CD_1,CRPIX,CRVAL

	real temp
	integer i,j,ierror
	character names*80,lit1*80

        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10

	do i=1,9	 
  	  read(10,*) 
	enddo

        read(10,*) names,lit1,CRVAL(1)
        read(10,*) names,lit1,CRVAL(2)
        read(10,*) names,lit1,CRPIX(1)
        read(10,*) names,lit1,CRPIX(2)

	do i=1,2
	  do j=1,2
    	    read(10,*) names,lit1,CD(i,j)
	  enddo 
	enddo
	do i=1,2
	  do j=0,2
    	    read(10,*) names,lit1,PV(i,j)
	  enddo 
	  PV(i,3)=0d0
	  do j=4,10
    	    read(10,*) names,lit1,PV(i,j)
	  enddo 
	enddo
	close(10)

	temp=CD(1,1)*CD(2,2)-CD(1,2)*CD(2,1)
	temp=1./temp

	CD_1(1,1)=CD(2,2)*temp
	CD_1(2,2)=CD(1,1)*temp
	CD_1(1,2)=-CD(1,2)*temp
	CD_1(2,1)=-CD(2,1)*temp

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mark_image()
	implicit none
	include 'para.inc'

	integer nx,ny
	real array(npx,npy)
	common /image_pass/ array,nx,ny
	real weight(npx,npy)
	common /weight_pass/ weight
	real mark(npx,npy)
	common /mark_pass/ mark

	real pix(np)	
	common /noise_sample_pass/ pix
	real med,sig
	common /noise_para_pass/ med,sig

	integer i,j,u,v,ix,iy,dx,dy,n1,n2,icount
	real ran1,temp

	icount=0
	do i=1,np
10	  ix=int(ran1()*(nx-1)+1)
	  iy=int(ran1()*(ny-1)+1)
	  if (weight(ix,iy).eq.0) then
	    icount=icount+1 
	    if (icount.ge.np) pause 'The masked area is too large!'
	    goto 10
	  endif
	  pix(i)=array(ix,iy)
	enddo
	call sort(np,np,pix)


	n1=np/6
	n2=(np*5)/6	

	med=pix(np/2)
	sig=(pix(n2)-pix(n1))*0.5

	do i=1,nx
	  do j=1,ny
	    if (weight(i,j).eq.0) then
              mark(i,j)=thresh_low
	    else
	      mark(i,j)=int((array(i,j)-med)/sig)
	      if (mark(i,j).le.thresh_low) then
	        mark(i,j)=thresh_low
	      elseif (mark(i,j).le.0) then
	        mark(i,j)=0
	      elseif (mark(i,j).ge.thresh_high) then
	        mark(i,j)=thresh_high
	      endif
	    endif
	  enddo
	enddo

	do i=1,nx
	  do j=1,ny
	    if (mark(i,j).eq.thresh_low) then
	      do u=max(1,i-2),min(nx,i+2)
	        do v=max(1,j-2),min(ny,j+2)
	          mark(u,v)=-5
	        enddo
	      enddo 
	    endif
	  enddo
	enddo

	do i=1,nx
	  do j=1,ny
	    if (mark(i,j).eq.-5) mark(i,j)=thresh_low
	  enddo
	enddo


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coordinate_transfer(a,d,x,y,direc)
	implicit none
	
	integer direc
	real a,d,x,y,xx,yy,xxx,yyy,rr
        real CRPIX(2),CD(2,2),CD_1(2,2),CRVAL(2)
	real PV(2,0:10)
	common /coor_pass/ PV,CD,CD_1,CRPIX,CRVAL
	real temp1,temp2,da,dd,const1,const2,ds
	real cosda,tandc,tandd,tanda
	real pi
	parameter (pi=3.1415926)

	const1=pi/180.
	tandc=tan(CRVAL(2)*const1)

	if (direc.eq.1) then
	  xx=CD(1,1)*(x-CRPIX(1))+CD(1,2)*(y-CRPIX(2))
	  yy=CD(2,1)*(x-CRPIX(1))+CD(2,2)*(y-CRPIX(2))	

	  xxx=PV(1,0)+PV(1,1)*xx+PV(1,2)*yy+PV(1,4)*xx**2
     .+PV(1,5)*xx*yy+PV(1,6)*yy**2+PV(1,7)*xx**3+PV(1,8)*xx**2*yy
     .+PV(1,9)*xx*yy**2+PV(1,10)*yy**3

	  yyy=PV(2,0)+PV(2,1)*yy+PV(2,2)*xx+PV(2,4)*yy**2
     .+PV(2,5)*xx*yy+PV(2,6)*xx**2+PV(2,7)*yy**3+PV(2,8)*yy**2*xx
     .+PV(2,9)*yy*xx**2+PV(2,10)*xx**3

	  xxx=xxx*const1
	  yyy=yyy*const1

	  da=xxx/(cos(CRVAL(2)*const1)*(1.-yyy*tandc))
	  da=da-da*da*da*0.33333333333
	  a=da/const1+CRVAL(1)
	  cosda=1.-da*da*0.5+da**4/24.

	  dd=(yyy*(cosda+tandc**2)+tandc*(cosda-1.))
     ./(yyy*tandc*(cosda-1.)+1.+cosda*tandc**2)	  
	  dd=dd-dd*dd*dd*0.3333333333	  	 
	  d=dd/const1+CRVAL(2)
 
	else
	
	  da=(a-CRVAL(1))*const1
	  dd=(d-CRVAL(2))*const1
	  tandd=tan(dd)
	  cosda=cos(da)
	  tanda=tan(da)

	  yyy=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
	  yyy=yyy/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
	  xxx=tanda*(cos(CRVAL(2)*const1)*(1.-yyy*tandc))

	  xxx=xxx/const1
	  yyy=yyy/const1

	  xx=xxx
	  yy=yyy

	  temp1=(xxx-PV(1,0)-PV(1,2)*yy-PV(1,4)*xx**2
     .-PV(1,5)*xx*yy-PV(1,6)*yy**2-PV(1,7)*xx**3
     .-PV(1,8)*xx**2*yy-PV(1,9)*xx*yy**2-PV(1,10)*yy**3)/PV(1,1)
	  temp2=(yyy-PV(2,0)-PV(2,2)*xx-PV(2,4)*yy**2
     .-PV(2,5)*xx*yy-PV(2,6)*xx**2-PV(2,7)*yy**3
     .-PV(2,8)*yy**2*xx-PV(2,9)*yy*xx**2-PV(2,10)*xx**3)/PV(2,1)
	  xx=temp1
	  yy=temp2

	  x=xx*CD_1(1,1)+yy*CD_1(1,2)+CRPIX(1)
	  y=xx*CD_1(2,1)+yy*CD_1(2,2)+CRPIX(2)

	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine field_distortion(x,y,distort_matrix)
	implicit none
	
	real x,y,xx,yy,xxx,yyy
        real CRPIX(2),CD(2,2),CD_1(2,2),CRVAL(2)
	real PV(2,0:10)
	common /coor_pass/ PV,CD,CD_1,CRPIX,CRVAL
	real distort_matrix(2,2)
	real dxx_dx,dxx_dy,dyy_dx,dyy_dy,temp
	real dxxx_dxx,dxxx_dyy,dyyy_dxx,dyyy_dyy

        xx=CD(1,1)*(x-CRPIX(1))+CD(1,2)*(y-CRPIX(2))
        yy=CD(2,1)*(x-CRPIX(1))+CD(2,2)*(y-CRPIX(2))	

	dxx_dx=CD(1,1)
        dyy_dx=CD(2,1)	
	dxx_dy=CD(1,2)
        dyy_dy=CD(2,2)	

        dxxx_dxx=PV(1,1)+2.*PV(1,4)*xx+PV(1,5)*yy+3.*PV(1,7)*xx**2
     .+2.*PV(1,8)*xx*yy+PV(1,9)*yy**2
        dxxx_dyy=PV(1,2)+PV(1,5)*xx+2.*PV(1,6)*yy+PV(1,8)*xx**2
     .+2.*PV(1,9)*xx*yy+3.*PV(1,10)*yy**2
	dyyy_dyy=PV(2,1)+2.*PV(2,4)*yy+PV(2,5)*xx+3.*PV(2,7)*yy**2
     .+2.*PV(2,8)*yy*xx+PV(2,9)*xx**2
	dyyy_dxx=PV(2,2)+PV(2,5)*yy+2.*PV(2,6)*xx+PV(2,8)*yy**2
     .+2.*PV(2,9)*yy*xx+3.*PV(2,10)*xx**2

	distort_matrix(1,1)=dxxx_dxx*dxx_dx+dxxx_dyy*dyy_dx
	distort_matrix(1,2)=dxxx_dxx*dxx_dy+dxxx_dyy*dyy_dy
	distort_matrix(2,1)=dyyy_dxx*dxx_dx+dyyy_dyy*dyy_dx
	distort_matrix(2,2)=dyyy_dxx*dxx_dy+dyyy_dyy*dyy_dy

	distort_matrix(1,1)=-distort_matrix(1,1)
	distort_matrix(1,2)=-distort_matrix(1,2)
	

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine draw_distortion(x,y,distort_matrix
     .,nx,ny,npx,npy,array)
	implicit none
	
	integer npx,npy
	integer nx,ny,i,j,ix,iy,u,v
	real array(npx,npy)	
	real x,y,distort_matrix(2,2),dx,dy
	real ksi,eta,phi,dphi,g1,g2,magn,temp
	real twopi,const,inten,exag,a,b
	parameter (twopi=3.1415926*2.)
	parameter (const=3e2*5e-5)
	parameter (inten=1e2)
	parameter (exag=20.)
	integer nstep,thick
	parameter (nstep=5000)
	parameter (thick=30)

	dphi=twopi/nstep

	temp=2./(distort_matrix(1,1)+distort_matrix(2,2))
	magn=const*temp

	g1=(distort_matrix(1,1)*temp-1.)*exag
	g2=0.5*(distort_matrix(2,1)+distort_matrix(1,2))*temp*exag

	
	do i=1,nstep
	  phi=dphi*i
	  ksi=cos(phi)
	  eta=sin(phi)

          a=((1.-g1)*ksi-g2*eta)*magn
          b=(-g2*ksi+(1.+g1)*eta)*magn
	  
	  do u=340,400
	    temp=u/400.	
	    dx=a*temp
	    dy=b*temp	  
	    ix=int(x+dx+0.5)
	    iy=int(y+dy+0.5)	 
            array(ix,iy)=inten
    	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine transfer_to_standard_coor(x,y,xs,ys)
	implicit none
	
	real x,y,xx,yy,xs,ys
        real CRPIX(2),CD(2,2),CD_1(2,2),CRVAL(2)
	real PV(2,0:10)
	common /coor_pass/ PV,CD,CD_1,CRPIX,CRVAL

        xx=CD(1,1)*(x-CRPIX(1))+CD(1,2)*(y-CRPIX(2))
        yy=CD(2,1)*(x-CRPIX(1))+CD(2,2)*(y-CRPIX(2))	

        xs=PV(1,0)+PV(1,1)*xx+PV(1,2)*yy+PV(1,4)*xx**2
     .+PV(1,5)*xx*yy+PV(1,6)*yy**2+PV(1,7)*xx**3+PV(1,8)*xx**2*yy
     .+PV(1,9)*xx*yy**2+PV(1,10)*yy**3

        ys=PV(2,0)+PV(2,1)*yy+PV(2,2)*xx+PV(2,4)*yy**2
     .+PV(2,5)*xx*yy+PV(2,6)*xx**2+PV(2,7)*yy**3+PV(2,8)*yy**2*xx
     .+PV(2,9)*yy*xx**2+PV(2,10)*xx**3

	xs=-xs

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
