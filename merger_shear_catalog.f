        subroutine merger_shear_catalog(fieldname)
        implicit none
        include 'para.inc'

	real gcat(ngfieldmax,nexpo,npara)
	integer ngcat(n_RA,n_Dec),mngcat(n_RA,n_Dec,nglimit,2),ngtot
	common /gcat_pass/ gcat,mngcat,ngcat,ngtot	

	integer i,j,k,u,v
        character filename*80,fieldname*6
	
        filename='./step3/shear_'//fieldname//'.dat'
        open(unit=10,file=filename)
        rewind 10
	do i=1,n_RA
	  do j=1,n_Dec
	    if (ngcat(i,j).ge.1) then
	      do k=1,ngcat(i,j)
	        do u=1,mngcat(i,j,k,2)
                write(10,*) i,j,k,u,gcat(mngcat(i,j,k,1),u,1)
     .,gcat(mngcat(i,j,k,1),u,2),gcat(mngcat(i,j,k,1),u,3)
     .,gcat(mngcat(i,j,k,1),u,4),gcat(mngcat(i,j,k,1),u,5)
     .,gcat(mngcat(i,j,k,1),u,6),gcat(mngcat(i,j,k,1),u,7)
     .,gcat(mngcat(i,j,k,1),u,8),gcat(mngcat(i,j,k,1),u,9)
     .,gcat(mngcat(i,j,k,1),u,10),gcat(mngcat(i,j,k,1),u,11)
     .,gcat(mngcat(i,j,k,1),u,12),gcat(mngcat(i,j,k,1),u,13)
     .,gcat(mngcat(i,j,k,1),u,14)
	        enddo
	      enddo
	    endif
	  enddo
	enddo
	close(10)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shrink_catalog(fieldname)
        implicit none
        include 'para.inc'

	real gcat(nexpo,npara),cat(npara),a(npara)	
	integer i,j,k,u,v,iexpo,ngal,ii,jj,ierror
        character filename*80,fieldname*6

        filename='./step3/final_shear_'//fieldname//'.dat'
        open(unit=20,file=filename)
        rewind 20
	
	ngal=0
	iexpo=0

        filename='./step3/shear_'//fieldname//'.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          pause 'field catalog reading error!!'
        endif
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) i,j,k,u,a(1),a(2),a(3),a(4),a(5)
     .,a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14)
          if (ierror.ge.0) then	    
	    if (u.eq.1.and.iexpo.gt.0) then
	      call collapse_cat(iexpo,gcat,cat)
	      write(20,*) ii,jj,cat(1),cat(2),cat(3),cat(4),cat(5)
     .,cat(6),cat(7),cat(8),cat(9),cat(10),cat(11),cat(12),cat(13)
     .,cat(14)
	      iexpo=0
	      ngal=ngal+1
	    endif
	    iexpo=iexpo+1
            do k=1,npara
	      gcat(iexpo,k)=a(k)
	    enddo
	    ii=i	
	    jj=j
	  endif
	enddo
	close(10)

        call collapse_cat(iexpo,gcat,cat)
        write(20,*) ii,jj,cat(1),cat(2),cat(3),cat(4),cat(5)
     .,cat(6),cat(7),cat(8),cat(9),cat(10),cat(11),cat(12),cat(13)
     .,cat(14)
        ngal=ngal+1
	
	close(20)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine collapse_cat(iexpo,gcat,cat)
	implicit none
        include 'para.inc'

	real gcat(nexpo,npara),cat(npara)
	integer iexpo,i
	real func_weight,w

	cat(1)=gcat(1,1)
	cat(2)=gcat(1,2)
	cat(3)=0.
	cat(4)=0.
	cat(5)=0.
	cat(6)=gcat(1,6)
	cat(7)=gcat(1,7)
	cat(8)=gcat(1,8)
	cat(9)=gcat(1,9)
	cat(10)=0.
	cat(11)=0.
	cat(12)=0.
	cat(13)=0.
	cat(14)=0.
	do i=1,iexpo
	  cat(3)=cat(3)+gcat(i,3)
	  cat(4)=cat(4)+gcat(i,4)
	  cat(5)=cat(5)+gcat(i,5)
	  w=func_weight(gcat(i,3))
	  cat(10)=cat(10)+gcat(i,10)*w
	  cat(11)=cat(11)+gcat(i,11)*w
	  cat(12)=cat(12)+gcat(i,12)*w
	  cat(13)=cat(13)+gcat(i,13)*gcat(i,12)*w
	  cat(14)=cat(14)+gcat(i,14)*gcat(i,12)*w
	enddo
	cat(13)=cat(13)/cat(12)
	cat(14)=cat(14)/cat(12)
        cat(3)=cat(3)/iexpo
        cat(4)=cat(4)/iexpo
        cat(5)=cat(5)/iexpo
        cat(10)=cat(10)/iexpo
        cat(11)=cat(11)/iexpo
        cat(12)=cat(12)/iexpo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine merge_exposures(fieldname,exponame)
        implicit none
        include 'para.inc'
	
        character exponame*6,s*2,filename*80,fieldname*6
	integer iC,ierror,aa,i,nr,nd,found
	real a(npara)

	real gcat(ngfieldmax,nexpo,npara)
	integer ngcat(n_RA,n_Dec),mngcat(n_RA,n_Dec,nglimit,2),ngtot
	common /gcat_pass/ gcat,mngcat,ngcat,ngtot	
	
	do iC=1,36
          write(s,30) iC
          filename='../'//fieldname//'/step2/shear_info'//exponame//'_'
     .//s//'.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
	    pause 'Gal catalog reading error!!'
          endif
	  do i=1,22
	    read(10,*)
	  enddo
          do while (ierror.ge.0)
  	    read(10,*,iostat=ierror) aa,a(1),a(2),a(3),a(4),a(5)
     .,a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15)
     .,a(16),a(17),a(18),a(19),a(20),a(21)

	    if (ierror.ge.0) then	    
	      nr=int((a(3)-w1_RA_min)/dangle)
	      nd=int((a(4)-w1_Dec_min)/dangle)
	
	      i=1
	      found=0
	      do while (i.le.ngcat(nr,nd).and.found.eq.0)
	        if (gcat(mngcat(nr,nd,i,1),1,1).eq.a(3)
     ..and.gcat(mngcat(nr,nd,i,1),1,2).eq.a(4)) then
                  found=i
	        else
	          i=i+1
	        endif	      
	      enddo
  	      if (found.eq.0) then
	        ngcat(nr,nd)=ngcat(nr,nd)+1
	        found=ngcat(nr,nd)
	        ngtot=ngtot+1
	        mngcat(nr,nd,found,1)=ngtot
	      endif
	      mngcat(nr,nd,found,2)=mngcat(nr,nd,found,2)+1
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),1)=a(3)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),2)=a(4)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),3)=a(5)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),4)=a(6)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),5)=a(7)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),6)=a(8)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),7)=a(10)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),8)=a(11)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),9)=a(12)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),10)=a(17)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),11)=a(18)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),12)=a(19)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),13)=a(20)
	      gcat(mngcat(nr,nd,found,1)
     .,mngcat(nr,nd,found,2),14)=a(21)
	    endif
	  enddo
	  close(10)
	enddo
        write(*,*) 'Exposure: ',fieldname,exponame
     .,' Merging Shear Catalogs'

c   '#1 gal RA'
c   '#2 gal Dec'
c   '#3 gal SNR'
c   '#4 gal area'
c   '#5 gal flux'
c   '#6 gal Z_B'
c   '#7  lensfit e1'
c   '#8 lensfit e2'
c   '#9 lensfit weight'
c   '#10 g1'
c   '#11 g2'
c   '#12 de'
c   '#13 g1_field'
c   '#14 g2_field'

	
30      format(I2.2)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine map_FD_distribution(fieldname,stat)
        implicit none
        include 'para.inc'

	integer stat
	real gcat(nexpo,npara),a(npara)	
	integer i,j,k,u,v,iexpo,ierror,ii,jj
        character filename*80,fieldname*6
	
	integer ndata,idata
	parameter (ndata=7000000)
	real dfd1(ndata),dfd2(ndata)
	common /num_fdd_pass/ dfd1,dfd2,idata

	if (stat.eq.0) then
	  idata=0
	  return
	elseif (stat.eq.2) then
	  filename='FD1_distribution.fits'
	  call data_binning(ndata,idata,dfd1,0,0.002,0.0002,filename)
	  filename='FD2_distribution.fits'
	  call data_binning(ndata,idata,dfd2,0,0.002,0.0002,filename)
	  return
	endif

	iexpo=0

        filename='./step3/shear_'//fieldname//'.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          pause 'field catalog reading error!!'
        endif
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) i,j,k,u,a(1),a(2),a(3),a(4),a(5)
     .,a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14)
          if (ierror.ge.0) then	    
	    ii=i
	    jj=j
	    if (u.eq.1) then
	      if (iexpo.gt.1) then
	        do i=1,iexpo-1
	          do j=i+1,iexpo
	            idata=idata+1
	            dfd1(idata)=abs(gcat(i,13)-gcat(j,13))
	            dfd2(idata)=abs(gcat(i,14)-gcat(j,14))
	          enddo
	        enddo
	      endif
	      iexpo=0
	    endif
	    iexpo=iexpo+1
            do u=1,npara
	      gcat(iexpo,u)=a(u)
	    enddo
	  endif
	enddo
	close(10)

        if (iexpo.gt.1) then
          do i=1,iexpo-1
            do j=i+1,iexpo
              idata=idata+1
              dfd1(idata)=abs(gcat(i,13)-gcat(j,13))
              dfd2(idata)=abs(gcat(i,14)-gcat(j,14))
	    enddo
          enddo
        endif


        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

