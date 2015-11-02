        subroutine shear_FD_test(fieldname)
        implicit none
        include 'para.inc'
	
        character fieldname*6,filename*80
	integer ierror,i,j,k,u,v,ngal	
	real a(npara),gal_para(ngfieldmax,npara),dg1,dg2
	real g1(nbin,4),g2(nbin,4),yeqx,gf1max,gf1min,gf2max,gf2min
	real de1,de2
	external yeqx

        ngal=0
        filename='./step3/final_shear_'//fieldname//'.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          pause 'Shear catalog reading error!!'
        endif
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) i,j,a(1),a(2),a(3),a(4),a(5)
     .,a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14)
	  if (ierror.ge.0) then
	    ngal=ngal+1
            do v=1,14
	      gal_para(ngal,v)=a(v)
	    enddo
	  endif
	enddo
        close(10)

	gf1max=gal_para(1,13)
	gf1min=gf1max

	gf2max=gal_para(1,14)
	gf2min=gf2max

	do i=2,ngal
	  if (gal_para(i,13).gt.gf1max) gf1max=gal_para(i,13)
	  if (gal_para(i,13).lt.gf1min) gf1min=gal_para(i,13)
	  if (gal_para(i,14).gt.gf2max) gf2max=gal_para(i,14)
	  if (gal_para(i,14).lt.gf2min) gf2min=gal_para(i,14)
	enddo

	gf1max=min(gf1max,0.006)
	gf1min=max(gf1min,-0.006)
	gf2max=min(gf2max,0.006)
	gf2min=max(gf2min,-0.006)
	
	dg1=(gf1max-gf1min)/nbin
	dg2=(gf2max-gf2min)/nbin
	
	do i=1,nbin
	  g1(i,1)=gf1min+dg1*(i-0.5)
	  g2(i,1)=gf2min+dg2*(i-0.5)
	  g1(i,2)=0.
	  g2(i,2)=0.
	  g1(i,3)=0.
	  g2(i,3)=0.
	  g1(i,4)=0.
	  g2(i,4)=0.
	  de1=0.
	  de2=0.
	  do j=1,ngal
	    if (gal_para(j,13).ge.g1(i,1)-dg1*0.5
     ..and.gal_para(j,13).le.g1(i,1)+dg1*0.5) then
	      g1(i,2)=g1(i,2)+gal_para(j,10)
	      de1=de1+gal_para(j,12)
	      g1(i,3)=g1(i,3)+gal_para(j,10)**2
	      g1(i,4)=g1(i,4)+1.
	    endif
	    if (gal_para(j,14).ge.g2(i,1)-dg2*0.5
     ..and.gal_para(j,14).le.g2(i,1)+dg2*0.5) then
	      g2(i,2)=g2(i,2)+gal_para(j,11)
	      de2=de2+gal_para(j,12)
	      g2(i,3)=g2(i,3)+gal_para(j,11)**2
	      g2(i,4)=g2(i,4)+1.
	    endif
	  enddo
	  if (g1(i,4).gt.0) then
            g1(i,2)=g1(i,2)/de1
            g1(i,3)=sqrt(g1(i,3)/de1/de1)
	  endif
	  if (g2(i,4).gt.0) then	
            g2(i,2)=g2(i,2)/de2
            g2(i,3)=sqrt(g2(i,3)/de2/de2)
  	  endif
	  write(*,*) g1(i,1),g1(i,2),g1(i,3),g1(i,4)
	enddo


	filename='./step3/shear_FD_test_g1_'//fieldname//'.fits'
	call map_function(yeqx,nbin,nbin,g1,500,500,filename)

	filename='./step3/shear_FD_test_g2_'//fieldname//'.fits'
	call map_function(yeqx,nbin,nbin,g2,500,500,filename)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shear_FD_test_multiple_field(nf,fieldname)
        implicit none
        include 'para.inc'
	
	integer nf
        character fieldname(nf)*6,filename*80
	integer ierror,i,j,k,u,v,ngal,jf	
	real a(npara),gal_para(ngfieldmax*72,npara),dg1,dg2
	real g1(nbin,4),g2(nbin,4),yeqx,gf1max,gf1min,gf2max,gf2min
	real de1,de2
	external yeqx

        ngal=0
	
	do jf=1,nf
          filename='./step3/final_shear_'//fieldname(jf)//'.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
            pause 'Shear catalog reading error!!'
          endif
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) i,j,a(1),a(2),a(3),a(4),a(5)
     .,a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14)
	    if (ierror.ge.0) then
	      ngal=ngal+1
              do v=1,14
	        gal_para(ngal,v)=a(v)
	      enddo
	    endif
	  enddo
          close(10)
	enddo

	gf1max=gal_para(1,13)
	gf1min=gf1max

	gf2max=gal_para(1,14)
	gf2min=gf2max

	do i=2,ngal
	  if (gal_para(i,13).gt.gf1max) gf1max=gal_para(i,13)
	  if (gal_para(i,13).lt.gf1min) gf1min=gal_para(i,13)
	  if (gal_para(i,14).gt.gf2max) gf2max=gal_para(i,14)
	  if (gal_para(i,14).lt.gf2min) gf2min=gal_para(i,14)
	enddo

	gf1max=min(gf1max,0.005)
	gf1min=max(gf1min,-0.005)
	gf2max=min(gf2max,0.005)
	gf2min=max(gf2min,-0.005)
	
	dg1=(gf1max-gf1min)/nbin
	dg2=(gf2max-gf2min)/nbin
	
	do i=1,nbin
	  g1(i,1)=gf1min+dg1*(i-0.5)
	  g2(i,1)=gf2min+dg2*(i-0.5)
	  g1(i,2)=0.
	  g2(i,2)=0.
	  g1(i,3)=0.
	  g2(i,3)=0.
	  g1(i,4)=0.
	  g2(i,4)=0.
	  de1=0.
	  de2=0.
	  do j=1,ngal
	    if (gal_para(j,13).ge.g1(i,1)-dg1*0.5
     ..and.gal_para(j,13).le.g1(i,1)+dg1*0.5) then
	      g1(i,2)=g1(i,2)+gal_para(j,10)
	      de1=de1+gal_para(j,12)
	      g1(i,3)=g1(i,3)+gal_para(j,10)**2
	      g1(i,4)=g1(i,4)+1.
	    endif
	    if (gal_para(j,14).ge.g2(i,1)-dg2*0.5
     ..and.gal_para(j,14).le.g2(i,1)+dg2*0.5) then
	      g2(i,2)=g2(i,2)+gal_para(j,11)
	      de2=de2+gal_para(j,12)
	      g2(i,3)=g2(i,3)+gal_para(j,11)**2
	      g2(i,4)=g2(i,4)+1.
	    endif
	  enddo
	  if (g1(i,4).gt.0) then
            g1(i,2)=g1(i,2)/de1
            g1(i,3)=sqrt(g1(i,3)/de1/de1)
	  endif
	  if (g2(i,4).gt.0) then	
            g2(i,2)=g2(i,2)/de2
            g2(i,3)=sqrt(g2(i,3)/de2/de2)
  	  endif
	  write(*,*) g1(i,1),g1(i,2),g1(i,3),g1(i,4)
	enddo


	filename='./step3/shear_FD_test_g1_multi.fits'
	call map_function(yeqx,nbin,nbin,g1,500,500,filename)

	filename='./step3/shear_FD_test_g2_multi.fits'
	call map_function(yeqx,nbin,nbin,g2,500,500,filename)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shear_field_test_accu(exponame)
        implicit none
        include 'para.inc'
	
        character exponame(nexpo)*6,s*2,filename*80
	integer iC,ierror,aa,i,j,u,v,nn	
	real a(npara)
	integer ngal_t
	real gal_para_t(ngal_max*nC*10,npara),dg1,dg2

	integer ngal
	real gal_para(ngal_max*nC,npara)
	common /gal_para_pass/ gal_para,ngal	

	real g1(nbin,4),g2(nbin,4),yeqx,gf1max,gf1min,gf2max,gf2min
	real de1,de2,w,func_weight
	external yeqx,func_weight

        ngal_t=0
	
	do i=1,nexpo	
	  call shear_field_test(exponame(i))
	  do j=1,ngal
            ngal_t=ngal_t+1
            do u=1,15
              gal_para_t(ngal_t,u)=gal_para(j,u)
            enddo
	  enddo
	  write(*,*) 'Exposure:',exponame,'Collecting Shears'
	enddo

	gf1max=gal_para_t(1,14)
	gf1min=gf1max

	gf2max=gal_para_t(1,15)
	gf2min=gf2max

	do i=2,ngal_t
	  if (gal_para_t(i,14).gt.gf1max) gf1max=gal_para_t(i,14)
	  if (gal_para_t(i,14).lt.gf1min) gf1min=gal_para_t(i,14)
	  if (gal_para_t(i,15).gt.gf2max) gf2max=gal_para_t(i,15)
	  if (gal_para_t(i,15).lt.gf2min) gf2min=gal_para_t(i,15)
	enddo

	dg1=(gf1max-gf1min)/nbin
	dg2=(gf2max-gf2min)/nbin
	
	do i=1,nbin
	  g1(i,1)=gf1min+dg1*(i-0.5)
	  g2(i,1)=gf2min+dg2*(i-0.5)
	  g1(i,2)=0.
	  g2(i,2)=0.
	  g1(i,3)=0.
	  g2(i,3)=0.
	  g1(i,4)=0.
	  g2(i,4)=0.
	  de1=0.
	  de2=0.
	  do j=1,ngal_t
	    if (gal_para_t(j,14).ge.g1(i,1)-dg1*0.5
     ..and.gal_para_t(j,14).le.g1(i,1)+dg1*0.5) then
	      w=func_weight(gal_para_t(j,5))
	      g1(i,2)=g1(i,2)+gal_para_t(j,11)*w
	      de1=de1+gal_para_t(j,13)*w
	      g1(i,3)=g1(i,3)+(gal_para_t(j,11)*w)**2
	      g1(i,4)=g1(i,4)+1.
	    endif
	    if (gal_para_t(j,15).ge.g2(i,1)-dg2*0.5
     ..and.gal_para_t(j,15).le.g2(i,1)+dg2*0.5) then
	      w=func_weight(gal_para_t(j,5))
	      g2(i,2)=g2(i,2)+gal_para_t(j,12)*w
	      de2=de2+gal_para_t(j,13)*w
	      g2(i,3)=g2(i,3)+(gal_para_t(j,12)*w)**2
	      g2(i,4)=g2(i,4)+1.
	    endif
	  enddo
	  if (g1(i,4).gt.0) then
            g1(i,2)=g1(i,2)/de1
            g1(i,3)=sqrt(g1(i,3)/de1/de1)
	  endif
	  if (g2(i,4).gt.0) then	
            g2(i,2)=g2(i,2)/de2
            g2(i,3)=sqrt(g2(i,3)/de2/de2)
  	  endif
	enddo


	filename='./step3/shear_field_test_g1_all.fits'
	call map_function(yeqx,nbin,nbin,g1,500,500,filename)

	filename='./step3/shear_field_test_g2_all.fits'
	call map_function(yeqx,nbin,nbin,g2,500,500,filename)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shear_FD_test_multiple_field_CFHTlens(nf,fieldname)
        implicit none
        include 'para.inc'
	
	integer nf
        character fieldname(nf)*6,filename*80
	integer ierror,i,j,k,u,v,ngal,jf	
	real a(npara),gal_para(ngfieldmax*72,npara),dg1,dg2
	real g1(nbin,4),g2(nbin,4),yeqx,gf1max,gf1min,gf2max,gf2min
	real de1,de2
	external yeqx

        ngal=0
	
	do jf=1,nf
          filename='./step3/final_shear_'//fieldname(jf)//'.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
            pause 'Shear catalog reading error!!'
          endif
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) i,j,a(1),a(2),a(3),a(4),a(5)
     .,a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14)
	    if (ierror.ge.0) then
	      ngal=ngal+1
              do v=1,14
	        gal_para(ngal,v)=a(v)
	      enddo
	    endif
	  enddo
          close(10)
	enddo

	gf1max=gal_para(1,13)
	gf1min=gf1max

	gf2max=gal_para(1,14)
	gf2min=gf2max

	do i=2,ngal
	  if (gal_para(i,13).gt.gf1max) gf1max=gal_para(i,13)
	  if (gal_para(i,13).lt.gf1min) gf1min=gal_para(i,13)
	  if (gal_para(i,14).gt.gf2max) gf2max=gal_para(i,14)
	  if (gal_para(i,14).lt.gf2min) gf2min=gal_para(i,14)
	enddo

	gf1max=min(gf1max,0.005)
	gf1min=max(gf1min,-0.005)
	gf2max=min(gf2max,0.005)
	gf2min=max(gf2min,-0.005)
	
	dg1=(gf1max-gf1min)/nbin
	dg2=(gf2max-gf2min)/nbin
	
	do i=1,nbin
	  g1(i,1)=gf1min+dg1*(i-0.5)
	  g2(i,1)=gf2min+dg2*(i-0.5)
	  g1(i,2)=0.
	  g2(i,2)=0.
	  g1(i,3)=0.
	  g2(i,3)=0.
	  g1(i,4)=0.
	  g2(i,4)=0.
	  de1=0.
	  de2=0.
	  do j=1,ngal
	    if (gal_para(j,13).ge.g1(i,1)-dg1*0.5
     ..and.gal_para(j,13).le.g1(i,1)+dg1*0.5) then
	      g1(i,2)=g1(i,2)+gal_para(j,7)*gal_para(j,9)
	      de1=de1+gal_para(j,9)
	      g1(i,3)=g1(i,3)+(gal_para(j,7)*gal_para(j,9))**2
	      g1(i,4)=g1(i,4)+1.
	    endif
	    if (gal_para(j,14).ge.g2(i,1)-dg2*0.5
     ..and.gal_para(j,14).le.g2(i,1)+dg2*0.5) then
	      g2(i,2)=g2(i,2)+gal_para(j,8)*gal_para(j,9)
	      de2=de2+gal_para(j,9)
	      g2(i,3)=g2(i,3)+(gal_para(j,8)*gal_para(j,9))**2
	      g2(i,4)=g2(i,4)+1.
	    endif
	  enddo
	  if (g1(i,4).gt.0) then
            g1(i,2)=g1(i,2)/de1
            g1(i,3)=sqrt(g1(i,3)/de1/de1)
	  endif
	  if (g2(i,4).gt.0) then	
            g2(i,2)=g2(i,2)/de2
            g2(i,3)=sqrt(g2(i,3)/de2/de2)
  	  endif
	  write(*,*) g1(i,1),g1(i,2),g1(i,3),g1(i,4)
	enddo


	filename='./step3/shear_FD_test_g1_multi_CFHTlens.fits'
	call map_function(yeqx,nbin,nbin,g1,500,500,filename)

	filename='./step3/shear_FD_test_g2_multi_CFHTlens.fits'
	call map_function(yeqx,nbin,nbin,g2,500,500,filename)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shear_field_test(exponame)
        implicit none
        include 'para.inc'
	
        character exponame*6,s*2,filename*80
	integer iC,ierror,aa,i,j,u,v,nn	
	real a(npara)
	integer ngal
	real gal_para(ngal_max*nC,npara),dg1,dg2
	common /gal_para_pass/ gal_para,ngal	

	real g1(nbin,4),g2(nbin,4),yeqx,gf1max,gf1min,gf2max,gf2min
	real de1,de2,w,func_weight
	external yeqx,func_weight

        ngal=0
	do iC=1,36
          write(s,30) iC
          filename='./step2/shear_info'//exponame//'_'//s//'.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
	    pause 'Gal catalog reading error!!'
          endif
	  do i=1,16
	    read(10,*)
	  enddo
          do while (ierror.ge.0)
  	    read(10,*,iostat=ierror) aa,a(1),a(2),a(3),a(4),a(5)
     .,a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15)
	    if (ierror.ge.0) then
	      ngal=ngal+1
	      do i=1,15
	        gal_para(ngal,i)=a(i)
	      enddo
	    endif
	  enddo
	  close(10)
	enddo
	
30      format(I2.2)

	gf1max=gal_para(1,14)
	gf1min=gf1max

	gf2max=gal_para(1,15)
	gf2min=gf2max

	do i=2,ngal
	  if (gal_para(i,14).gt.gf1max) gf1max=gal_para(i,14)
	  if (gal_para(i,14).lt.gf1min) gf1min=gal_para(i,14)
	  if (gal_para(i,15).gt.gf2max) gf2max=gal_para(i,15)
	  if (gal_para(i,15).lt.gf2min) gf2min=gal_para(i,15)
	enddo

	dg1=(gf1max-gf1min)/nbin
	dg2=(gf2max-gf2min)/nbin
	
	do i=1,nbin
	  g1(i,1)=gf1min+dg1*(i-0.5)
	  g2(i,1)=gf2min+dg2*(i-0.5)
	  g1(i,2)=0.
	  g2(i,2)=0.
	  g1(i,3)=0.
	  g2(i,3)=0.
	  g1(i,4)=0.
	  g2(i,4)=0.
	  de1=0.
	  de2=0.
	  do j=1,ngal
	    if (gal_para(j,14).ge.g1(i,1)-dg1*0.5
     ..and.gal_para(j,14).le.g1(i,1)+dg1*0.5) then
	      w=func_weight(gal_para(j,5))
	      g1(i,2)=g1(i,2)+gal_para(j,11)*w
	      de1=de1+gal_para(j,13)*w
	      g1(i,3)=g1(i,3)+(gal_para(j,11)*w)**2
	      g1(i,4)=g1(i,4)+1.
	    endif
	    if (gal_para(j,15).ge.g2(i,1)-dg2*0.5
     ..and.gal_para(j,15).le.g2(i,1)+dg2*0.5) then
	      w=func_weight(gal_para(j,5))
	      g2(i,2)=g2(i,2)+gal_para(j,12)*w
	      de2=de2+gal_para(j,13)*w
	      g2(i,3)=g2(i,3)+(gal_para(j,12)*w)**2
	      g2(i,4)=g2(i,4)+1.
	    endif
	  enddo
	  if (g1(i,4).gt.0) then
            g1(i,2)=g1(i,2)/de1
            g1(i,3)=sqrt(g1(i,3)/de1/de1)
	  endif
	  if (g2(i,4).gt.0) then	
            g2(i,2)=g2(i,2)/de2
            g2(i,3)=sqrt(g2(i,3)/de2/de2)
  	  endif
	enddo


	filename='./step3/shear_field_test_g1_'//exponame//'.fits'
	call map_function(yeqx,nbin,nbin,g1,500,500,filename)

	filename='./step3/shear_field_test_g2_'//exponame//'.fits'
	call map_function(yeqx,nbin,nbin,g2,500,500,filename)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function yeqx(x)
	implicit none
	real yeqx,x	
	yeqx=x
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

