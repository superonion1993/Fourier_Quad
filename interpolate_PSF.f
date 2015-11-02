        subroutine prepare_star_kriging(exponame)
        implicit none
	include 'para.inc'
	
	character exponame*6
	integer iC
	real corr_ss(npcorr,2)
	integer n_ss
	common /corr_pass/ corr_ss,n_ss
	character imagename*80

	real a,b
	real dist(nC,nstar_max,nstar_max),cov(nC,nstar_max,nstar_max)
	common /kriging_pass/ dist,cov,a,b
	
	integer u,v
	
	n_ss=0
	do iC=1,36
	  call accumulate_corr(iC)
	enddo
	
	imagename='./step2/variogram'//exponame//'.fits'	

	call make_scatter_map(npcorr,n_ss,corr_ss,10,500,500
     .,imagename,a,b)

	a=exp(a)*0.5

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine accumulate_corr(iC)
	implicit none
	include 'para.inc'

	integer iC
	integer nPSF(nC)
	real PSF(nC,nstar_max,ns,ns),PSF_para(nC,nstar_max,npara)
	common /PSF_pass/ PSF,PSF_para,nPSF

	real corr_ss(npcorr,2)
	integer n_ss
	common /corr_pass/ corr_ss,n_ss
	
	integer ncorr,i,j,u,v,nstar

	real a,b
	real dist(nC,nstar_max,nstar_max),cov(nC,nstar_max,nstar_max)
	common /kriging_pass/ dist,cov,a,b

	nstar=nPSF(iC)

	do i=1,nstar-1
	  do j=i+1,nstar
	    n_ss=n_ss+1

	    dist(iC,i,j)
     .=sqrt((PSF_para(iC,i,1)-PSF_para(iC,j,1))**2
     .+(PSF_para(iC,i,2)-PSF_para(iC,j,2))**2)
	    dist(iC,j,i)=dist(iC,i,j)

	    corr_ss(n_ss,1)=dist(iC,i,j)

	    corr_ss(n_ss,2)=0.
	    do u=1,ns
	      do v=1,ns
	        corr_ss(n_ss,2)=corr_ss(n_ss,2)
     .+(PSF(iC,i,u,v)-PSF(iC,j,u,v))**2
	      enddo
	    enddo
	    corr_ss(n_ss,1)=log(corr_ss(n_ss,1))
	    corr_ss(n_ss,2)=log(corr_ss(n_ss,2))
	  enddo
	enddo

	do i=1,nstar
	  dist(iC,i,i)=0.
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function variogram(r,a,b)
	implicit none

	real a,b,r,variogram

	variogram=a*r**b

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine kriging_setup(iC)
	implicit none
	include 'para.inc'

	integer iC
	integer nPSF(nC)
	real PSF(nC,nstar_max,ns,ns),PSF_para(nC,nstar_max,npara)
	common /PSF_pass/ PSF,PSF_para,nPSF

	real a,b
	real dist(nC,nstar_max,nstar_max),cov(nC,nstar_max,nstar_max)
	common /kriging_pass/ dist,cov,a,b

	integer nstar,u,v
	real variogram

	nstar=nPSF(iC)
        do u=1,nstar
          do v=1,nstar
            cov(iC,u,v)=variogram(dist(iC,u,v),a,b)
          enddo
          cov(iC,u,nstar+1)=1.
          cov(iC,nstar+1,u)=1.
        enddo
        cov(iC,nstar+1,nstar+1)=0.

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine interpolate_kriging_PSF(iC,x,y,model)
	implicit none
	include 'para.inc'

	integer iC
	real model(ns,ns),x,y,variogram
	real r,vec(nstar_max),covv(nstar_max,nstar_max)
	
	integer nPSF(nC)
	real PSF(nC,nstar_max,ns,ns),PSF_para(nC,nstar_max,npara)
	common /PSF_pass/ PSF,PSF_para,nPSF

	integer i,j,u,v,nstar

	real a,b
	real dist(nC,nstar_max,nstar_max),cov(nC,nstar_max,nstar_max)
	common /kriging_pass/ dist,cov,a,b

	nstar=nPSF(iC)

        do u=1,nstar
	  r=sqrt((PSF_para(iC,u,1)-x)**2+(PSF_para(iC,u,2)-y)**2)
          vec(u)=variogram(r,a,b)
        enddo
	vec(nstar+1)=1

	do u=1,nstar+1
	  do v=1,nstar+1
	    covv(u,v)=cov(iC,u,v)
	  enddo
	enddo

	call gaussj(covv,nstar+1,nstar_max,vec,1,1)

	do u=1,ns
	  do v=1,ns
	    model(u,v)=0.
	    do i=1,nstar
	      model(u,v)=model(u,v)+PSF(iC,i,u,v)*vec(i)
	    enddo
	  enddo
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interpolate_radial_PSF(iC,x,y,model)
	implicit none
	include 'para.inc'

	integer iC
	real model(ns,ns),x,y,temp,r2,weight
	integer nPSF(nC)
	real PSF(nC,nstar_max,ns,ns),PSF_para(nC,nstar_max,npara)
	common /PSF_pass/ PSF,PSF_para,nPSF
	integer i,j,u,v,nstar

	nstar=nPSF(iC)

	do u=1,ns
	  do v=1,ns
	    model(u,v)=0.
	  enddo
	enddo
	
	temp=0.
	do i=1,nstar
	  r2=(PSF_para(iC,i,1)-x)**2+(PSF_para(iC,i,2)-y)**2
          if (r2.le.rc*rc) then
	    weight=r2**alpha
            temp=temp+weight
            do u=1,ns
	      do v=1,ns
	        model(u,v)=model(u,v)+PSF(iC,i,u,v)*weight
	      enddo
	    enddo
          endif
	enddo

	temp=1./temp

        do u=1,ns
          do v=1,ns
            model(u,v)=model(u,v)*temp
          enddo
        enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine poly_PSF_setup(iC)
	implicit none
	include 'para.inc'

	integer iC
	integer nPSF(nC)
	real PSF(nC,nstar_max,ns,ns),PSF_para(nC,nstar_max,npara)
	common /PSF_pass/ PSF,PSF_para,nPSF

	integer nstar,u,v,i,j
	real arr(nstar_max,4),c3(10),c2(6)

	real poly_coe(nC,10,ns,ns)
	common /poly_pass/ poly_coe

	nstar=nPSF(iC)

        do i=1,nstar
          arr(i,1)=PSF_para(iC,i,1)
          arr(i,2)=PSF_para(iC,i,2)
          arr(i,4)=1.	
        enddo

	do u=1,ns
	  do v=1,ns
	    do i=1,nstar
	      arr(i,3)=PSF(iC,i,u,v)
	    enddo
	    if (norder.eq.2) then
	      call fit_poly_2D_2(nstar_max,nstar,arr,c2)
	      do i=1,6
	        poly_coe(iC,i,u,v)=c2(i)
	      enddo
	    elseif (norder.eq.3) then
	      call fit_poly_2D_3(nstar_max,nstar,arr,c3)
	      do i=1,10
	        poly_coe(iC,i,u,v)=c3(i)
	      enddo
	    endif	   
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interpolate_poly_PSF(iC,x,y,model)
	implicit none
	include 'para.inc'

	integer iC
	real x,y,model(ns,ns),c2(6),c3(10)
	real poly_coe(nC,10,ns,ns)
	common /poly_pass/ poly_coe

	integer i,j,u,v

	do u=1,ns
	  do v=1,ns
	    if (norder.eq.2) then
	      do i=1,6
	        c2(i)=poly_coe(iC,i,u,v)
	      enddo
	      model(u,v)=c2(1)+c2(2)*x+c2(3)*y
     .+c2(4)*x*x+c2(5)*x*y+c2(6)*y*y

	    elseif (norder.eq.3) then
	      do i=1,10
	        c3(i)=poly_coe(iC,i,u,v)
	      enddo
	      model(u,v)=c3(1)+c3(2)*x+c3(3)*y
     .+c3(4)*x*x+c3(5)*x*y+c3(6)*y*y+c3(7)*x**3+c3(8)*x*x*y
     .+c3(9)*x*y*y+c3(10)*y**3
	    endif
	  enddo
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mean_PSF_setup(iC)
        implicit none
	include 'para.inc'

	integer iC
	integer nPSF(nC)
	real PSF(nC,nstar_max,ns,ns),PSF_para(nC,nstar_max,npara)
	common /PSF_pass/ PSF,PSF_para,nPSF

	real mean_PSF(nC,ns,ns)
	common /mean_PSF_pass/ mean_PSF

	real temp
	integer i,j,u,v,nstar

	nstar=nPSF(iC)

        do u=1,ns
          do v=1,ns
            mean_PSF(iC,u,v)=0.
          enddo
        enddo

	do i=1,nstar	  
          do u=1,ns
            do v=1,ns
              mean_PSF(iC,u,v)=mean_PSF(iC,u,v)+PSF(iC,i,u,v)
            enddo
          enddo
	enddo

	temp=1./nstar

	do u=1,ns
	  do v=1,ns
            mean_PSF(iC,u,v)=mean_PSF(iC,u,v)*temp
	  enddo
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine average_PSF(iC,x,y,model)
        implicit none
	include 'para.inc'

	integer iC
	real x,y,model(ns,ns)

	real mean_PSF(nC,ns,ns)
	common /mean_PSF_pass/ mean_PSF

	integer u,v

        do u=1,ns
          do v=1,ns
            model(u,v)=mean_PSF(iC,u,v)
          enddo
        enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc		
