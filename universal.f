	subroutine normalize_gal(nx,ny,gal,noise)
	implicit none

	integer nx,ny
	real gal(nx,ny),noise(nx,ny),temp
	integer i,j,cx,cy

	cx=nx/2+1
	cy=ny/2+1

	temp=1./gal(cx,cy)
	do i=1,nx
	  do j=1,ny
	    gal(i,j)=gal(i,j)*temp
	    noise(i,j)=noise(i,j)*temp
	  enddo
	enddo
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine trim_power(nx,ny,power)
	implicit none

	integer nx,ny
	real power(nx,ny)
	integer i,j
	real temp

	temp=0.
	do i=1,nx
	  temp=temp+power(i,1)+power(i,ny)
	enddo	
	do i=2,ny-1
	  temp=temp+power(1,i)+power(nx,i)
	enddo	
			
	temp=temp/(2.*(nx+ny)-4.)

	do i=1,nx
	  do j=1,ny
	    power(i,j)=power(i,j)-temp
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine regularize_power(nx,ny,power)
	implicit none

	integer nx,ny
	real power(nx,ny),temp
	integer i,j,cx,cy

	cx=nx/2+1
	cy=ny/2+1

	call trim_power(nx,ny,power)

	temp=1./power(cx,cy)
	do i=1,nx
	  do j=1,ny
	    power(i,j)=power(i,j)*temp
	  enddo
	enddo
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine matrix_inverse2(ma,n,np,vec)
	implicit none

	integer n,np
	real ma(np,np),vec(np),ma_1(n,n),d,b(n),tmp(n)
	integer indx(n),i,j

        call ludcmp(ma,n,np,indx,d)

	do i=1,n
	  do j=1,n
	    if (j.eq.i) then
	      b(j)=1.
	    else 
	      b(j)=0.
	    endif
	  enddo	
          call lubksb(ma,n,np,indx,b)
	  do j=1,n
            ma_1(j,i)=b(j)
	  enddo	
	enddo
		
	do i=1,n
	  tmp(i)=vec(i)
	enddo

	do i=1,n
	  vec(i)=0
	  do j=1,n
	    vec(i)=vec(i)+ma_1(i,j)*tmp(j)
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine matrix_inverse(ma,n,ma_1)
	implicit none

	integer n
	real ma(n,n),ma_1(n,n),d,b(n)
	integer indx(n),i,j

        call ludcmp(ma,n,n,indx,d)

	do i=1,n
	  do j=1,n
	    if (j.eq.i) then
	      b(j)=1.
	    else 
	      b(j)=0.
	    endif
	  enddo	

          call lubksb(ma,n,n,indx,b)

	  do j=1,n
            ma_1(j,i)=b(j)
	  enddo	
	enddo
		
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_poly_2D_3(np,n,arr,c)
	implicit none

c fit a 2D function with f=c1+c2x+c3y+c4x^2+c5xy+c6y^2+c7x^3+c8x^2y+c9xy^2+c10y^3
	integer np,n
	real arr(np,4),c(10),s_2
	real coe(10,10),coe_1(10,10),vec(10),temp(10)
	real x,y,f
	integer i,j,u,v

	do i=1,10
	  do j=1,10
	    coe(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	do i=1,n
	  x=arr(i,1)
	  y=arr(i,2)
	  f=arr(i,3)
	  s_2=arr(i,4)

	  temp(1)=1.
	  temp(2)=x
	  temp(3)=y
	  temp(4)=x*x
	  temp(5)=x*y
	  temp(6)=y*y
	  temp(7)=x*x*x
	  temp(8)=x*x*y
	  temp(9)=x*y*y
	  temp(10)=y*y*y

	  do u=1,10
	    do v=1,10
	      coe(u,v)=coe(u,v)+temp(u)*temp(v)*s_2
	    enddo
	    vec(u)=vec(u)+f*temp(u)*s_2		  
	  enddo
    
	enddo

        call matrix_inverse(coe,10,coe_1)

	do j=1,10
	  c(j)=0.
	  do i=1,10
	    c(j)=c(j)+coe_1(j,i)*vec(i)
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_poly_2D_2(np,n,arr,c)
	implicit none

c fit a 2D function with f=c1+c2x+c3y+c4x^2+c5xy+c6y^2
	integer np,n
	real arr(np,4),c(6),s_2
	real coe(6,6),coe_1(6,6),vec(6),temp(6)
	real x,y,f
	integer i,j,u,v

	do i=1,6
	  do j=1,6
	    coe(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	do i=1,n
	  x=arr(i,1)
	  y=arr(i,2)
	  f=arr(i,3)
	  s_2=arr(i,4)

	  temp(1)=1.
	  temp(2)=x
	  temp(3)=y
	  temp(4)=x*x
	  temp(5)=x*y
	  temp(6)=y*y

	  do u=1,6
	    do v=1,6
	      coe(u,v)=coe(u,v)+temp(u)*temp(v)*s_2
	    enddo
	    vec(u)=vec(u)+f*temp(u)*s_2		  
	  enddo
    
	enddo

        call matrix_inverse(coe,6,coe_1)

	do j=1,6
	  c(j)=0.
	  do i=1,6
	    c(j)=c(j)+coe_1(j,i)*vec(i)
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find_slope_2D(np,n,arr,aa,bb,cc)
	implicit none

	integer np,n
	real arr(np,3),aa,bb,cc
	integer i,j
	real c(3,3),vec(3),c_1(3,3)

	do i=1,3
	  do j=1,3
	    c(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	do i=1,n
	  c(1,1)=c(1,1)+1.
	  c(1,2)=c(1,2)+arr(i,1)
	  c(1,3)=c(1,3)+arr(i,2)
	  vec(1)=vec(1)+arr(i,3)		  

	  c(2,1)=c(2,1)+arr(i,1)
	  c(2,2)=c(2,2)+arr(i,1)*arr(i,1)
	  c(2,3)=c(2,3)+arr(i,2)*arr(i,1)
	  vec(2)=vec(2)+arr(i,3)*arr(i,1)		  

	  c(3,1)=c(3,1)+arr(i,2)
	  c(3,2)=c(3,2)+arr(i,1)*arr(i,2)
	  c(3,3)=c(3,3)+arr(i,2)*arr(i,2)
	  vec(3)=vec(3)+arr(i,3)*arr(i,2)		  
	enddo

        call matrix_inverse(c,3,c_1)

	aa=0.
	bb=0.
	cc=0.
	do i=1,3
	  aa=aa+c_1(1,i)*vec(i)
	  bb=bb+c_1(2,i)*vec(i)
	  cc=cc+c_1(3,i)*vec(i)
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Gauss_map(sig_2,nx,ny,map)
      implicit none
      
      real sig_2
      integer nx,ny
      real map(nx,ny)

      integer i,j,cx,cy
      real x,y,x2,y2
		
      cx=nx/2+1
      cy=ny/2+1
	
      do i=1,nx
	x=i-cx
	x2=x*x
	do j=1,ny
	  y=j-cy
	  y2=y*y
	  map(i,j)=exp(-(x2+y2)*sig_2*0.5)
	enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reconvolution(sig_2,nx,ny,power,map)
      implicit none
      
      real sig_2
      integer nx,ny
      real power(nx,ny),map(nx,ny)

      integer i,j
		
      call Gauss_map(sig_2,nx,ny,map)

      do i=1,nx
	do j=1,ny
	  map(i,j)=map(i,j)/power(i,j)
	enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_power(nx,ny,power,nbp,n,npx,npy,basis
     .,np,coe,model,diff)
	implicit none

	integer np,n,nx,ny,nbp,npx,npy
	real power(nx,ny),basis(nbp,npx,npy),coe(np)
	real model(nx,ny),diff(nx,ny)

	real cov(n,n),cov_1(n,n),vec(n)
	real temp(nx,ny)
	integer i,j,u,v

	integer cx,cy
	
	cx=nx/2+1
	cy=ny/2+1

	do u=1,nx
	  do v=1,ny
	    temp(u,v)=power(u,v)**(-1)
	  enddo
	enddo

	do i=1,n
	  do j=1,n	    
	    cov(i,j)=0.
	    do u=1,nx
	      do v=1,ny
	        if (u.ne.cx.or.v.ne.cy) then	    
	          cov(i,j)=cov(i,j)
     .+basis(i,u,v)*basis(j,u,v)*temp(u,v)
	        endif
	      enddo
	    enddo
	  enddo
	enddo

	do i=1,n
          vec(i)=0.
          do u=1,nx
            do v=1,ny
              if (u.ne.cx.or.v.ne.cy) then	    
                vec(i)=vec(i)+power(u,v)*basis(i,u,v)*temp(u,v)
	      endif
            enddo
          enddo
	enddo
		
	call matrix_inverse(cov,n,cov_1)
	
	do i=1,n
	  coe(i)=0
	  do j=1,n
	    coe(i)=coe(i)+cov_1(i,j)*vec(j)
	  enddo
	enddo
	
	do u=1,nx
	  do v=1,ny
	    model(u,v)=0.
	    do i=1,n
	      model(u,v)=model(u,v)+basis(i,u,v)*coe(i)	
	    enddo
	    diff(u,v)=power(u,v)-model(u,v)
	  enddo
	enddo
	  	

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find_PCA_weighted(np,n,nx,ny,map,weight,PCA)
	implicit none

	integer np,n,nx,ny
	real map(np,nx,ny),weight(np),PCA(nx,ny)
	integer i,j,k,u,v,ip
	real mean(nx,ny),sig(nx,ny)

	integer NMAX,p,q
	parameter (NMAX=1000)
	real matx(NMAX,NMAX),val(NMAX),vec(NMAX,NMAX)
	real temp
	
	do i=1,nx	
	  do j=1,ny
	    mean(i,j)=0.
	    temp=0.
	    do k=1,n
	      mean(i,j)=mean(i,j)+map(k,i,j)*weight(k)
	      temp=temp+weight(k)
	    enddo
            mean(i,j)=mean(i,j)/temp
	    sig(i,j)=0.
	    do k=1,n
	      map(k,i,j)=map(k,i,j)-mean(i,j)
	      sig(i,j)=sig(i,j)+map(k,i,j)**2
	    enddo
            sig(i,j)=sqrt(sig(i,j)/n)	
	    do k=1,n
	      map(k,i,j)=map(k,i,j)/sig(i,j)
	    enddo
	  enddo
	enddo

	p=0
	do i=1,nx
	  do j=1,ny
	    p=p+1
	    q=0
	    do u=1,nx
	      do v=1,ny
	        q=q+1
	        matx(p,q)=0.
	        do k=1,n
	          matx(p,q)=matx(p,q)+map(k,i,j)*map(k,u,v)*weight(k)
	        enddo
	      enddo
	    enddo
	  enddo
	enddo

	call find_eigen_vec(NMAX,nx*ny,matx,val,vec)
	
	ip=1
	temp=val(ip)
	do j=2,nx*ny
	  if (val(j).gt.temp) then
	    temp=val(j)
	    ip=j
	  endif
	enddo

	p=0
	do i=1,nx
	  do j=1,ny
	    p=p+1
	    PCA(i,j)=vec(p,ip)
	  enddo
	enddo  

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find_PCA(np,n,nx,ny,map,PCA)
	implicit none

	integer np,n,nx,ny
	real map(np,nx,ny),PCA(nx,ny)
	integer i,j,k,u,v,ip
	real mean(nx,ny),sig(nx,ny)

	integer NMAX,p,q
	parameter (NMAX=1000)
	real matx(NMAX,NMAX),val(NMAX),vec(NMAX,NMAX)
	real temp
	
	do i=1,nx	
	  do j=1,ny
	    mean(i,j)=0.
	    do k=1,n
	      mean(i,j)=mean(i,j)+map(k,i,j)
	    enddo
            mean(i,j)=mean(i,j)/n
	    sig(i,j)=0.
	    do k=1,n
	      map(k,i,j)=map(k,i,j)-mean(i,j)
	      sig(i,j)=sig(i,j)+map(k,i,j)**2
	    enddo
            sig(i,j)=sqrt(sig(i,j)/n)	
	    do k=1,n
	      map(k,i,j)=map(k,i,j)/sig(i,j)
	    enddo
	  enddo
	enddo

	p=0
	do i=1,nx
	  do j=1,ny
	    p=p+1
	    q=0
	    do u=1,nx
	      do v=1,ny
	        q=q+1
	        matx(p,q)=0.
	        do k=1,n
	          matx(p,q)=matx(p,q)+map(k,i,j)*map(k,u,v)
	        enddo
	      enddo
	    enddo
	  enddo
	enddo

	call find_eigen_vec(NMAX,nx*ny,matx,val,vec)
	
	ip=1
	temp=val(ip)
	do j=2,nx*ny
	  if (val(j).gt.temp) then
	    temp=val(j)
	    ip=j
	  endif
	enddo

	p=0
	do i=1,nx
	  do j=1,ny
	    p=p+1
	    PCA(i,j)=vec(p,ip)
	  enddo
	enddo  

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find_eigen_vec(np,n,matx,val,vec)
	implicit none

	integer np,n
	real matx(np,np),val(np),vec(np,np)
	real d(np),e(np)	
	integer i,j,k
	
	do i=1,n
	  do j=1,n
	    vec(i,j)=matx(i,j)
	  enddo
	enddo

	
 	call tred2(vec,n,np,d,e)

	call tqli(d,e,n,np,vec)

	do i=1,n
          val(i)=d(i)
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gen_orthonormal_basis(np,n,nx,ny,stamp,nmax,basis)
	implicit none

	integer np,n,nx,ny,nmax
	real stamp(np,nx,ny),basis(np,nx,ny),coe(np)
	real sample(nx,ny),model(nx,ny),diff(nx,ny)
	real temp,weight(np),norm,temp1
	integer i,j,u,v,k,cx,cy,ipca,m

	character filename*80

	cx=nx/2+1
	cy=ny/2+1

	do i=1,n
	  weight(i)=stamp(i,cx,cy)**0.5
	  temp=1./stamp(i,cx,cy)
	  do u=1,nx
	    do v=1,ny
	      stamp(i,u,v)=stamp(i,u,v)*temp
	    enddo
	  enddo
	enddo
	
	temp=0.
	do i=1,n
	  temp=temp+weight(i)
	enddo
	temp=1./temp
	do i=1,n
	  weight(i)=weight(i)*temp
	enddo

        do u=1,nx
          do v=1,ny
            model(u,v)=0.
  	    do i=1,n 	
              model(u,v)=model(u,v)+weight(i)*stamp(i,u,v)
            enddo
          enddo
        enddo
	call normalize_image(nx,ny,model)
        do u=1,nx
          do v=1,ny
            basis(1,u,v)=model(u,v)
          enddo
        enddo

	do k=2,nmax
	  do i=1,n 	
	    temp=0.
            do u=1,nx
              do v=1,ny
                temp=temp+basis(k-1,u,v)*stamp(i,u,v)
              enddo
            enddo
            do u=1,nx
              do v=1,ny
                stamp(i,u,v)=stamp(i,u,v)-temp*basis(k-1,u,v)
              enddo
            enddo
	    ipca=k-1
	  enddo

          do u=1,nx
            do v=1,ny
	      model(u,v)=stamp(ipca,u,v)
            enddo
          enddo
		
c	  do m=1,5  
c	    call smooth_image55(nx,ny,model)
c	  enddo

c	  do i=1,k-1 	
c	    temp=0.
c           do u=1,nx
c              do v=1,ny
c                temp=temp+basis(i,u,v)*model(u,v)
c              enddo
c            enddo
c            do u=1,nx
c              do v=1,ny
c                model(u,v)=model(u,v)-temp*basis(i,u,v)
c              enddo
c            enddo
c	  enddo
	  call normalize_image(nx,ny,model)
          do u=1,nx
            do v=1,ny
	      basis(k,u,v)=model(u,v)
            enddo
          enddo
	 
c        filename='star_basis.fits'		
c	call show_stamps(np,k,nx,ny,basis,nx*6,ny*6,filename)
c       filename='star_residual.fits'		
c	call show_stamps(np,n,nx,ny,stamp,nx*6,ny*6,filename)
c	pause

	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gen_orthonormal_basis_PCA(np,n,nx,ny,stamp
     .,nmax,basis)
	implicit none

	integer np,n,nx,ny,nmax
	real stamp(np,nx,ny),basis(np,nx,ny),coe(np)
	real PCAstamp(np,nx,ny),PCAweight(np)

	real sample(nx,ny),model(nx,ny),diff(nx,ny)
	real temp,weight(np),norm,temp1
	integer i,j,u,v,k,cx,cy,ipca,m

	character filename*80

	cx=nx/2+1
	cy=ny/2+1

	do i=1,n
	  weight(i)=stamp(i,cx,cy)**0.5
	  temp=1./stamp(i,cx,cy)
	  do u=1,nx
	    do v=1,ny
	      stamp(i,u,v)=stamp(i,u,v)*temp
	    enddo
	  enddo
	enddo
	
	temp=0.
	do i=1,n
	  temp=temp+weight(i)
	enddo
	temp=1./temp
	do i=1,n
	  weight(i)=weight(i)*temp
	enddo

        do u=1,nx
          do v=1,ny
            model(u,v)=0.
  	    do i=1,n 	
              model(u,v)=model(u,v)+weight(i)*stamp(i,u,v)
            enddo
          enddo
        enddo
	call normalize_image(nx,ny,model)
        do u=1,nx
          do v=1,ny
            basis(1,u,v)=model(u,v)
          enddo
        enddo

	do k=2,nmax
	  do i=1,n 	
	    temp=0.
            do u=1,nx
              do v=1,ny
                temp=temp+basis(k-1,u,v)*stamp(i,u,v)
              enddo
            enddo
            do u=1,nx
              do v=1,ny
                stamp(i,u,v)=stamp(i,u,v)-temp*basis(k-1,u,v)
              enddo
            enddo	    
	  enddo

	  do i=k-1,n
            do u=1,nx
              do v=1,ny
                PCAstamp(i-k+2,u,v)=stamp(i,u,v)
              enddo
            enddo	    
            PCAweight(i-k+2)=weight(i)
	  enddo

c	  call find_PCA(np,n-k+2,nx,ny,PCAstamp,model)
	  call find_PCA_weighted(np,n-k+2,nx,ny,PCAstamp
     .,PCAweight,model)
	  
	  do i=1,k-1 	
	    temp=0.
            do u=1,nx
              do v=1,ny
                temp=temp+basis(i,u,v)*model(u,v)
              enddo
            enddo
            do u=1,nx
              do v=1,ny
                model(u,v)=model(u,v)-temp*basis(i,u,v)
              enddo
            enddo
	  enddo
	  call normalize_image(nx,ny,model)
          do u=1,nx
            do v=1,ny
	      basis(k,u,v)=model(u,v)
            enddo
          enddo

	write(*,*) k	 
c        filename='star_basis.fits'		
c	call show_stamps(np,k,nx,ny,basis,nx*6,ny*6,filename)
c       filename='star_residual.fits'		
c	call show_stamps(np,n,nx,ny,stamp,nx*6,ny*6,filename)
c	pause
	  
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gen_orthonormal_basis_mix(np,n,nx,ny,stamp
     .,nmax,basis)
	implicit none

	integer np,n,nx,ny,nmax
	real stamp(np,nx,ny),basis(np,nx,ny),coe(np)
	real PCAstamp(np,nx,ny),PCAweight(np)

	real sample(nx,ny),model(nx,ny),diff(nx,ny)
	real temp,weight(np),norm,temp1
	integer i,j,u,v,k,cx,cy,ipca,m

	character filename*80

	cx=nx/2+1
	cy=ny/2+1

	do i=1,n
	  weight(i)=stamp(i,cx,cy)**0.5
	  temp=1./stamp(i,cx,cy)
	  do u=1,nx
	    do v=1,ny
	      stamp(i,u,v)=stamp(i,u,v)*temp
	    enddo
	  enddo
	enddo
	
	temp=0.
	do i=1,n
	  temp=temp+weight(i)
	enddo
	temp=1./temp
	do i=1,n
	  weight(i)=weight(i)*temp
	enddo

        do u=1,nx
          do v=1,ny
            model(u,v)=0.
  	    do i=1,n 	
              model(u,v)=model(u,v)+weight(i)*stamp(i,u,v)
            enddo
          enddo
        enddo
	call normalize_image(nx,ny,model)
        do u=1,nx
          do v=1,ny
            basis(1,u,v)=model(u,v)
          enddo
        enddo

	do k=2,nmax
	  do i=1,n 	
	    temp=0.
            do u=1,nx
              do v=1,ny
                temp=temp+basis(k-1,u,v)*stamp(i,u,v)
              enddo
            enddo
            do u=1,nx
              do v=1,ny
                stamp(i,u,v)=stamp(i,u,v)-temp*basis(k-1,u,v)
              enddo
            enddo	    
	    ipca=k-1
	  enddo

	  if (k.gt.5) then	
 	    do i=k-1,n
              do u=1,nx
                do v=1,ny
                  PCAstamp(i-k+2,u,v)=stamp(i,u,v)
                enddo
              enddo	    
              PCAweight(i-k+2)=weight(i)
c              PCAweight(i-k+2)=1.
	    enddo
c	    call find_PCA(np,n-k+2,nx,ny,PCAstamp,model)
	    call find_PCA_weighted(np,n-k+2,nx,ny,PCAstamp
     .,PCAweight,model)
	  else
            do u=1,nx
              do v=1,ny
	        model(u,v)=stamp(ipca,u,v)
              enddo
            enddo
	  endif
	  
	  do i=1,k-1 	
	    temp=0.
            do u=1,nx
              do v=1,ny
                temp=temp+basis(i,u,v)*model(u,v)
              enddo
            enddo
            do u=1,nx
              do v=1,ny
                model(u,v)=model(u,v)-temp*basis(i,u,v)
              enddo
            enddo
	  enddo
	  call normalize_image(nx,ny,model)
          do u=1,nx
            do v=1,ny
	      basis(k,u,v)=model(u,v)
            enddo
          enddo
	 
	write(*,*) k
c        filename='star_basis.fits'		
c	call show_stamps(np,k,nx,ny,basis,nx*6,ny*6,filename)
c       filename='star_residual.fits'		
c	call show_stamps(np,n,nx,ny,stamp,nx*6,ny*6,filename)
c	pause
	  
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine normalize_image(nx,ny,image)
	implicit none
	
	integer nx,ny
	real image(nx,ny)
	integer u,v
	real norm

	  norm=0.
          do u=1,nx
            do v=1,ny
	      norm=norm+image(u,v)**2
            enddo
          enddo
	  norm=1./sqrt(norm)

          do u=1,nx
            do v=1,ny
	      image(u,v)=image(u,v)*norm
            enddo
          enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image55(nx,ny,map)
	implicit none

	integer nx,ny
	real map(nx,ny),temp(nx,ny),f(5,5)
	integer i,j,u,v

	do i=3,nx-2
	  do j=3,ny-2
	    do u=1,5
	      do v=1,5
	        f(u,v)=map(i+u-3,j+v-3)
	      enddo
	    enddo
	    call smooth_grid55(f)
	    temp(i,j)=f(3,3)
	    if (i.eq.3) then
	      temp(1,j)=f(1,3)
	      temp(2,j)=f(2,3)
	    endif
	    if (i.eq.nx-2) then
	      temp(nx,j)=f(5,3)
	      temp(nx-1,j)=f(4,3)
	    endif
	    if (j.eq.3) then
	      temp(i,1)=f(3,1)
	      temp(i,2)=f(3,2)
	    endif
	    if (j.eq.ny-2) then
	      temp(i,ny)=f(3,5)
	      temp(i,ny-1)=f(3,4)
	    endif
	    if (i.eq.3.and.j.eq.3) then
	      temp(1,1)=f(1,1)
	      temp(1,2)=f(1,2)
	      temp(2,1)=f(2,1)
	      temp(2,2)=f(2,2)
	    endif
	    if (i.eq.3.and.j.eq.ny-2) then
	      temp(1,ny)=f(1,5)
	      temp(1,ny-1)=f(1,4)
	      temp(2,ny)=f(2,5)
	      temp(2,ny-1)=f(2,4)
	    endif
	    if (i.eq.nx-2.and.j.eq.3) then
	      temp(nx,1)=f(5,1)
	      temp(nx-1,1)=f(4,1)
	      temp(nx,2)=f(5,2)
	      temp(nx-1,2)=f(4,2)
	    endif
	    if (i.eq.nx-2.and.j.eq.ny-2) then
	      temp(nx,ny)=f(5,5)
	      temp(nx-1,ny)=f(4,5)
	      temp(nx,ny-1)=f(5,4)
	      temp(nx-1,ny-1)=f(4,4)
	    endif
	  enddo
	enddo
	
	do i=1,nx
	  do j=1,ny
	    map(i,j)=temp(i,j)
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image55ln(nx,ny,map,n)
	implicit none

	integer nx,ny,n
	real map(nx,ny)
	integer i,j

	do i=1,nx
	  do j=1,ny
	    map(i,j)=log(map(i,j))
	  enddo
	enddo

	do i=1,n
	  call smooth_image55(nx,ny,map)
	enddo

	do i=1,nx
	  do j=1,ny
	    map(i,j)=exp(map(i,j))
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image33ln(nx,ny,map,n)
	implicit none

	integer nx,ny,n
	real map(nx,ny)
	integer i,j

	do i=1,nx
	  do j=1,ny
	    map(i,j)=log(map(i,j))
	  enddo
	enddo

	do i=1,n
	  call smooth_image33(nx,ny,map)
	enddo

	do i=1,nx
	  do j=1,ny
	    map(i,j)=exp(map(i,j))
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image33(nx,ny,map)
	implicit none

	integer nx,ny
	real map(nx,ny),temp(nx,ny),f(3,3)
	integer i,j,u,v

	do i=2,nx-1
	  do j=2,ny-1
	    do u=1,3
	      do v=1,3
	        f(u,v)=map(i+u-2,j+v-2)
	      enddo
	    enddo
	    call smooth_grid33(f)
	    temp(i,j)=f(2,2)
	    if (i.eq.2) temp(1,j)=f(1,2)
	    if (i.eq.nx-1) temp(nx,j)=f(3,2)
	    if (j.eq.2) temp(i,1)=f(2,1)
	    if (j.eq.ny-1) temp(i,ny)=f(2,3)
	    if (i.eq.2.and.j.eq.2) temp(1,1)=f(1,1)
	    if (i.eq.2.and.j.eq.ny-1) temp(1,ny)=f(1,3)
	    if (i.eq.nx-1.and.j.eq.2) temp(nx,1)=f(3,1)
	    if (i.eq.nx-1.and.j.eq.ny-1) temp(nx,ny)=f(3,3)
	  enddo
	enddo
	
	do i=1,nx
	  do j=1,ny
	    map(i,j)=temp(i,j)
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_grid33(f)
	implicit none

	real f(3,3)
	real vec(6),matx(6,6),x(3),y(3)
	real matx_1(6,6),a(6)

	integer i,j

	integer first
        SAVE first,matx_1
        DATA first /0/
	

	if (first.eq.0) then
 
	  first=1
          matx(1,1)=9. 
          matx(1,2)=0.    
          matx(1,3)=0.    
          matx(1,4)=6.    
          matx(1,5)=0.    
          matx(1,6)=6.    

          matx(2,1)=0. 
          matx(2,2)=6.    
          matx(2,3)=0.    
          matx(2,4)=0.    
          matx(2,5)=0.    
          matx(2,6)=0.    

          matx(3,1)=0. 
          matx(3,2)=0.    
          matx(3,3)=6.    
          matx(3,4)=0.    
          matx(3,5)=0.    
          matx(3,6)=0.    

          matx(4,1)=6. 
          matx(4,2)=0.    
          matx(4,3)=0.    
          matx(4,4)=6.    
          matx(4,5)=0.    
          matx(4,6)=4.    

          matx(5,1)=0. 
          matx(5,2)=0.    
          matx(5,3)=0.    
          matx(5,4)=0.    
          matx(5,5)=4.    
          matx(5,6)=0.    

          matx(6,1)=6. 
          matx(6,2)=0.    
          matx(6,3)=0.    
          matx(6,4)=4.    
          matx(6,5)=0.    
          matx(6,6)=6.    
	
   	  call matrix_inverse(matx,6,matx_1)
	endif

	do i=1,3	
	  x(i)=-2+i	
	  y(i)=-2+i	
	enddo

	do i=1,6
	  vec(i)=0.
	enddo

	do i=1,3
	  do j=1,3
	    vec(1)=vec(1)+f(i,j)	    
	    vec(2)=vec(2)+f(i,j)*x(i)	    
	    vec(3)=vec(3)+f(i,j)*y(j)	    
	    vec(4)=vec(4)+f(i,j)*x(i)*x(i)	    
	    vec(5)=vec(5)+f(i,j)*x(i)*y(j)	    
	    vec(6)=vec(6)+f(i,j)*y(j)*y(j)	
	  enddo
	enddo

	do i=1,6
	  a(i)=0.
	  do j=1,6
	    a(i)=a(i)+matx_1(i,j)*vec(j)
	  enddo
	enddo

	do i=1,3
	  do j=1,3
	    f(i,j)=a(1)+a(2)*x(i)+a(3)*y(j)+a(4)*x(i)*x(i)
     .+a(5)*x(i)*y(j)+a(6)*y(j)*y(j)	    
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_grid55(f)
	implicit none

	real f(5,5)
	real vec(6),matx(6,6),x(5),y(5)
	real matx_1(6,6),a(6)

	integer i,j

	integer first
        SAVE first,matx_1
        DATA first /0/
	

	if (first.eq.0) then
 
	  first=1
          matx(1,1)=25. 
          matx(1,2)=0.    
          matx(1,3)=0.    
          matx(1,4)=50.    
          matx(1,5)=0.    
          matx(1,6)=50.    

          matx(2,1)=0. 
          matx(2,2)=50.    
          matx(2,3)=0.    
          matx(2,4)=0.    
          matx(2,5)=0.    
          matx(2,6)=0.    

          matx(3,1)=0. 
          matx(3,2)=0.    
          matx(3,3)=50.    
          matx(3,4)=0.    
          matx(3,5)=0.    
          matx(3,6)=0.    

          matx(4,1)=50. 
          matx(4,2)=0.    
          matx(4,3)=0.    
          matx(4,4)=170.    
          matx(4,5)=0.    
          matx(4,6)=100.    

          matx(5,1)=0. 
          matx(5,2)=0.    
          matx(5,3)=0.    
          matx(5,4)=0.    
          matx(5,5)=100.    
          matx(5,6)=0.    

          matx(6,1)=50. 
          matx(6,2)=0.    
          matx(6,3)=0.    
          matx(6,4)=100.    
          matx(6,5)=0.    
          matx(6,6)=170.    
	
   	  call matrix_inverse(matx,6,matx_1)
	endif

	do i=1,5	
	  x(i)=-3+i	
	  y(i)=-3+i	
	enddo

	do i=1,6
	  vec(i)=0.
	enddo

	do i=1,5
	  do j=1,5
	    vec(1)=vec(1)+f(i,j)	    
	    vec(2)=vec(2)+f(i,j)*x(i)	    
	    vec(3)=vec(3)+f(i,j)*y(j)	    
	    vec(4)=vec(4)+f(i,j)*x(i)*x(i)	    
	    vec(5)=vec(5)+f(i,j)*x(i)*y(j)	    
	    vec(6)=vec(6)+f(i,j)*y(j)*y(j)	
	  enddo
	enddo

	do i=1,6
	  a(i)=0.
	  do j=1,6
	    a(i)=a(i)+matx_1(i,j)*vec(j)
	  enddo
	enddo

	do i=1,5
	  do j=1,5
	    f(i,j)=a(1)+a(2)*x(i)+a(3)*y(j)+a(4)*x(i)*x(i)
     .+a(5)*x(i)*y(j)+a(6)*y(j)*y(j)	    
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_grid55_sig(f)
	implicit none

	real f(5,5)
	real vec(6),matx(6,6),x(5),y(5)
	real matx_1(6,6),a(6)

	integer i,j

	integer first
        SAVE first,matx_1
        DATA first /0/
	

	if (first.eq.0) then
 
	  first=1
          matx(1,1)=25. 
          matx(1,2)=0.    
          matx(1,3)=0.    
          matx(1,4)=50.    
          matx(1,5)=0.    
          matx(1,6)=50.    

          matx(2,1)=0. 
          matx(2,2)=50.    
          matx(2,3)=0.    
          matx(2,4)=0.    
          matx(2,5)=0.    
          matx(2,6)=0.    

          matx(3,1)=0. 
          matx(3,2)=0.    
          matx(3,3)=50.    
          matx(3,4)=0.    
          matx(3,5)=0.    
          matx(3,6)=0.    

          matx(4,1)=50. 
          matx(4,2)=0.    
          matx(4,3)=0.    
          matx(4,4)=170.    
          matx(4,5)=0.    
          matx(4,6)=100.    

          matx(5,1)=0. 
          matx(5,2)=0.    
          matx(5,3)=0.    
          matx(5,4)=0.    
          matx(5,5)=100.    
          matx(5,6)=0.    

          matx(6,1)=50. 
          matx(6,2)=0.    
          matx(6,3)=0.    
          matx(6,4)=100.    
          matx(6,5)=0.    
          matx(6,6)=170.    
	
   	  call matrix_inverse(matx,6,matx_1)
	endif

	do i=1,5	
	  x(i)=-3+i	
	  y(i)=-3+i	
	enddo

	do i=1,6
	  vec(i)=0.
	enddo

	do i=1,5
	  do j=1,5
	    vec(1)=vec(1)+f(i,j)	    
	    vec(2)=vec(2)+f(i,j)*x(i)	    
	    vec(3)=vec(3)+f(i,j)*y(j)	    
	    vec(4)=vec(4)+f(i,j)*x(i)*x(i)	    
	    vec(5)=vec(5)+f(i,j)*x(i)*y(j)	    
	    vec(6)=vec(6)+f(i,j)*y(j)*y(j)	
	  enddo
	enddo

	do i=1,6
	  a(i)=0.
	  do j=1,6
	    a(i)=a(i)+matx_1(i,j)*vec(j)
	  enddo
	enddo

	do i=1,5
	  do j=1,5
	    f(i,j)=a(1)+a(2)*x(i)+a(3)*y(j)+a(4)*x(i)*x(i)
     .+a(5)*x(i)*y(j)+a(6)*y(j)*y(j)	    
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine test1()

	real f(3,3)
	integer i,j
	real ran1

	do i=1,3
	  do j=1,3
	    f(i,j)=4.+i**2+j
	    write(*,*) i,j,f(i,j)
	  enddo
	enddo

	call smooth_grid33(f)

	do i=1,3
	  do j=1,3
	    write(*,*) i,j,f(i,j)
	  enddo
	enddo

	do i=1,3
	  do j=1,3
	    f(i,j)=4.+i**2+j+(ran1()-0.5)*0.5
	    write(*,*) i,j,f(i,j)
	  enddo
	enddo
	call smooth_grid33(f)

	do i=1,3
	  do j=1,3
	    write(*,*) i,j,f(i,j)
	  enddo
	enddo

	pause

	return
	end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine test()

	real matx(3,3),val(3),vec(3,3)
	integer i,j

	matx(1,1)=1.
	matx(2,2)=1.
	matx(3,3)=1.
	matx(1,2)=0.5
	matx(1,3)=0.5
	matx(2,3)=0.5
	matx(2,1)=matx(1,2)
	matx(3,1)=matx(1,3)
	matx(3,2)=matx(2,3)
	
	call find_eigen_vec(3,3,matx,val,vec)

	do i=1,3
	  write(*,*) i,val(i),vec(1,i),vec(2,i),vec(3,i)
	enddo

	return
	end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine check_linear_dependence(np,n,nx,ny,stamp
     .,coe,stacked)
	implicit none

	integer np,n,nx,ny
	real stamp(np,nx,ny),coe(np)
	real cov(n-1,n-1),cov_1(n-1,n-1),vec(n-1)
	real stacked(nx,ny),temp(nx,ny)
	integer i,j,u,v

	integer cx,cy
	
	cx=nx/2+1
	cy=ny/2+1

	do u=1,nx
	  do v=1,ny
	    temp(u,v)=stamp(1,u,v)**(-1)
	  enddo
	enddo

	do i=1,n-1
	  do j=1,n-1	    
	    cov(i,j)=0.
	    do u=1,nx
	      do v=1,ny
	        if (u.ne.cx.or.v.ne.cy) then	    
	          cov(i,j)=cov(i,j)
     .+stamp(i+1,u,v)*stamp(j+1,u,v)*temp(u,v)
	        endif
	      enddo
	    enddo
	  enddo
	enddo

	do i=1,n-1
          vec(i)=0.
          do u=1,nx
            do v=1,ny
              if (u.ne.cx.or.v.ne.cy) then	    
                vec(i)=vec(i)-stamp(1,u,v)*stamp(i+1,u,v)*temp(u,v)
	      endif
            enddo
          enddo
	enddo
		
	call matrix_inverse(cov,n-1,cov_1)
	
	do i=1,n-1
	  coe(i)=0
	  do j=1,n-1
	    coe(i)=coe(i)+cov_1(i,j)*vec(j)
	  enddo
	enddo
	
	do u=1,nx
	  do v=1,ny
	    stacked(u,v)=stamp(1,u,v)
	    do i=1,n-1
	      stacked(u,v)=stacked(u,v)+stamp(i+1,u,v)*coe(i)	
	    enddo
	  enddo
	enddo
	  	

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

