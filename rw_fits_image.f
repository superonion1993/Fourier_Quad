        subroutine map_points(np,n,x,xc,yc,nx,ny,imagename)
        implicit none

        integer np,n,nx,ny
        real x(np,2),ratio,xc,yc

        real ymin,ymax,xmin,xmax,dx,dy,temp
        real map(nx,ny),yratio,yoffset,xratio,xoffset

        integer i,j
	real x1,x2,y1,y2,xx,yy
	character imagename*80

        xmax=x(1,1)
	xmin=xmax
	ymax=x(1,2)
	ymin=ymax

        do i=2,n
          if (x(i,1).gt.xmax) xmax=x(i,1)
          if (x(i,1).lt.xmin) xmin=x(i,1)
          if (x(i,2).gt.ymax) ymax=x(i,2)
          if (x(i,2).lt.ymin) ymin=x(i,2)
        enddo

c	dx=xmax-xmin
c	dy=ymax-ymin

	ymax=ymax/3.

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        do i=1,n	
          xx=(x(i,1)-xoffset)*xratio+nx/2
          yy=(x(i,2)-yoffset)*yratio+ny/2
          call draw_dot(nx,ny,map,xx,yy,100.,1.)
        enddo	

        x1=(xc-xoffset)*xratio+nx/2
        x2=x1
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(yc-yoffset)*yratio+ny/2
        y2=(yc-yoffset)*yratio+ny/2
        call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        call reverse_color(nx,ny,map)
        call writeimage(imagename,nx,ny,nx,ny,map)	  

        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine map_function(func,np,n,x,nx,ny,imagename)
        implicit none

c Plot monotonic function & its data

        integer np,n,nx,ny
        real x(np,4),func
	character imagename*80
	external func

	real ratio,ymin,ymax,xmin,xmax,dx,dy,x1,x2,y1,y2,xc,yc,sig
	integer nm,i
	parameter (nm=100)
        real map(nx,ny),model(nm,2)
        real yratio,yoffset,xratio,xoffset

        xmax=x(1,1)
	xmin=x(1,1)
        do i=1,n
	  if (x(i,4).gt.0) then
            if (x(i,1).gt.xmax) xmax=x(i,1)
            if (x(i,1).lt.xmin) xmin=x(i,1)
	  endif
        enddo
	dx=xmax-xmin
	  
        ymax=x(1,2)
	ymin=x(1,2)
        do i=1,n
	  if (x(i,4).gt.0) then
            if (x(i,2)+x(i,3).gt.ymax) ymax=x(i,2)+x(i,3)
            if (x(i,2)-x(i,3).lt.ymin) ymin=x(i,2)-x(i,3)
	  endif
        enddo
	ymax=max(ymax,func(xmax))
	ymin=min(ymin,func(xmin))

	dy=ymax-ymin

	xmax=xmax+dx*0.1
	xmin=xmin-dx*0.1
	
	ymax=ymax+dy*0.1
	ymin=ymin-dy*0.1

	dx=(xmax-xmin)/(nm-1.)	
	do i=1,nm
	  model(i,1)=xmin+dx*(i-1.)
	  model(i,2)=func(model(i,1))
	enddo

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        do i=1,n	
	  if (x(i,4).gt.0) then
            xc=(x(i,1)-xoffset)*xratio+nx/2
            yc=(x(i,2)-yoffset)*yratio+ny/2
            sig=x(i,3)*yratio
            call draw_error(nx,ny,map,xc,yc,sig,180.,1.)
	  endif
        enddo	

        do i=1,nm-1
          x1=(model(i,1)-xoffset)*xratio+nx/2
          x2=(model(i+1,1)-xoffset)*xratio+nx/2
          y1=(model(i,2)-yoffset)*yratio+ny/2
          y2=(model(i+1,2)-yoffset)*yratio+ny/2
          call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      
        enddo	

        x1=(0.-xoffset)*xratio+nx/2
        x2=x1
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(0.-yoffset)*yratio+ny/2
        y2=y1
        call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      
	 
        call reverse_color(nx,ny,map)
        call writeimage(imagename,nx,ny,nx,ny,map)	  
	   

        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine map_stars(np,n,n1,n2,star,x,nx,ny,imagename)
        implicit none

        integer np,n,nx,ny,n1,n2
        real x(np,2),star(np,n1,n2)

        real map(nx,ny)

        integer i,j,u,v,x1,x2,y1,y2,n1_2,n2_2
	character imagename*80

	do i=1,nx
	  do j=1,ny
	    map(i,j)=1e-6
	  enddo
	enddo

	n1_2=n1/2
	n2_2=n2/2

	do i=1,n
	  x1=int(x(i,1)+1)
	  x2=x1+n1-1
	  y1=int(x(i,2)+1)
	  y2=y1+n2-1
	  do u=x1,x2
	    do v=y1,y2
	      map(u,v)=star(i,u-x1+1,v-y1+1)
	    enddo
	  enddo
	enddo

        call writeimage(imagename,nx,ny,nx,ny,map)	  
	  

        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine data_binning(np,n,x,a,b,dx,imagename)
	implicit none

	integer np,n,nbin
	real x(np),a,b,dx

	character imagename*80
	integer NMAX
	parameter (NMAX=200)
	real num(NMAX),xbin(NMAX)

        real ymin,ymax,xmin,xmax,ratio
	
	integer nx,ny
	parameter (nx=500)
	parameter (ny=500)
        real map(nx,ny)

        real yratio,yoffset,xratio,xoffset
        real intensity,thickness

        integer i,ibin
	real x1,x2,y1,y2

	xmin=a
	xmax=b

	nbin=int((xmax-xmin)/dx+1)
	do ibin=1,nbin+1
	  xbin(ibin)=xmin+dx*(ibin-1.)
	enddo

	ymin=0
	ymax=0
	do ibin=1,nbin
	  num(ibin)=0.
          do i=1,n
	    if (x(i).ge.xbin(ibin).and.x(i).lt.xbin(ibin+1)) then
	      num(ibin)=num(ibin)+1
	    endif 
	    if (x(i).gt.xmax) then
	      pause 'need to increase box size to:'
	      write(*,*) x(i)
	    endif
          enddo
	  ymax=max(ymax,num(ibin))
	enddo

	ymax=ymax*1.2

	xmin=a
	xmax=b

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)
        intensity=100.
        thickness=0.1
       
        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        do ibin=1,nbin
          x1=(xbin(ibin)-xoffset)*xratio+nx/2
          x2=(xbin(ibin+1)-xoffset)*xratio+nx/2
          y1=(ymin-yoffset)*yratio+ny/2
          y2=(num(ibin)-yoffset)*yratio+ny/2
          call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      
        enddo	

        call reverse_color(nx,ny,map)
        call writeimage(imagename,nx,ny,nx,ny,map)	  

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
        subroutine make_scatter_map(np,n,x,nbin,nx,ny,imagename,a,b)
        implicit none

        integer np,n,nx,ny,nb,ibin,nbin
        real x(np,2),ratio,a,b
	integer NMAX5
	parameter (NMAX5=1000000)

	real xbin(nbin+1),y(NMAX5)
	real ymed(nbin),ylow(nbin),yhigh(nbin)
	integer counts(nbin)

	real num(nbin),mean(nbin),sig(nbin),fitx(nbin)
        real ymin,ymax,xmin,xmax,dx,dy,temp
        real map(nx,ny)
	real xxmin,yymin,xxmax,yymax

        real yratio,yoffset,xratio,xoffset
        real intensity,thickness,total

        integer i,ix,iy,ix1,ix2,iy1,iy2,j
	real x1,x2,y1,y2,xx,yy
	character imagename*80

        xmax=x(1,1)
	xmin=x(1,1)

        do i=1,n
          if (x(i,1).gt.xmax) xmax=x(i,1)
          if (x(i,1).lt.xmin) xmin=x(i,1)
        enddo

	dx=(xmax-xmin)/nbin
	do ibin=1,nbin+1
	  xbin(ibin)=xmin+dx*(ibin-1.)
	enddo


	do ibin=1,nbin

  	  goto 40

	  counts(ibin)=0
          do i=1,n
	    if (x(i,1).ge.xbin(ibin).and.x(i,1).le.xbin(ibin+1)) then
	      counts(ibin)=counts(ibin)+1
	      y(counts(ibin))=x(i,2)
	    endif 
          enddo
c	  write(*,*) ibin,counts(ibin) 
	  call sort(counts(ibin),NMAX5,y)	
	  ymed(ibin)=y(int(counts(ibin)/2+0.5))
	  ylow(ibin)=y(int(counts(ibin)*0.16+0.5))
	  yhigh(ibin)=y(int(counts(ibin)*0.84+0.5))

	  write(*,*) ibin,counts(ibin),0.5*(xbin(ibin)+xbin(ibin+1))
     .,ymed(ibin),yhigh(ibin)-ylow(ibin) 

	  mean(ibin)=ymed(ibin)
	  sig(ibin)=(yhigh(ibin)-ylow(ibin))*0.5

	  goto 30

40	  num(ibin)=0.
	  mean(ibin)=0.
	  sig(ibin)=0.
          do i=1,n
	    if (x(i,1).ge.xbin(ibin).and.x(i,1).le.xbin(ibin+1)) then
	      num(ibin)=num(ibin)+1
	      mean(ibin)=mean(ibin)+x(i,2)
	      sig(ibin)=sig(ibin)+x(i,2)*x(i,2)
	    endif 
          enddo
	  mean(ibin)=mean(ibin)/num(ibin)
	  sig(ibin)
     .=sqrt((sig(ibin)/num(ibin)-mean(ibin)*mean(ibin))/num(ibin))

c	  write(*,*) ibin,num(ibin),exp(0.5*(xbin(ibin)+xbin(ibin+1)))
c     .,mean(ibin),sig(ibin) 
	  write(*,*) ibin,num(ibin),0.5*(xbin(ibin)+xbin(ibin+1))
     .,mean(ibin),sig(ibin) 

	  ymed(ibin)=mean(ibin)
	  ylow(ibin)=mean(ibin)-sig(ibin)
	  yhigh(ibin)=mean(ibin)+sig(ibin)

30	enddo

	do i=1,nbin
	  fitx(i)=0.5*(xbin(i)+xbin(i+1))
	enddo

	call fit_line2(nbin,fitx,mean,sig,a,b)

	ymin=ylow(1)
	ymax=yhigh(1)	
	do ibin=1,nbin
	  if (ylow(ibin).lt.ymin) ymin=ylow(ibin)
	  if (yhigh(ibin).gt.ymax) ymax=yhigh(ibin)
	enddo

	dy=ymax-ymin

C	if (ymax.lt.0) ymax=dy*0.1
C	if (ymin.gt.0) ymin=-dy*0.1

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)
        intensity=100.
        thickness=0.1
       

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

         do ibin=1,nbin-1
           x1=((xbin(ibin)+xbin(ibin+1))/2.-xoffset)*xratio+nx/2
           x2=((xbin(ibin+1)+xbin(ibin+2))/2.-xoffset)*xratio+nx/2
           y1=(ylow(ibin)-yoffset)*yratio+ny/2
           y2=(ylow(ibin+1)-yoffset)*yratio+ny/2
           call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
           y1=(yhigh(ibin)-yoffset)*yratio+ny/2
           y2=(yhigh(ibin+1)-yoffset)*yratio+ny/2
           call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
           y1=(ymed(ibin)-yoffset)*yratio+ny/2
           y2=(ymed(ibin+1)-yoffset)*yratio+ny/2
           call draw_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
         enddo	

         x1=(0.-xoffset)*xratio+nx/2
 	 x2=x1
         y1=(ymin-yoffset)*yratio+ny/2
         y2=(ymax-yoffset)*yratio+ny/2
         call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      

         x1=(xmin-xoffset)*xratio+nx/2
         x2=(xmax-xoffset)*xratio+nx/2
         y1=(0.-yoffset)*yratio+ny/2
         y2=y1
         call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      
	 
         x1=(xmin-xoffset)*xratio+nx/2
         x2=(xmax-xoffset)*xratio+nx/2
         y1=(a+xmin*b-yoffset)*yratio+ny/2
         y2=(a+xmax*b-yoffset)*yratio+ny/2
         call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      


         call reverse_color(nx,ny,map)
         call writeimage(imagename,nx,ny,nx,ny,map)	  
	   

         return
         END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine make_scatter_map2(np,n,x,nbin,nx,ny,imagename)
        implicit none

        integer np,n,nx,ny,nb,ibin,nbin
        real x(np,2),ratio
	integer NMAX5
	parameter (NMAX5=1000000)

	real xbin(nbin+1),y(NMAX5)
	real ymed(nbin),ylow(nbin),yhigh(nbin)
	integer counts(nbin)

	real num(nbin),mean(nbin),sig(nbin),fitx(nbin)
        real ymin,ymax,xmin,xmax,dx,dy,temp
        real map(nx,ny)
	real xxmin,yymin,xxmax,yymax

        real yratio,yoffset,xratio,xoffset
        real intensity,thickness,total

        integer i,ix,iy,ix1,ix2,iy1,iy2,j
	real x1,x2,y1,y2,xx,yy
	character imagename*80

        xmax=x(1,1)
	xmin=x(1,1)

        do i=1,n
          if (x(i,1).gt.xmax) xmax=x(i,1)
          if (x(i,1).lt.xmin) xmin=x(i,1)
        enddo

	dx=(xmax-xmin)/nbin
	do ibin=1,nbin+1
	  xbin(ibin)=xmin+dx*(ibin-1.)
	enddo


	do ibin=1,nbin

c  	  goto 40

	  counts(ibin)=0
          do i=1,n
	    if (x(i,1).ge.xbin(ibin).and.x(i,1).le.xbin(ibin+1)) then
	      counts(ibin)=counts(ibin)+1
	      y(counts(ibin))=x(i,2)
	    endif 
          enddo
	  call sort(counts(ibin),NMAX5,y)	
	  ymed(ibin)=y(int(counts(ibin)/2+0.5))
	  ylow(ibin)=y(int(counts(ibin)*0.16+0.5))
	  yhigh(ibin)=y(int(counts(ibin)*0.84+0.5))

	  write(*,*) ibin,counts(ibin),0.5*(xbin(ibin)+xbin(ibin+1))
     .,ymed(ibin),yhigh(ibin)-ylow(ibin) 

	  goto 30

40	  num(ibin)=0.
	  mean(ibin)=0.
	  sig(ibin)=0.
          do i=1,n
	    if (x(i,1).ge.xbin(ibin).and.x(i,1).le.xbin(ibin+1)) then
	      num(ibin)=num(ibin)+1
	      mean(ibin)=mean(ibin)+x(i,2)
	      sig(ibin)=sig(ibin)+x(i,2)*x(i,2)
	    endif 
          enddo
	  mean(ibin)=mean(ibin)/num(ibin)
	  sig(ibin)
     .=sqrt((sig(ibin)/num(ibin)-mean(ibin)*mean(ibin))/num(ibin))

c	  write(*,*) ibin,num(ibin),exp(0.5*(xbin(ibin)+xbin(ibin+1)))
c     .,mean(ibin),sig(ibin) 
	  write(*,*) ibin,num(ibin),0.5*(xbin(ibin)+xbin(ibin+1))
     .,mean(ibin),sig(ibin) 

	  ymed(ibin)=mean(ibin)
	  ylow(ibin)=mean(ibin)-sig(ibin)
	  yhigh(ibin)=mean(ibin)+sig(ibin)

30	enddo

	ymin=ylow(1)
	ymax=yhigh(1)	
	do ibin=1,nbin
	  if (ylow(ibin).lt.ymin) ymin=ylow(ibin)
	  if (yhigh(ibin).gt.ymax) ymax=yhigh(ibin)
	enddo

	dy=ymax-ymin
	if (ymax.lt.0) ymax=dy*0.3
	if (ymin.gt.0) ymin=-dy*0.3

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)
        intensity=100.
        thickness=0.1
       

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        do i=1,n	
           xx=(x(i,1)-xoffset)*xratio+nx/2
           yy=(x(i,2)-yoffset)*yratio+ny/2
           call draw_dot(nx,ny,map,xx,yy,100.,1.)
        enddo	

         do ibin=1,nbin-1
           x1=((xbin(ibin)+xbin(ibin+1))/2.-xoffset)*xratio+nx/2
           x2=((xbin(ibin+1)+xbin(ibin+2))/2.-xoffset)*xratio+nx/2
           y1=(ylow(ibin)-yoffset)*yratio+ny/2
           y2=(ylow(ibin+1)-yoffset)*yratio+ny/2
           call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
           y1=(yhigh(ibin)-yoffset)*yratio+ny/2
           y2=(yhigh(ibin+1)-yoffset)*yratio+ny/2
           call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
           y1=(ymed(ibin)-yoffset)*yratio+ny/2
           y2=(ymed(ibin+1)-yoffset)*yratio+ny/2
           call draw_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
         enddo	

         x1=(0.-xoffset)*xratio+nx/2
 	 x2=x1
         y1=(ymin-yoffset)*yratio+ny/2
         y2=(ymax-yoffset)*yratio+ny/2
         call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      

         x1=(xmin-xoffset)*xratio+nx/2
         x2=(xmax-xoffset)*xratio+nx/2
         y1=(0.-yoffset)*yratio+ny/2
         y2=y1
         call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      
	 
         call reverse_color(nx,ny,map)
         call writeimage(imagename,nx,ny,nx,ny,map)	  
	   

         return
         END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine make_scatter_map3(np,n,x,f,nbin,nx,ny,imagename)
        implicit none

        integer np,n,nx,ny,nb,ibin,nbin
        integer x(np)
	real f(np),ratio
	integer NMAX5
	parameter (NMAX5=1000000)

	real y(NMAX5)
	real ymed(nbin),ylow(nbin),yhigh(nbin)
	integer counts(nbin)

	real num(nbin),mean(nbin),sig(nbin),fitx(nbin)
        real ymin,ymax,xmin,xmax,dx,dy,temp
        real map(nx,ny)
	real xxmin,yymin,xxmax,yymax

        real yratio,yoffset,xratio,xoffset
        real intensity,thickness,total

        integer i,ix,iy,ix1,ix2,iy1,iy2,j
	real x1,x2,y1,y2,xx,yy
	character imagename*80


	do ibin=2,nbin

	  counts(ibin)=0
          do i=1,n
	    if (x(i).eq.ibin) then
	      counts(ibin)=counts(ibin)+1
	      y(counts(ibin))=f(i)
	    endif 
          enddo
	  call sort(counts(ibin),NMAX5,y)	
	  ymed(ibin)=y(int(counts(ibin)/2+0.5))
	  ylow(ibin)=y(int(counts(ibin)*0.16+0.5))
	  yhigh(ibin)=y(int(counts(ibin)*0.84+0.5))

	  write(*,*) ibin,counts(ibin),ymed(ibin)
     .,yhigh(ibin)-ylow(ibin) 

	enddo

	ymin=ylow(1)
	ymax=yhigh(1)	
	do ibin=2,nbin
	  if (ylow(ibin).lt.ymin) ymin=ylow(ibin)
	  if (yhigh(ibin).gt.ymax) ymax=yhigh(ibin)
	enddo

	dy=ymax-ymin
c	ymax=ymax+dy*0.2
c	ymin=-dy*0.2

	xmin=0
	xmax=nbin+1

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)
        intensity=100.
        thickness=0.1
       

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        do i=1,n	
           xx=(x(i)-xoffset)*xratio+nx/2
           yy=(f(i)-yoffset)*yratio+ny/2
           call draw_dot(nx,ny,map,xx,yy,100.,1.)
        enddo	

         do ibin=2,nbin-1
           x1=(ibin-xoffset)*xratio+nx/2
           x2=(ibin+1.-xoffset)*xratio+nx/2
           y1=(ylow(ibin)-yoffset)*yratio+ny/2
           y2=(ylow(ibin+1)-yoffset)*yratio+ny/2
           call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
           y1=(yhigh(ibin)-yoffset)*yratio+ny/2
           y2=(yhigh(ibin+1)-yoffset)*yratio+ny/2
           call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
           y1=(ymed(ibin)-yoffset)*yratio+ny/2
           y2=(ymed(ibin+1)-yoffset)*yratio+ny/2
           call draw_line(nx,ny,map,x1,y1,x2,y2,255.,2.)      
         enddo	

         x1=(0.-xoffset)*xratio+nx/2
 	 x2=x1
         y1=(ymin-yoffset)*yratio+ny/2
         y2=(ymax-yoffset)*yratio+ny/2
         call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      

         x1=(xmin-xoffset)*xratio+nx/2
         x2=(xmax-xoffset)*xratio+nx/2
         y1=(0.-yoffset)*yratio+ny/2
         y2=y1
         call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      
	 
         call reverse_color(nx,ny,map)
         call writeimage(imagename,nx,ny,nx,ny,map)	  
	   

         return
         END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fit_line(np,n,x,a,b)
        implicit none

        integer np,n
        real x(np,2),a,b
	integer i

	real sumy,sumx,sumxx,sumxy

	sumy=0
	sumx=0
	sumxx=0
	sumxy=0

	do i=1,n
	  sumy=sumy+x(i,2)
	  sumx=sumx+x(i,1)
   	  sumxx=sumxx+x(i,1)*x(i,1)
 	  sumxy=sumxy+x(i,2)*x(i,2)
	enddo  	  
	
	b=(n*sumxy-sumx*sumy)/(n*sumxx-sumx*sumx)
	a=(sumy-b*sumx)/n

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fit_line2(nbin,x,y,sig,a,b)
        implicit none

        integer nbin
        real x(nbin),y(nbin),sig(nbin),a,b
	integer i

	real sumy,sumx,sumxx,sumxy,temp,sumn

	sumy=0
	sumx=0
	sumxx=0
	sumxy=0
	sumn=0.

	do i=1,nbin
	  temp=sig(i)**(-2)
	  sumn=sumn+temp
	  sumy=sumy+y(i)*temp
	  sumx=sumx+x(i)*temp
   	  sumxx=sumxx+x(i)*x(i)*temp
 	  sumxy=sumxy+x(i)*y(i)*temp
	enddo  	  
	
	b=(sumn*sumxy-sumx*sumy)/(sumn*sumxx-sumx*sumx)
	a=(sumy-b*sumx)/sumn

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine map_distribution(np,n,x,nx,ny,imagename)
        implicit none

        integer np,n,nx,ny
        real x(np,2),ratio

        real ymin,ymax,xmin,xmax,xx,yy
        real map(nx,ny)

        real yratio,yoffset,xratio,xoffset
        real intensity,thickness

        integer i,j
	real x1,x2,y1,y2
	character imagename*80

        xmax=x(1,1)
	xmin=x(1,1)
	ymax=x(1,2)
	ymin=x(1,2)
        do i=1,n
          if (x(i,1).gt.xmax) xmax=x(i,1)
          if (x(i,1).lt.xmin) xmin=x(i,1)
          if (x(i,2).gt.ymax) ymax=x(i,2)
          if (x(i,2).lt.ymin) ymin=x(i,2)
        enddo

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)
        intensity=100.
       
        do i=1,n	
           xx=(x(i,1)-xoffset)*xratio+nx/2
           yy=(x(i,2)-yoffset)*yratio+ny/2
           call draw_dot(nx,ny,map,xx,yy,100.,2.)
        enddo	

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        call reverse_color(nx,ny,map)
        call writeimage(imagename,nx,ny,nx,ny,map)	  
	  

        return
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc			
	subroutine read_stamps(np,nstart,n,nx,ny,stamp,n1,n2
     .,filename)
	implicit none

	integer np,nstart,n,nx,ny,n1,n2
	real stamp(np,nx,ny),large_stamp(2000,4000)
	character filename*80
	integer i,j,offx,offy,k

        call readimage(filename,n1,n2,2000,4000,large_stamp)

	offx=0	
	offy=0
	
	do k=nstart,n
	  do i=1,nx
	    do j=1,ny
	      stamp(k,i,j)=large_stamp(i+offx,j+offy)
	    enddo
	  enddo
	  offx=offx+nx
	  if (offx+nx.gt.n1) then
	    offx=0
	    offy=offy+ny
	  endif
	enddo


	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine write_stamps(np,nstart,n,nx,ny,stamp,n1,n2
     .,filename)
	implicit none

	integer np,nstart,n,nx,ny,n1,n2
	real stamp(np,nx,ny),large_stamp(2000,4000)
	character filename*80
	integer i,j,offx,offy,k

	do i=1,n1
	  do j=1,n2
	    large_stamp(i,j)=0.
	  enddo
	enddo

	offx=0	
	offy=0
	
	do k=nstart,n
          if (offy+ny.gt.n2) pause 'large_stamp is too small!!'
	  do i=1,nx
	    do j=1,ny
	      large_stamp(i+offx,j+offy)=stamp(k,i,j)
	    enddo
	  enddo
	  offx=offx+nx
	  if (offx+nx.gt.n1) then
	    offx=0
	    offy=offy+ny
	  endif
	enddo

        call writeimage(filename,n1,n2,2000,4000,large_stamp)

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine clean_map(nx,ny,map)
       implicit none
	   
       integer nx,ny,i,j
       real map(nx,ny)

       do i=1,nx
         do j=1,ny
           map(i,j)=0.
         enddo
       enddo
	   
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 subroutine draw_rectangle(nx,ny,map,x1,y1,x2,y2,intensity
     .,thickness)
	 implicit none
	
	 integer nx,ny
        real map(nx,ny)
        real x1,y1,x2,y2,intensity,thickness

	 call draw_line(nx,ny,map,x1,y1,x1,y2,intensity,thickness)
	 call draw_line(nx,ny,map,x2,y1,x2,y2,intensity,thickness)
	 call draw_line(nx,ny,map,x1,y1,x2,y1,intensity,thickness)
	 call draw_line(nx,ny,map,x1,y2,x2,y2,intensity,thickness)

	 return
	 end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	 subroutine reverse_color(nx,ny,map)
	 implicit none

	 integer nx,ny,i,j
	 real map(nx,ny)
	 
	 do i=1,nx
	   do j=1,ny
	     map(i,j)=255.-map(i,j)
	   enddo
        enddo
	 
	 return
	 end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
         subroutine draw_error(nx,ny,map,x,y,sig,intensity,thickness)
         implicit none
	   
  	 integer nx,ny
         real map(nx,ny)
         real x,y,sig,intensity,thickness
	 real x1,x2,y1,y2

	 x1=x
	 x2=x
	 y1=y-sig
	 y2=y+sig
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 x1=x-3.*thickness
	 x2=x+3.*thickness
	 y1=y
	 y2=y
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 y1=y+sig
	 y2=y+sig
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 y1=y-sig
	 y2=y-sig
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)
	   
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 subroutine draw_line(nx,ny,map,x1,y1,x2,y2,intensity
     .,thickness)
         implicit none
	   
  	  integer nx,ny
         real map(nx,ny)
         real x1,y1,x2,y2,intensity,thickness
	   
	 real r,dr,cosx,cosy,x,y	   
         integer ndots,i
	   
         r=sqrt((x2-x1)**2+(y2-y1)**2)	   	  

	   
         if (r.gt.0.) then
	   
           cosx=(x2-x1)/r
           cosy=(y2-y1)/r	
           ndots=int(r)	   
           dr=r/ndots	    
	      
           x=x1
           y=y1
           call draw_dot(nx,ny,map,x,y,intensity,thickness)
           do i=1,ndots
   	     x=x+dr*cosx
   	     y=y+dr*cosy
   	     call draw_dot(nx,ny,map,x,y,intensity,thickness)
 	   enddo	   	  	   
         endif	   
	   
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine draw_dashed_line(nx,ny,map,x1,y1,x2,y2
     .,intensity,thickness)
         implicit none
	   
  	  integer nx,ny
         real map(nx,ny)
         real x1,y1,x2,y2,intensity,thickness
	   
	  real r,dr,cosx,cosy,x,y	   
         integer ndots,i
	   
         r=sqrt((x2-x1)**2+(y2-y1)**2)	   	  

	   
         if (r.gt.0.) then
	   
           cosx=(x2-x1)/r
           cosy=(y2-y1)/r	
           ndots=int(r)	   
           dr=r/ndots	    
	      
           x=x1
           y=y1
           call draw_dot(nx,ny,map,x,y,intensity,thickness)
           do i=1,ndots
   	      x=x+dr*cosx
   	      y=y+dr*cosy
	      if (mod(i,20).le.10) then
     	        call draw_dot(nx,ny,map,x,y,intensity,thickness)
	      endif
 	    enddo	   	  	   
         endif	   
	   
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   subroutine draw_dot(nx,ny,map,x,y,intensity,thickness)
	   implicit none
	   
	   integer nx,ny
	   real map(nx,ny)
	   real x,y,intensity,thickness
	   
	   integer ix,iy,i,j
	   
	   ix=int(x+0.5)
	   iy=int(y+0.5)
	   
           do i=ix-thickness,ix+thickness
             do j=iy-thickness,iy+thickness
		if (i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) 
     .map(i,j)=intensity
	     enddo
	   enddo
	   	  
	   
	   return
	   end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   subroutine projection(x,y,z,theta,phi,psi,xx,yy)
	   implicit none
	
	   real x,y,z,xx,yy
	   real theta,phi,psi
	   real costheta,sintheta,cosphi,sinphi,cospsi,sinpsi

	   costheta=cos(theta)
	   sintheta=sin(theta)

	   cosphi=cos(phi)
	   sinphi=sin(phi)
	
	   cospsi=cos(psi)
	   sinpsi=sin(psi)
	  
	   xx=x*(cosphi*cospsi-sinphi*sinpsi*costheta)
     .-y*(cosphi*sinpsi*costheta+sinphi*cospsi)
     .+z*sinpsi*sintheta

	   yy=x*(cosphi*sinpsi+sinphi*cospsi*costheta)
     .+y*(cosphi*cospsi*costheta-sinphi*sinpsi)
     .-z*cospsi*sintheta

	   
	   return
	   end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This is a simple module contains subroutine that can read and write 2D fits images
! with no extensions to the primary array. It is stored in single precision.
!IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE readimage_CFHT(filename,nx,ny,npx,npy,array
     .,CRPIX,CDELT,CRVAL)
! read a 2D fits images from a fits file.
	  IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER group,firstpix,nbuffer,npixels,i,status2
      INTEGER naxes(2)
      INTEGER nx,ny,npx,npy
      double precision CRPIX(2),CDELT(2),CRVAL(2),temp(2)
      REAL array(npx,npy)
      REAL nullval,anyf
      LOGICAL anynull
      CHARACTER filename*80
      
!  The STATUS parameter must always be initialized.
      status=0
      nullval=0.
      
!  Get an unused Logical Unit Number to use to open the FITS file.
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)


!  Determine the coordinate parameters of the image.
      CALL ftgknd(unit,'CRPIX',1,2,temp,nfound,status)
!  Check that it found both CRPIX1 and CRPIX2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the CRPIXn keywords.'
          return
      end if
      CRPIX(1)=temp(1)
      CRPIX(2)=temp(2)

      CALL ftgknd(unit,'CRVAL',1,2,temp,nfound,status)
!  Check that it found both CRVAL1 and CRVAL2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the CRVALn keywords.'
          return
      end if
      CRVAL(1)=temp(1)
      CRVAL(2)=temp(2)

      CALL ftgknd(unit,'CDELT',1,2,temp,nfound,status)
!  Check that it found both CDELT1 and CDELT2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the CDELTn keywords.'
          return
      end if
      CDELT(1)=temp(1)
      CDELT(2)=temp(2)

!      CALL ftgknd(unit,'CD1_',1,2,temp,nfound,status)
!  Check that it found both CD1_1 and CD1_2 keywords.
!      if (nfound .ne. 2)then
!          print *,'READIMAGE failed to read the CD1_n keywords.'
!          return
!      end if
!      CD(1,1)=temp(1)
!      CD(1,2)=temp(2)

!      CALL ftgknd(unit,'CD2_',1,2,temp,nfound,status)
!  Check that it found both CD2_1 and CD2_2 keywords.
!      if (nfound .ne. 2)then
!          print *,'READIMAGE failed to read the CD2_n keywords.'
!          return
!      end if
!      CD(2,1)=temp(1)
!      CD(2,2)=temp(2)


!  Determine the size of the image.
      CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

!  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords.'
          return
      end if

      group=1
      nx=naxes(1)
      ny=naxes(2)

      status2=0
	  
      CALL FTG2DE(unit,group,nullval,npx,nx,ny,array,anyf,status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      END SUBROUTINE readimage_CFHT
! *************************************************************************
      SUBROUTINE writeimage(filename,nx,ny,npx,npy,array)

!  Create a FITS primary array containing a 2-D image
      IMPLICIT NONE
      
      INTEGER status,unit,blocksize,bitpix,naxis
      INTEGER i,j,group
      INTEGER nx,ny,npx,npy
      INTEGER naxes(2)
      REAL array(npx,npy)
      CHARACTER filename*80      
      LOGICAL simple,extend
      
!  The STATUS parameter must be initialized before using FITSIO.  A
!  positive value of STATUS is returned whenever a serious error occurs.
!  FITSIO uses an `inherited status' convention, which means that if a
!  subroutine is called with a positive input value of STATUS, then the
!  subroutine will exit immediately, preserving the status value. For 
!  simplicity, this program only checks the status value at the end of 
!  the program, but it is usually better practice to check the status 
!  value more frequently.

      status=0

!  Delete the file if it already exists, so we can then recreate it.
!  The deletefile subroutine is listed at the end of this file.
      CALL deletefile(filename,status)

!  Get an unused Logical Unit Number to use to open the FITS file.
!  This routine is not required;  programmers can choose any unused
!  unit number to open the file.
      CALL ftgiou(unit,status)

!  Create the new empty FITS file.  The blocksize parameter is a
!  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      CALL ftinit(unit,filename,blocksize,status)

!  Initialize parameters about the FITS image.
!  The size of the image is given by the NAXES values. 
!  The EXTEND = TRUE parameter indicates that the FITS file
!  may contain extensions following the primary array.
      simple=.true.
      bitpix=-32
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.false.

!  Write the required header keywords to the file
      CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!  Write the array to the FITS file.
!  The last letter of the subroutine name defines the datatype of the
!  array argument; in this case the 'J' indicates that the array has an
!  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
!  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
!  total number of pixels.  GROUP is seldom used parameter that should
!  almost always be set = 1.
      group=1
      CALL FTP2DE(unit,group,npx,nx,ny,array,status)

!  Write another optional keyword to the header
!  The keyword record will look like this in the FITS file:
!
!  EXPOSURE=                 1500 / Total Exposure Time
!
!     CALL ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

!  Check for any errors, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      END SUBROUTINE writeimage
! *************************************************************************
      SUBROUTINE readimage(filename,nx,ny,npx,npy,array)
! read a 2D fits images from a fits file.
	  IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER group,firstpix,nbuffer,npixels,i,status2
      INTEGER naxes(2)
      INTEGER nx,ny,npx,npy
      REAL array(npx,npy)
      REAL nullval,anyf
      LOGICAL anynull
      CHARACTER filename*80
      
!  The STATUS parameter must always be initialized.
      status=0
      nullval=0.
      
!  Get an unused Logical Unit Number to use to open the FITS file.
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)

!  Determine the size of the image.
      CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

!  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords.'
          return
      end if

      group=1
      nx=naxes(1)
      ny=naxes(2)

      status2=0
	  
      CALL FTG2DE(unit,group,nullval,npx,nx,ny,array,anyf,status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      END SUBROUTINE readimage
! *************************************************************************
      SUBROUTINE printerror(status)

!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.
      IMPLICIT NONE
      INTEGER status
      CHARACTER errtext*30
      CHARACTER errmessage*80

!  Check if status is OK (no error); if so, simply return
      IF (status .le. 0) RETURN

!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      CALL ftgerr(status,errtext)
      WRITE(*,*) 'FITSIO Error Status =',status,': ',errtext

!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      CALL ftgmsg(errmessage)
      DO WHILE (errmessage .ne. ' ')
          WRITE(*,*) errmessage
          CALL ftgmsg(errmessage)
      END DO
      END SUBROUTINE
! *************************************************************************
      SUBROUTINE deletefile(filename,status)

!  A simple little routine to delete a FITS file
      IMPLICIT NONE
      
      INTEGER status,unit,blocksize
      CHARACTER*(*) filename

!  Simply return if status is greater than zero
      IF (status .gt. 0) RETURN

!  Get an unused Logical Unit Number to use to open the FITS file
      CALL ftgiou(unit,status)

!  Try to open the file, to see if it exists
      CALL ftopen(unit,filename,1,blocksize,status)

      IF (status .eq. 0) THEN
!         file was opened;  so now delete it 
          CALL ftdelt(unit,status)
      ELSE IF (status .eq. 103) THEN
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
          CALL ftcmsg
      ELSE
!         there was some other error opening the file; delete the file anyway
          status=0
          CALL ftcmsg
          CALL ftdelt(unit,status)
      END IF

!  Free the unit number for later reuse
      CALL ftfiou(unit, status)
      END SUBROUTINE
!*****************************************************************************
      subroutine write_copyhdu(infilename,outfilename
     .,nx,ny,npx,npy,array)

C     copy the 1st and 3rd HDUs from the input file to a new FITS file

      integer status,inunit,outunit,readwrite,blocksize,morekeys,hdutype
      character infilename*80,outfilename*80

      INTEGER group
      INTEGER nx,ny,npx,npy
      REAL array(npx,npy)

 1    status=0

C     Delete the file if it already exists, so we can then recreate it
 2    call deletefile(outfilename,status)

C     Get  unused Logical Unit Numbers to use to open the FITS files
 3    call ftgiou(inunit,status)
      call ftgiou(outunit,status)

C     open the input FITS file, with readonly access
      readwrite=0
 4    call ftopen(inunit,infilename,readwrite,blocksize,status)

C     create the new empty FITS file with the standard block size
      blocksize=1
 5    call ftinit(outunit,outfilename,blocksize,status)

C     copy the primary array from the input file to the output file
      morekeys=0
 6    call ftcopy(inunit,outunit,morekeys,status)

      group=1
      CALL FTP2DE(outunit,group,npx,nx,ny,array,status)

C     append/create a new empty extension on the end of the output file
c 7    call ftcrhd(outunit,status)
C     skip to the 3rd extension in the input file
c 8    call ftmahd(inunit,3,hdutype,status)
C     copy this extension from the input file to the output file
c 9    call ftcopy(inunit,outunit,morekeys,status)  

C     close the FITS file and free the unit numbers
 10   call ftclos(inunit, status)
      call ftclos(outunit, status)
 11   call ftfiou(-1, status)

C     check for any error, and if so print out error messages
 12   if (status .gt. 0)call printerror(status)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_serial_stamps(filename,nx,ny,npx,npy,array
     .,sign_order)

C     copy a series of stamp images to a .fits file.

      integer status,unit,readwrite,blocksize,morekeys,hdutype
      save status,unit
      character filename*80

      INTEGER group,sign_order
      INTEGER nx,ny,npx,npy
      REAL array(npx,npy)

      INTEGER bitpix,naxis
      INTEGER i,j
      INTEGER naxes(2)
c      LOGICAL simple,extend
	
      character item_name*8,item_explanation*24
      integer item_value
      common /item_pass/ item_name,item_value,item_explanation

        bitpix=-32
        naxis=2
        group=1

      if (sign_order.eq.0) then
        status=0
        call deletefile(filename,status)
        call ftgiou(unit,status)
        blocksize=1
 5      call ftinit(unit,filename,blocksize,status)

      elseif (sign_order.eq.1) then
        naxes(1)=nx
        naxes(2)=ny
        CALL ftphps(unit,bitpix,naxis,naxes,status)
        CALL FTP2DE(unit,group,npx,nx,ny,array,status)

      elseif (sign_order.eq.2) then
        call ftcrhd(unit,status)
        naxes(1)=nx
        naxes(2)=ny
        CALL ftphps(unit,bitpix,naxis,naxes,status)
        CALL FTP2DE(unit,group,npx,nx,ny,array,status)

      elseif (sign_order.eq.3) then
        call ftpkyj(unit,item_name,item_value,item_explanation,status)

      elseif (sign_order.eq.-1) then
        call ftclos(unit, status)
        call ftfiou(-1, status)

      endif

C     check for any error, and if so print out error messages
 12   if (status .gt. 0)call printerror(status)
      end
