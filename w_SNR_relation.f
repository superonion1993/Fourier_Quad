        subroutine w_SNR_relation(fieldname)
        implicit none
        include 'para.inc'
	
        character fieldname*6,filename*80
	integer ierror,aa,bb,cc,dd,i,j,u,v,nn	
	real a(npara)

	integer ngal
	real gal_para(ngfieldmax*nexpo,npara)

	real ssnr(nbin,4),func_weight,SNRmax,SNRmin,dSNR,snr
	external func_weight

        ngal=0

        filename='./step3/shear_'//fieldname//'.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          pause 'Shear catalog reading error!!'
        endif
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) aa,bb,cc,dd,a(1),a(2),a(3),a(4)
     .,a(5),a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14)
	  if (ierror.ge.0) then
            ngal=ngal+1 
            do j=1,15
              gal_para(ngal,j)=a(j)
            enddo
          endif
        enddo
        close(10)


	SNRmax=gal_para(1,3)
	
	SNRmin=SNRmax

	do i=2,ngal
	  if (gal_para(i,3).gt.SNRmax) SNRmax=gal_para(i,3)
	  if (gal_para(i,3).lt.SNRmin) SNRmin=gal_para(i,3)
	enddo
	
c	SNRmax=log(SNRmax)
c	SNRmin=log(SNRmin)

	SNRmax=160
	SNRmin=5

	dSNR=(SNRmax-SNRmin)/nbin
	
	do i=1,nbin
	  ssnr(i,1)=SNRmin+dSNR*(i-0.5)
	  ssnr(i,2)=0.
	  ssnr(i,3)=0.
	  ssnr(i,4)=0.
	  do j=1,ngal
c	    snr=log(gal_para(j,3))
	    snr=gal_para(j,3)
	    if (snr.ge.ssnr(i,1)-dSNR*0.5
     ..and.snr.le.ssnr(i,1)+dSNR*0.5) then
	      ssnr(i,2)=ssnr(i,2)+gal_para(j,10)**2+gal_para(j,11)**2
	      ssnr(i,3)=ssnr(i,3)+gal_para(j,12)
	      ssnr(i,4)=ssnr(i,4)+1.
	    endif
	  enddo
	  if (ssnr(i,4).gt.0) then
            ssnr(i,2)=ssnr(i,3)/ssnr(i,2)*0.1
            ssnr(i,3)=0.0002
	  endif
c	  write(*,*) ssnr(i,1),ssnr(i,2),ssnr(i,3),ssnr(i,4)
	enddo

30      format(I2.2)

	filename='./step3/w_SNR_'//fieldname//'.fits'
	call map_function(func_weight,nbin,nbin,ssnr,500,500,filename)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function func_weight(SNR)
	implicit none

	real func_weight,SNR,s2
	s2=SNR**2
	func_weight=SNR**(-2)*(s2/(s2+120.))
c	func_weight=SNR**(-2)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
