	  module MyModules
	  contains 
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Subroutine CRoswf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CRoswf(c,m,lnum,x,r1,r1d,r2,r2d,dcoef,maxd)
      
!     Purpose: To calculate the first and second kind oblate radial
!     functions, their first derivatives with respect to xi,
!     and the so-called d-coefficients, as a function of the
!     size parameter C, the order M and for a specified number LNUM 
!     of values for the degree l = M, M+1, ..., M+LNUM-1.
!
!
!
!     Input
!     c    : size parameter
!     m    : order
!     lnum : specified number LNUM of values for 
!     the degree l = M, M+1, ..., M+LNUM-1.
!     xi
!     
!     output 
!     r1 : Radial oblate spheroidal function of the first kind
!		   array of size lnum
!     r1d : derivative of R1 with respect to xi
!		   array of size lnum
!     r2 : Radial  oblate spheroidal function of the second kind
!     r2d : derivative of the Radial spheroidal function of the second kind
! 
!     dcoef: matrix of size (lnum,maxd) of the d-coefficients dcoef(i,j)=d_r(c,ml)
!     The correpondance between m,l,r and the matrix indices i,j is given by
!            l=i-1, r=2(j-1)+E(m,l)
!      or inversely
!            i=l+1 and j= (r-E(m,l))/2+1
!      with E(m,l) =O if l-m is even or E(m,l)=1 if l-m is odd  
!		maxd : maximum number of d-coefficients calculated    
		use param
		use oblfcn_mod
		
		implicit none
		
		integer,parameter :: narg=1
		
		integer m,lnum,i,j,l,maxd
		integer iopang,ioprad,iopnorm
        real(knd) c,x,ten,N0
        real(knd) arg(1)
		
        real(knd), dimension(1:lnum):: r1c,r1dc,r2c,r2dc
        real(knd), dimension(1:lnum):: r1,r1d,r2,r2d
		integer, dimension(1:lnum)::  ir1e,ir1de,ir2e,ir2de
		
        real(knd), dimension(1:lnum,1:narg):: s1c,s1dc		
        integer, dimension(1:lnum,1:narg)::	  is1e,is1de
        
        real(knd), dimension(1:lnum,1:1000):: enrl
        
        real(knd), dimension(1:lnum):: dcl
        integer, dimension(1:lnum)::   idcl
        
        real(knd), dimension(1:lnum):: dcln
        
        !real(knd), dimension(:,:), allocatable :: dcoef
        real(knd), dimension(:,:):: dcoef        

!  input data
		iopang=1 !attention : pour le calcul de normalisation  (pour coefficients d) iopang doit être = à 1
		ioprad=2 ! Radial functions and derivatives
		iopnorm=1
		arg(1)=1

		call oblfcn(c,m,lnum,ioprad,x,r1c,ir1e,r1dc,ir1de,r2c,&
     	&             ir2e,r2dc,ir2de,iopang,iopnorm,narg,arg,s1c,&
     	&             is1e,s1dc,is1de,enrl,dcl,idcl,maxd)

!       Spheroidal radial functions of first and second kind and their derivatives
    	ten=10.0e0_knd
		do l=1,lnum
		  r1(l)=r1c(l)*ten**ir1e(l)
		  r1d(l)=r1dc(l)*ten**ir1de(l)		  
		  r2(l)=r2c(l)*ten**ir2e(l)
		  r2d(l)=r2dc(l)*ten**ir2de(l)				
		end do
!       Normalisation of the first d-coefficient d_0 (c,ml) or d_1(c,ml)
		do i=0, lnum-1
			if (m==0) then
				N0=sqrt(2.0/(2*i+1))
			else
				N0=sqrt(2.0*fac(i+m)/((2*i+1)*fac(i-m)))
			end if
			dcln(i+1)=dcl(i+1)/N0
		end do
	
		!computation of all the d-coefficients by
		! recursion using enrl matrix
	    if (size(dcoef,2)< maxd) then
	    	write(*,*) '2nd dim of dcoef is too small. Number of d-coefficients is ', maxd
	    end if
		do l=1,lnum
	    	dcoef(l,1)=dcln(l)*ten**idcl(l)
	    	do j=1,maxd-1
	    		dcoef(l,j+1)= dcoef(l,j)*enrl(l,j)
	    	end do	
	    end do
	    
	    
		end subroutine CRoswf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function fac(n)  ! Factorial(n)

        use param 
        integer, intent(in):: n     
        real(knd) fac
        	fac=gamma(1.0*(n+1))
        return
        end	

		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Subroutine TransMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		subroutine TransMatrix(MC,dcoef1,dcoef2,SH,P12,LD,nr)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine calculate the transition matrix C
!
!
!	Calling parameters
!		IN	
!			dcoef1, dcoef2, SH, Ndd, P12, LD,nr 
!		OUT	
!			matrix MC			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		use param
		implicit none
		
		complex(knd), dimension(:,:), intent(out):: MC
		real(knd), dimension(:,:), intent(in):: dcoef1, dcoef2
		complex(knd), dimension(:), intent(in):: SH
		real(knd), dimension(:), intent(in) :: P12
		integer, intent(in) :: LD,nr
 		!real(knd), dimension(0:120,0:120) :: PREC

		integer :: s,r,p,l,n,i,j,ns,np,jmin,jmax,pl,pn,u,v
		
		real(knd),dimension(1:2*(nr+1)+1) :: w3j
		complex(knd) :: im, Sum, A,T,B
		real(knd), parameter :: tol=100_knd
		real(knd) :: precis

	    ns=nr
		!PREC=0
	       
	    im=complex(0.0,1.0)
	    
	    MC=0

		BL:do l=0,LD-1
			pl=mod(l,2)	! parity of l
	   		BN:do n=0,LD-1
	    		pn=mod(n,2) 	! parity of n
	    		Sum=0
	    		precis=0
	    		w3j=0
	  			BU:do u=0,ns,2
	  				B=0
	    			BV:do v=0,nr,2
	    				T=0
	    				call Wigner3j(w3j, jmin, jmax, u+pl, v+pn, 0, 0, 0)
	    				BP:do p=abs(u+pl-v-pn),u+v+pl+pn,2
	    					np=p-abs(u+pl-v-pn)+1  ! place du coefficient dans le vecteur w3j
	    					i=u/2+1
	    					j=v/2+1								
	    				    A=im**(p)*2*((2*p+1)*w3j(np)**2)*SH(p+1)*P12(p+1)&
	    				    	&*dcoef2(n+1,j)*dcoef1(l+1,i)
 	    					T=T+A  
	    				end do BP
	    				B=B+T
	    			end do BV
	    			Sum=Sum+B
	    		end do BU
	    		MC(l+1,n+1)=Sum*im**(l-n)
	    	end do BN
	    end do BL

	    
! 	    
! 	BL:	do l=0,LD-1
! 			pl=mod(l,2)	! parity of l
! 	   BN: 	do n=0,LD-1
! 	    		pn=mod(n,2) 	! parity of n
! 	    		Sum=0
! 	    		precis=0
! 	    		w3j=0
! 	  		BU:	do u=0,ns,2
! 	  				B=0
! 	    		BV:	do v=0,nr,2
! 	    				T=0
! 	    				call Wigner3j(w3j, jmin, jmax, u+pl, v+pn, 0, 0, 0)
! 	    			BP:	do p=abs(u+pl-v-pn),u+v+pl+pn,2
! 	    					np=p-abs(u+pl-v-pn)+1  ! place du coefficient dans le vecteur w3j
! 	    					i=u/2+1
! 	    					j=v/2+1								
! 	    				    A=im**(p)*2*((2*p+1)*w3j(np)**2)*SH(p+1)*P12(p+1)&
! 	    				    	&*dcoef2(n+1,j)*dcoef1(l+1,i)
!  	    					T=T+A  
! 	    				end do BP
! 	    				B=B+T
! 	    			end do BV
! 	    			Sum=Sum+B
! 	    		end do BU
! 	    		precis=abs(B/Sum)
! 	    		!PREC(l,n)=precis
! 	    		if ( precis < tol) then
! 	    			MC(l+1,n+1)=Sum*im**(l-n)
! 	    		else
! 	    			MC(l+1,n+1)=0
! 	    			write(30,*) 'Element l=',l, 'n=',n,'fails to satistfy ﻿the minimum accuracy'
! 	    			if (n.eq.0) then
! 	    				exit BL
! 	    			else
! 	    				exit BN
! 	    			end if
! 	    		end if
! 	    	end do BN
! 	    end do BL
	    


		end subroutine 	TransMatrix		
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Subroutine SPHJ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        	
        SUBROUTINE SPHJ(N,X,NM,SJ,DJ)

!      =======================================================
!      Purpose: Compute spherical Bessel functions jn(x) and
!               their derivatives
!      Input :  x --- Argument of jn(x)
!               n --- Order of jn(x)  ( n = 0,1,úúú )
!      Output:  SJ(n) --- jn(x)
!               DJ(n) --- jn'(x)
!               NM --- Highest order computed
!      Routines called:
!               MSTA1 and MSTA2 for computing the starting
!               point for backward recurrence
!      =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).EQ.1.0D-100) THEN
           DO 10 K=0,N
              SJ(K)=0.0D0
              DJ(K)=0.0D0
10		   CONTINUE
           SJ(0)=1.0D0
           DJ(1)=.3333333333333333D0
           RETURN
        ENDIF
        SJ(0)=DSIN(X)/X
        SJ(1)=(SJ(0)-DCOS(X))/X
        IF (N.GE.2) THEN
           SA=SJ(0)
           SB=SJ(1)
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) SJ(K)=F
              F0=F1
			  F1=F
15         CONTINUE
           IF (DABS(SA).GT.DABS(SB)) CS=SA/F
           IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
           DO 20 K=0,NM
              SJ(K)=CS*SJ(K)
20		   CONTINUE
        ENDIF      
        DJ(0)=(DCOS(X)-DSIN(X)/X)/X
        DO 25 K=1,NM
           DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
25		CONTINUE
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)

!      ===================================================
!      Purpose: Determine the starting point for backward  
!               recurrence such that the magnitude of    
!               Jn(x) at that point is about 10^(-MP)
!      Input :  x     --- Argument of Jn(x)
!               MP    --- Value of magnitude
!      Output:  MSTA1 --- Starting point   
!      ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
 10		CONTINUE
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)

!      ===================================================
!      Purpose: Determine the starting point for backward
!               recurrence such that all Jn(x) has MP
!               significant digits
!      Input :  x  --- Argument of Jn(x)
!               n  --- Order of Jn(x)
!               MP --- Significant digit
!      Output:  MSTA2 --- Starting point
!      ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
10		CONTINUE
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Subroutine SPHY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE SPHY(N,X,NM,SY,DY)

!       ======================================================
!       Purpose: Compute spherical Bessel functions yn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ò 0 )
!                n --- Order of yn(x) ( n = 0,1,úúú )
!       Output:  SY(n) --- yn(x)
!                DY(n) --- yn'(x)
!                NM --- Highest order computed
!       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SY(K)=-1.0D+300
              DY(K)=1.0D+300
10		CONTINUE
           RETURN
        ENDIF
        SY(0)=-DCOS(X)/X
        SY(1)=(SY(0)-DSIN(X))/X
        F0=SY(0)
        F1=SY(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X-F0
           SY(K)=F
           IF (DABS(F).GE.1.0D+300) GO TO 20              
           F0=F1
           F1=F
15		CONTINUE
20      NM=K-1
           DY(0)=(DSIN(X)+DCOS(X)/X)/X
           DO 25 K=1,NM
              DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
25		   CONTINUE
        RETURN
        END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Function hcat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function hcat( A, B ) result( X )
            use param 
    		complex(knd), dimension(:,:) :: A, B
    		complex(knd) :: X( size(A,1), size(A,2)+size(B,2) )

    		X = reshape( [ A, B], shape( X ) )
		end function
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Function vcat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		function vcat( A, B ) result( X )
		    use param 
    		complex(knd), dimension(:,:) :: A, B
    		complex(knd) :: X( size(A,1)+size(B,1), size(A,2) )

    		X = transpose( reshape( &
            [ transpose(A), transpose(B) ], &
            [ size(X,2), size(X,1) ] ) )
		end function	
 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Subroutine invertC_Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     	subroutine invertC_Matrix(R,A)
     	
     	use param
     	implicit none
     	
     	complex(knd), dimension(:,:) :: R,A
  		complex(knd),allocatable,dimension(:)::WORK
  		integer,allocatable,dimension(:)::IPIV
  		integer info,M,error

  		M=size(A,1)
  		allocate(WORK(M),IPIV(M))
		R=A

  		call ZGETRF(M,M,R,M,IPIV,info)
  		if(info .ne. 0) then
   			write(*,*)"ZGETRF failed"
  		end if

  		call ZGETRI(M,R,M,IPIV,WORK,M,info)
  		if(info .ne. 0) then
   			write(*,*)"ZGETRI failed"
  		end if


  		deallocate(IPIV,WORK,stat=error) 
  		if (error.ne.0)then
    		print *,"error:fail to release"
    	stop
  		end if 
		end subroutine invertC_matrix   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Subroutine Wigner3j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will calculate the Wigner 3j symbols
!
!		j  j2 j3
!		m1 m2 m3
!
!	for all allowable values of j. The returned values in the array j are 
!	calculated only for the limits
!
!		jmin = max(|j2-j3|, |m1|)
!		jmax = j2 + j3
!
!	To be non-zero, m1 + m2 + m3 = 0. In addition, it is assumed that all j and m are 
!	integers. Returned values have a relative error less than ~1.d-8 when j2 and j3 
!	are less than 103 (see below). In practice, this routine is probably usable up to 165.
!
!	This routine is based upon the stable non-linear recurence relations of Luscombe and 
!	Luban (1998) for the "non classical" regions near jmin and jmax. For the classical 
!	region, the standard three term recursion relationship is used (Schulten and Gordon 1975). 
!	Note that this three term recursion can be unstable and can also lead to overflows. Thus 
!	the values are rescaled by a factor "scalef" whenever the absolute value of the 3j coefficient 
!	becomes greater than unity. Also, the direction of the iteration starts from low values of j
!	to high values, but when abs(w3j(j+2)/w3j(j)) is less than one, the iteration will restart 
!	from high to low values. More efficient algorithms might be found for specific cases 
!	(for instance, when all m's are zero).
!
!	Verification: 
!
!	The results have been verified against this routine run in quadruple precision.
!	For 1.e7 acceptable random values of j2, j3, m2, and m3 between -200 and 200, the relative error
!	was calculated only for those 3j coefficients that had an absolute value greater than 
!	1.d-17 (values smaller than this are for all practical purposed zero, and can be heavily 
!	affected by machine roundoff errors or underflow). 853 combinations of parameters were found
!	to have relative errors greater than 1.d-8. Here I list the minimum value of max(j2,j3) for
!	different ranges of error, as well as the number of times this occured
!	
!	1.d-7 < error  <=1.d-8 = 103	# = 483
!	1.d-6 < error <= 1.d-7 =  116	# = 240
!	1.d-5 < error <= 1.d-6 =  165	# = 93
!	1.d-4 < error <= 1.d-5 = 167	# = 36
!
!	Many times (maybe always), the large relative errors occur when the 3j coefficient 
!	changes sign and is close to zero. (I.e., adjacent values are about 10.e7 times greater 
!	in magnitude.) Thus, if one does not need to know highly accurate values of the 3j coefficients
!	when they are almost zero (i.e., ~1.d-10) then this routine is probably usable up to about 160.
!
!	These results have also been verified for parameter values less than 100 using a code
!	based on the algorith of de Blanc (1987), which was originally coded by Olav van Genabeek, 
!	and modified by M. Fang (note that this code was run in quadruple precision, and
!	only calculates one coefficient for each call. I also have no idea if this code
!	was verified.) Maximum relative errors in this case were less than 1.d-8 for a large number
!	of values (again, only 3j coefficients greater than 1.d-17 were considered here).
!	
!	The biggest improvement that could be made in this routine is to determine when one should
!	stop iterating in the forward direction, and start iterating from high to low values. 
!
!	Calling parameters
!		IN	
!			j2, j3, m1, m2, m3 	Integer values.
!		OUT	
!			w3j			Array of length jmax - jmin + 1.
!			jmin, jmax		Minimum and maximum values
!						out output array.
!	Dependencies: None
!	
!	Written by Mark Wieczorek August (2004)
!
!	August 2009: Based on the suggestions of Roelof Rietbroek, the calculation of RS has been slightly
!	modified so that division by zero will not cause a run time crash (this behavior depends on how the 
!	compiler treats IEEE floating point exceptions). These values were never used in the original code 
!	when this did occur.
!
!	Copyright (c) 2005-2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	j2, j3, m1, m2, m3
	integer, intent(out) ::	jmin, jmax
	real*8, intent(out) ::	w3j(:)
	real*8 ::		wnmid, wpmid, scalef, denom, rs(j2+j3+1), &
				wl(j2+j3+1), wu(j2+j3+1), xjmin, yjmin, yjmax, zjmax, xj, zj
	integer :: 		j, jnum, jp, jn, k, flag1, flag2, jmid
	
	
	if (size(w3j) < j2+j3+1) then
		print*, "Error --- Wigner3j"
		print*, "W3J must be dimensioned (J2+J3+1) where J2 and J3 are ", j2, j3
		print*, "Input array is dimensioned ", size(w3j)
		stop
	endif
	
	w3j = 0.0d0
	
	flag1 = 0
	flag2 = 0
	
	scalef = 1.0d3
	
	jmin = max(abs(j2-j3), abs(m1))
	jmax = j2 + j3
	jnum = jmax - jmin + 1
	
	if (abs(m2) > j2 .or. abs(m3) > j3) then
		return
	elseif (m1 + m2 + m3 /= 0) then
		return
	elseif (jmax < jmin) then
		return
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Only one term is present
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (jnum == 1) then
		w3j(1) = 1.0d0 / sqrt(2.0d0*jmin+1.0d0)
		if ( (w3j(1) < 0.0d0 .and. (-1)**(j2-j3+m2+m3) > 0) .or. &
			(w3j(1) > 0.0d0 .and. (-1)**(j2-j3+m2+m3) < 0) ) &
			w3j(1) = -w3j(1)
		return	
	endif
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Calculate lower non-classical values for [jmin, jn]. If the second term
	!	can not be calculated because the recursion relationsips give rise to a
	!	1/0, then set flag1 to 1.  If all m's are zero, then this is not a problem 
	!	as all odd terms must be zero.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	rs = 0.0d0
	wl = 0.0d0
	
	xjmin = x(jmin)
	yjmin = y(jmin)
	
	if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then		! All m's are zero
	
		wl(jindex(jmin)) = 1.0d0
		wl(jindex(jmin+1)) = 0.0d0
		jn = jmin+1
		
	elseif (yjmin == 0.0d0) then				! The second terms is either zero
	
		if (xjmin == 0.0d0) then			! or undefined
			flag1 = 1
			jn = jmin
		else
			wl(jindex(jmin)) = 1.0d0
			wl(jindex(jmin+1)) = 0.0d0
			jn = jmin+1
		endif
		
	elseif ( xjmin * yjmin >= 0.0d0) then			! The second term is outside of the 
								! non-classical region 
		wl(jindex(jmin)) = 1.0d0
		wl(jindex(jmin+1)) = -yjmin / xjmin
		jn = jmin+1
		
	else							! Calculate terms in the non-classical region
	
		rs(jindex(jmin)) = -xjmin / yjmin
		
		jn = jmax
		do j=jmin + 1, jmax-1, 1
			denom =  y(j) + z(j)*rs(jindex(j-1))
			xj = x(j)
			if (abs(xj) > abs(denom) .or. xj * denom >= 0.0d0 .or. denom == 0.0d0) then
				jn = j-1
				exit
			else
				rs(jindex(j)) = -xj / denom
			endif
				
		enddo
		
		wl(jindex(jn)) = 1.0d0
		
		do k=1, jn - jmin, 1
			wl(jindex(jn-k)) = wl(jindex(jn-k+1)) * rs(jindex(jn-k))
		enddo
		
		if (jn == jmin) then					! Calculate at least two terms so that
			wl(jindex(jmin+1)) = -yjmin / xjmin		! these can be used in three term
			jn = jmin+1					! recursion
			
		endif

	endif
	
	if (jn == jmax) then					! All terms are calculated
	
		w3j(1:jnum) = wl(1:jnum)
		call normw3j
		call fixsign
		
		return

	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Calculate upper non-classical values for [jp, jmax].
	!	If the second last term can not be calculated because the
	!	recursion relations give a 1/0, then set flag2 to 1. 
	!	(Note, I don't think that this ever happens).
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	wu = 0.0d0
	
	yjmax = y(jmax)
	zjmax = z(jmax)
	
	if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then
	
		wu(jindex(jmax)) = 1.0d0
		wu(jindex(jmax-1)) = 0.0d0
		jp = jmax-1
		
	elseif (yjmax == 0.0d0) then
	
		if (zjmax == 0.0d0) then
			flag2 = 1
			jp = jmax
		else
			wu(jindex(jmax)) = 1.0d0
			wu(jindex(jmax-1)) = - yjmax / zjmax
			jp = jmax-1
		endif
		
	elseif (yjmax * zjmax >= 0.0d0) then
	
		wu(jindex(jmax)) = 1.0d0
		wu(jindex(jmax-1)) = - yjmax / zjmax
		jp = jmax-1

	else
		rs(jindex(jmax)) = -zjmax / yjmax

		jp = jmin
		do j=jmax-1, jn, -1
			denom = y(j) + x(j)*rs(jindex(j+1))
			zj = z(j)
			if (abs(zj) > abs(denom) .or. zj * denom >= 0.0d0 .or. denom == 0.0d0) then
				jp = j+1
				exit
			else
				rs(jindex(j)) = -zj / denom
			endif
		enddo	
		
		wu(jindex(jp)) = 1.0d0
		
		do k=1, jmax - jp, 1
			wu(jindex(jp+k)) = wu(jindex(jp+k-1))*rs(jindex(jp+k))
		enddo
		
		if (jp == jmax) then
			wu(jindex(jmax-1)) = - yjmax / zjmax
			jp = jmax-1
		endif
		
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Calculate classical terms for [jn+1, jp-1] using standard three
	! 	term rercusion relationship. Start from both jn and jp and stop at the
	! 	midpoint. If flag1 is set, then perform the recursion solely from high to
	! 	low values. If flag2 is set, then perform the recursion solely from low to high.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (flag1 == 0) then
	
		jmid = (jn + jp)/2
		
		do j=jn, jmid - 1, 1			
			wl(jindex(j+1)) = - (z(j)*wl(jindex(j-1)) +y(j)*wl(jindex(j))) / x(j)
			
			if (abs(wl(jindex(j+1))) > 1.0d0) then				! watch out for overflows.
				wl(jindex(jmin):jindex(j+1)) = wl(jindex(jmin):jindex(j+1)) / scalef
			endif
			
			if (abs(wl(jindex(j+1)) / wl(jindex(j-1))) < 1.0d0 .and. &	! if values are decreasing
				wl(jindex(j+1)) /= 0.0d0) then				! then stop upward iteration
				jmid = j+1						! and start with the downward
				exit							! iteration.
			endif
		enddo
		
		wnmid = wl(jindex(jmid))
		
		if (abs(wnmid/wl(jindex(jmid-1))) < 1.d-6 .and. &
			wl(jindex(jmid-1)) /= 0.0d0) then				! Make sure that the stopping
			wnmid = wl(jindex(jmid-1))					! midpoint value is not a zero,
			jmid = jmid - 1							! or close to it!
		endif
		
		
		do j=jp, jmid+1, -1
			wu(jindex(j-1)) = - (x(j)*wu(jindex(j+1)) + y(j)*wu(jindex(j)) ) / z(j)
			if (abs(wu(jindex(j-1))) > 1.0d0) then
				wu(jindex(j-1):jindex(jmax)) = wu(jindex(j-1):jindex(jmax)) / scalef
			endif
	
		enddo
		
		wpmid = wu(jindex(jmid))
		
		! rescale two sequences to common midpoint
		
		if (jmid == jmax) then
			w3j(1:jnum) = wl(1:jnum)
		elseif (jmid == jmin) then
			w3j(1:jnum) = wu(1:jnum)
		else
			w3j(1:jindex(jmid)) = wl(1:jindex(jmid)) * wpmid / wnmid 
			w3j(jindex(jmid+1):jindex(jmax)) = wu(jindex(jmid+1):jindex(jmax))
		endif
		
	elseif (flag1 == 1 .and. flag2 == 0) then	! iterature in downward direction only
		
		do j=jp, jmin+1, -1
			wu(jindex(j-1)) = - (x(j)*wu(jindex(j+1)) + y(j)*wu(jindex(j)) ) / z(j)
			if (abs(wu(jindex(j-1))) > 1) then
				wu(jindex(j-1):jindex(jmax)) = wu(jindex(j-1):jindex(jmax)) / scalef
			endif
		enddo
		
		w3j(1:jnum) = wu(1:jnum)
		
	elseif (flag2 == 1 .and. flag1 == 0) then	! iterature in upward direction only
		
		do j=jn, jp-1, 1
			wl(jindex(j+1)) = - (z(j)*wl(jindex(j-1)) +y(j)*wl(jindex(j))) / x(j)
			if (abs(wl(jindex(j+1))) > 1) then
				wl(jindex(jmin):jindex(j+1)) = wl(jindex(jmin):jindex(j+1))/ scalef
			endif
		enddo
		
		w3j(1:jnum) = wl(1:jnum)
		
	elseif (flag1 == 1 .and. flag2 == 1) then

		print*, "Fatal Error --- Wigner3j"
		print*, "Can not calculate function for input values, both flag1 and flag 2 are set."
		stop
	endif

	
	call normw3j
	call fixsign
		
	
	contains
	
		integer function jindex(j)
			integer :: j
			jindex = j-jmin+1
		end function jindex
	
		real*8 function a(j)
			integer :: j
			a = (dble(j)**2 - dble(j2-j3)**2) * (dble(j2+j3+1)**2 - dble(j)**2) * (dble(j)**2-dble(m1)**2)
			a = sqrt(a)
		end function a
		
		real*8 function y(j)
			integer :: j
			y = -dble(2*j+1) * &
				( dble(m1) * (dble(j2)*dble(j2+1) - dble(j3)*dble(j3+1) ) - dble(m3-m2)*dble(j)*dble(j+1) )
		end function y

		real*8 function x(j)	
			integer :: j
			x = dble(j) * a(j+1)
		end function x
		
		real*8 function z(j)
			integer :: j
			z = dble(j+1)*a(j)
		end function z
		
		subroutine normw3j
			real*8:: norm
			integer j
			
			norm = 0.0d0
			do j = jmin, jmax
				norm = norm + dble(2*j+1) * w3j(jindex(j))**2
			enddo
			
			w3j(1:jnum) = w3j(1:jnum) / sqrt(norm)
			
		end subroutine normw3j
		
		subroutine fixsign
		
			if ( (w3j(jindex(jmax)) < 0.0d0 .and. (-1)**(j2-j3+m2+m3) > 0) .or. &
				(w3j(jindex(jmax)) > 0.0d0 .and. (-1)**(j2-j3+m2+m3) < 0) ) then
				w3j(1:jnum) = -w3j(1:jnum)
			endif
			
		end subroutine fixsign

		
		end subroutine Wigner3j

	
	
		    	
      end module MyModules
      
      
