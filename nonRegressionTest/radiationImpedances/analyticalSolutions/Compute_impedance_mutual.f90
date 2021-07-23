		program Compute_impedance_mutual
		
		use param
		use MyModules
		
		implicit none

		integer, parameter :: lnum=50, LD=30, m=0
		integer, parameter :: ND=100   ! nombre max de d-coefficients
		real(knd), parameter :: xi=0.0_knd
		integer i,j,k,l,maxd_1,maxd_2,NM,nf,nr
        real(knd) c1,c2,a1,a2,r12,kr12,fmin,fmax,&
        		  & pi,cmax,wave_num,csound,df	
        complex :: im,Zc,Zsc,Zm
		
        real(knd), dimension(1:lnum) :: r1_1,r1d_1,r2_1,r2d_1
        real(knd), dimension(1:lnum) :: r1_2,r1d_2,r2_2,r2d_2
        complex(knd), dimension(1:lnum) :: r4_1,r4d_1,r4_2,r4d_2
        real(knd), dimension(1:lnum,1:ND):: dcoef_1,dcoef_2
        
        real(knd), dimension(0:ND)::SJ,DJ,SY,DY
		complex(knd), dimension(1:ND+1)::SH
		real(knd), dimension(1:ND+1)::P12,P21
		
		complex(knd), dimension(1:LD,1:LD):: C12,C21
		complex(knd), dimension(1:LD):: F4_1,F4_2,D_1,D_2
		complex(knd), dimension(1:LD,1:LD):: dD_1,dD_2,Id,H11,H21,BB,CC
		complex(knd), dimension(1:2*LD,1:2*LD):: DD,MM,HH,Minv
		
		real(knd), dimension(:), allocatable :: freq
		complex(knd), dimension(:), allocatable :: Z,Zs,Zmu

		!  input data
		OPEN(UNIT=1,FILE="data.txt",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
		!frequency & calculus parameters
		read(1,*) nf, cmax, nr
		!géométrie
		read(1,*) a1, a2, r12
		! dimension in cm
 		a1=a1/100   
 		a2=a2/100
 		r12=r12/100		
		
		! Legendre polynomial for the particular argument +/- 1
		
		do i=1,ND+1
			P12(i)=1.0_knd
			P21(i)=(-1.0_knd)**(i+1)
		end do
	 
		!fréquences
		pi=acos(-1.0_knd)
		csound=340.0
		if ((pi*LD/2)<cmax) then
			write(*,*) 'pi*LD/2<cmax'
			cmax=pi*LD/2
		end if

		fmin=.1*csound/(2*pi*r12)
		fmax=cmax*csound/(2*pi*r12)
		
		allocate(freq(1:nf))	! frequency vector
		allocate(Z(1:nf))		! Impedance vector 
		allocate(Zmu(1:nf))		! Impedance (mutual) vector 
		allocate(Zs(1:nf))		! Impedance vector 1 disk 
	
		if (nf==1) then
		    df=0
		    else
			df=(fmax-fmin)/((nf-1)*1.0)
		end if

    	open(20, file='fort.20')		
		open(30, file='fort.30')	
		write(30,*) a1,a2,r12	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !boucle principale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        boucle_principale: do 100 k=1,nf
	        !write(*,*) '*'
	        freq(k)=fmin+df*(k-1)
		
			!wave_num=omega/c  wave number
			wave_num=2*pi*freq(k)/csound

			c1=wave_num*a1
			c2=wave_num*a2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!compute the spheroidal radial functions and their derivatives
			! first, second and fourth kinds
			!compute the d-coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call CRoswf(c1,m,lnum,xi,r1_1,r1d_1,r2_1,r2d_1,dcoef_1,maxd_1)			
			call CRoswf(c2,m,lnum,xi,r1_2,r1d_2,r2_2,r2d_2,dcoef_2,maxd_2)
            

			! Oblate speheroidal of the 4th kind			
			do l=1,lnum
				R4_1(l)=complex(r1_1(l),-r2_1(l))
				R4d_1(l)=complex(r1d_1(l),-r2d_1(l))
				R4_2(l)=complex(r1_2(l),-r2_2(l))
				R4d_2(l)=complex(r1d_2(l),-r2d_2(l))	
			end	do	
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			TABLE OF SPHERICAL BESSEL FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			kr12=wave_num*r12
         	call SPHJ(ND,kr12,NM,SJ,DJ)
        	call SPHY(ND,kr12,NM,SY,DY)
        	! Spherical Bessel Hankel function of the second kind
        	do i=1,ND+1
	    		SH(i)=complex(SJ(i-1),-SY(i-1))
	    	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			compute the transition matrix C, and matrix M
!			compute the inverse of M
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	    	

! matrix C		
			if ((maxd_1 < nr).or.(maxd_2<nr)) then
				write(*,*) ' warning (maxd_1 < nr).or.(maxd_2<nr)'
			end if
			call TransMatrix(C12,dcoef_1,dcoef_2,SH,P12,LD,nr) 
		    !call TransMatrix(C21,dcoef_2,dcoef_1,SH,P21,LD,nr)	
		    C21=transpose(C12)	
		      
		    
! vecteur F, matrix D, M, 	
			DD=0
			dD_1=0
			dD_2=0
			Id=0	    
		    do l=1,LD
		    	!vectors
		    	F4_1(l)=r4_1(l)*dcoef_1(l,1)
		    	F4_2(l)=r4_2(l)*dcoef_2(l,1)
		    	D_1(l)=(R1d_1(l)/R4d_1(l))   
		    	D_2(l)=(R1d_2(l)/R4d_2(l)) 
		    	!matrix
		    	dD_1(l,l)=D_1(l)
		    	dD_2(l,l)=D_2(l) 
		    	DD(l,l)=D_1(l)
			    DD(l+LD,l+LD)=D_2(l)
			    Id(l,l)=1
			end do

			BB=MATMUL(dD_1, transpose(C21))
			CC=MATMUL(dD_2, transpose(C12))

			!M=[II BB;
			!    CC II];
           	MM = vcat( hcat( Id, BB ), hcat(CC,Id ) )
           	          	
			call invertC_Matrix(Minv,MM)
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			compute the impedance Z for the interaction of the two spheroids
!			compute the impedance in free-space Zs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
			HH=MATMUL(DD, transpose(Minv));
			H11(:,:)=HH(1:LD,1:LD);
			H21(:,:)=HH(LD+1:2*LD,1:LD);

	        Zc=0
	        Zsc=0
	        Zm=0
			do i=2,LD,2
				do j=2,LD,2 
				    Zc=Zc+F4_1(i)*H11(i,j)*F4_1(j)
				    Zm=Zm+F4_2(i)*H21(i,j)*F4_1(j)
				    Zsc=Zsc+F4_1(i)*DD(i,j)*F4_1(j)
				end do
			end do

		    Z(k)=-(4.0_knd/9.0_knd)*(c1**2)*Zc
		    Zmu(k)=-(4.0_knd/9.0_knd)*(c2**2)*Zm
		    Zs(k)=-(4.0_knd/9.0_knd)*(c1**2)*Zsc
		    
			write(*,*) k,c1,freq(k),Z(k),Zmu(k)
			write(30,*) freq(k),Z(k), Zmu(k)
			write(20,*) k,c1,freq(k),realpart(Z(k)),imagpart(Z(k)),realpart(Zmu(k)),imagpart(Zmu(k)),realpart(Zs(k)),imagpart(Zs(k))
		    

100     end do boucle_principale 		
		end program Compute_impedance_mutual
		
	