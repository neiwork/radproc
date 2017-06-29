	program N
	implicit none
	double precision :: logE(10000),logN,Ne(10000),r,t
	double precision :: a(10000,10001),suma
	integer i,j, iE, iR
	
	open(unit=10,file='eDist.txt')
	open(unit=20,file='eDist_mod.txt')
	
	iE = 100
	iR = 50
	read(10,*)
	
	do i=1,iE*iR
		read(10,*)logE(i),r,t,logN
		!if(logN.gt.0.0)then 
			Ne(i) = (10.0**logN)*(10**logE(i)*1.6e-12)**2 
		!else
		!	Ne(i)= 0.0 
		!endif
	enddo
	
	do i=1,iE
		a(i,1) = logE(i)
		do j=2,iR+1
			a(i,j) = log10(Ne((j-2)*iE+i))  !log10(Ne((i-1)*100+j))
		enddo
	enddo
	
	do i=1,iE		
		write(20,*) (a(i,j), j=1,iR+1)
		suma =0.0
	!	do j=2,iR+1
	!		suma = suma + 10**a(i,j)
	!	enddo
	!	write(30,*)logE(i),log10(suma)
	enddo


	end
	
		