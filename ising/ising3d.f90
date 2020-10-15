program ising3d

implicit none

!declare variables-----------------------
!parameters
real*8:: J_ising=1.0, KT=4.7, MU=0.33, B=0

!lattice
integer:: L, N, iter, x1, x2, y1, y2, z1, z2
integer, dimension(:,:,:), allocatable:: lat

!thermodynamic quanities
real*8:: E=0,M=0,Ei,Ef,dE,r,E1=0,E2=0,Cv,M1=0,M2=0,ksi

!loop variables
integer:: i, j, k, time, ii, jj, kk


!initialisation--------------------------
print*,"enter lattice side size: "
read*, L
N = L*L*L
!allocate array
allocate (lat(L,L,L))
!initialise lattice
do i=1,L
   do j=1,L
      do k=1,L
         lat(i,j,k)=1 !all spins UP
      end do
   end do
end do

!initial energy/magnetisation
do i=1,L
  do j=1,L
    do k=1,L
        !six-nearest neighbours
        x1=i-1; x2=i+1; y1=j-1; y2=j+1; z1=k-1; z2=k+1
        !periodic boundary conditions
        if (i==1) x1=L
        if (i==L) x2=1
        if (j==1) y1=L
        if (j==L) y2=1
        if (k==1) z1=L
        if (k==L) z2=1
      
        E=E-J_ising*float(lat(i,j,k)*(lat(x1,j,k)+lat(x2,j,k)+lat(i,y1,k)+lat(i,y2,k)+lat(i,j,z1)+lat(i,j,z2)))-MU*B*lat(i,j,k)
        M=M+float(lat(i,j,k)) 
    end do
  end do
end do

E=E*0.5
print*, "initial lattice energy: ",E
print*, "initial lattice magnetisation: ",M

!Metropolis loop------------------------------------
print*, "enter number of MCS iterations: "
read*, iter

!write to file
open(15,file='ising3d_Cv_data.csv') !file for Temperature data
write(15,'(*(g0.7,:,","))')"Temp","E","E/N","Cv","ksi"
do while(KT>=3.8)
  do time=1,iter !MCS loop
    do ii=1,L
      do jj=1,L
        do kk=1,L
          !random lattice point
          call RANDOM_NUMBER(r); i=int(r*float(L))+1
          call RANDOM_NUMBER(r); j=int(r*float(L))+1            
          call RANDOM_NUMBER(r); k=int(r*float(L))+1
          !six-nearest neighbours
          x1=i-1; x2=i+1; y1=j-1; y2=j+1; z1=k-1; z2=k+1
          !periodic boundary conditions
          if (i==1) x1=L
          if (i==L) x2=1
          if (j==1) y1=L
          if (j==L) y2=1
          if (k==1) z1=L
          if (k==L) z2=1
            
          !before spin-flip
          Ei=-J_ising*float(lat(i,j,k)*(lat(x1,j,k)+lat(x2,j,k)+lat(i,y1,k)+lat(i,y2,k)+lat(i,j,z1)+lat(i,j,z2)))-MU*B*lat(i,j,k)
          !spin-flip
          lat(i,j,k)=-lat(i,j,k)
          !after spin-flip
          Ef=-J_ising*float(lat(i,j,k)*(lat(x1,j,k)+lat(x2,j,k)+lat(i,y1,k)+lat(i,y2,k)+lat(i,j,z1)+lat(i,j,z2)))-MU*B*lat(i,j,k)
          !energy difference
          dE=Ef-Ei
            
          if(dE<=0)then   !accept spin-flip
            E=E+dE
            M=M+2*float(lat(i,j,k))
          else
            call RANDOM_NUMBER(r)
            if(r<exp(-dE/KT))then !accept spin-flip with relative probability
              E=E+dE
              M=M+2*float(lat(i,j,k))
            else
              lat(i,j,k)=-lat(i,j,k)
            end if
          end if
        end do
      end do
    end do
    
    !after equilibrium
    if(time>10000)then
      E1=E1+E
      E2=E2+(E*E)
      M1=M1+M
      M2=M2+(M*M)
    end if     
  end do
  
  E1=E1/float(time-10000)
  E2=E2/float(time-10000)
  M1=M1/float(time-10000)
  M2=M2/float(time-10000)
  !Specific Heat    
  Cv=(E2-(E1*E1))/(KT*KT)
  !Susceptibility
  ksi=(M2-(M1*M1))/(KT) 
  !write second data
  write(15,'(*(g0.7,:,","))')KT,E1,E1/float(N),M1,M1/float(N),Cv,ksi
  !change temperature
  KT=KT-0.02
  
end do   

close(15)                  
end program ising3d

      
      
      
        
    
    
    
    
    
    
    
    
    
    
    
    
    
