! Lennard-Jones potential, 'velocity' Verlet time integration algorithm.
! Computes kinetic and potential energy.

module Particles
!
!  Contains all the data structures containing informations about atom postion,velocity and accelertion
!
   integer, parameter :: DIM = 3     !the dimension of thw system
   logical :: VelAcc = .FALSE.       ! velocities and accelerations in file or not
   integer :: N=0!the number of atoms 
   double precision, dimension(DIM) :: BoxSize! real size of the system
!     the following arrays are allocated at run time when the number of atoms N is known.
   double precision, dimension(:,:), allocatable :: pos     ! positions
   double precision, dimension(:,:), allocatable :: vel     ! velocities
   double precision, dimension(:,:), allocatable :: acc     ! accelerations
   integer, dimension(:), allocatable :: slab_label! labeling the slabs 
   double precision, dimension(:), allocatable :: vel_sq   ! kinetic energies
   double precision :: volume, density         		! these are constants 
end module Particles

module Simulation_Control
!
!  Contains the parameters controlling the simulation, which can be changed by the user.
!
   character*80:: title         !  title of this program
   character*80:: SampIn        !  name of file containing input sample
   integer:: Nsteps=0, init=0        !  number of time steps to do
   double precision:: deltat         !  time steps (redueced units)
   double precision:: TRequested    !  desired temperature, or <0.d0 if constant E
   logical :: ConstantT     !  made .TRUE. when (TRequested >= 0)
   double precision:: e=0!the amplitude of flux curve
end module Simulation_Control

module Potential
!     
!     Contains parameters of the Lennard-Jones potential
!
   double precision,parameter:: Rcutoff=3.d0   ! cutoff distance
end module Potential

module Statistics
!
!  Contains statistical quantities accumulated during the run.
!  All quantities are resetted to zero.
!
   integer,parameter::SLAB=12,W=15,Mtime=7!!W is the step interval for the impose,Mtime is the step interval for measuring
   integer,dimension(0:SLAB-1):: number_of_atoms=0!!number of atoms in each slabs
   integer::Tsteps=0! counting the number of steps measuring the system
   double precision,dimension(0:SLAB-1):: slab_temperature_sum=0.d0,slab_temperature_SQsum=0.d0,slab_temperature=0.d0!!   
   double precision ::temperature=0.d0,temperature_sum=0.d0!
end module Statistics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CODE STARTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     HERE

program Main

use Particles
implicit none
call Initialize!Initialization procedure
call Evolve_Sample! This is the main subroutine we are using
call Terminate! Deallocating space 
end program Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     ROUTINES

subroutine Initialize
!
!  Initialization procedure (called once at the beginning, before the time evolution starts)
!
use Particles
use Statistics 
implicit none
!  Read the user directives controlling the simulation
call Read_Input
!
!  Read the input sample containing the initial particle coordinates
!  (and perhaps also velocities and accelerations)
!
call Read_Sample
!
!  Print informations on the run on the standard output file
!
call Print_Parameters
call Refold_Positions
call Grouping!! grouping the atoms into slabs

end subroutine Initialize


subroutine Read_Sample
!
!  Reads the initial sample from file unit 1
!
use Statistics
use Particles
use Simulation_Control
implicit none
double precision, dimension(DIM) :: PosAtomReal,VelAtomReal,AccAtomReal!temporary variables
double precision, dimension(DIM) :: Mass_center=0.d0, displacement=0.d0
double precision,dimension(DIM,4) :: class=0.0d0 
integer :: i,j,k,l,lisin
logical :: setpos=.FALSE.
lisin=len_trim(SampIn)
open(unit=1,file=SampIn(1:lisin),status='old',action='read')!!open the file here read the input

!!! VelAcc judges whether there are velocities and accelerations in file : VelAcc = T means there is Velocity and Acceleration in input.dat
!!! N contains number of atoms
!!! Boxsize contatins height, width and depth of simulation supercell
!!! init contains the number of steps up till now
read(1,'(1X,L2,I7,3E23.15,2I7)') VelAcc,N,BoxSize,init,Tsteps
read(1,'(24E23.15)')slab_temperature_sum,slab_temperature_SQsum
if ( N <= 0 ) then
   print*,'Read_Sample: FATAL: N is',N
   stop
endif
!
!  compute volume and density once for all (they do not change in the run)
!
volume  = product(BoxSize)
density = real(N) / volume
! DEBUG
write(*,*) 'density', density
!
! allocate space here
!
allocate( pos(DIM,N) )
allocate( vel(DIM,N) )
allocate( acc(DIM,N) )
allocate( vel_sq(N) )!record the kinetic energy here
allocate( slab_label(N) )!used for labeling each atoms
!
!  read the coordinates from the file (one line per atom), normalize them to the box size along each direction and store them.Energies are set initially to zero.
!
!
!In the case of requiring certain velocities and accelerations
!
do i=1,N
   read(1,'(1X,3E23.15)') PosAtomReal
   pos(:,i) = PosAtomReal/BoxSize
   read(1,'(1X,3E23.15)') VelAtomReal
   vel(:,i) = VelAtomReal/BoxSize
   read(1,'(1X,3E23.15)') AccAtomReal
   acc(:,i) = AccAtomReal/BoxSize
enddo
  vel_sq = 0.d0

!setpos=.TRUE.

if (setpos) then
        class(1,1)=0.0d0
        class(2,1)=0.0d0
        class(3,1)=0.0d0
 
        class(1,2)=0.0d0
        class(2,2)=7.8d-1
        class(3,2)=7.8d-1
 
        class(1,3)=7.8d-1
        class(2,3)=0.0d0
        class(3,3)=7.8d-1

        class(1,4)=7.8d-1
        class(2,4)=7.8d-1
        class(3,4)=0.0d0

i=1
 do l=0,11
   do j=0,5
     do k=0,5
       displacement(1)=j*1.56
       displacement(2)=k*1.56
       displacement(3)=l*1.56
       !!class 1
       pos(:,i)=(class(:,1)+displacement)/BoxSize
       i=i+1
       !!class 2
       pos(:,i)=(class(:,2)+displacement)/BoxSize
       i=i+1
       !!class 3
       pos(:,i)=(class(:,3)+displacement)/BoxSize
       i=i+1
       !!class 4
       pos(:,i)=(class(:,4)+displacement)/BoxSize
       i=i+1
     enddo
   enddo
enddo
endif
!
!  translate atoms so that center of mass is at the origin

do k=1,DIM
   Mass_Center(k)= sum(pos(k,:))/real(N)!  compute center of mass coordinates
   pos(k,:) = pos(k,:) - Mass_Center(k)!  translate atoms so that center of mass is at the origin
enddo
!
!  all coordinates read successfully if we get to this point
!
 close(unit=1)!! close the file
return
!
!  handling of various kinds of errors
!

end subroutine Read_Sample


subroutine Read_Input
!
!  Read the parameters controlling the simulation from the standard input.
!
use Simulation_Control
implicit none
!
!  Read the input parameters controlling the simulation from the standard
!  input.
!
title='First Program' !name of the program
SampIn='input.dat'    ! input file
! Nsteps=2000*10      ! This is the number of step in this run
Nsteps = 100 ! DEBUG
deltat=6.965d-3       ! the time step
TRequested= -1.0d0    ! if TRequested<0. then at constant Energy 
! TRequested=6.6d-1 ! DEBUG
ConstantT = ( TRequested >= 0 ) ! if ConstantT = True if constant T
e = 5.866 ! the amplitude of the flux : unit
return   ! successful exit
200 continue
   print*,'Read_Input: FATAL: premature end-of-file in standard input'
   stop
800 continue
   print*,'Read_Input: FATAL: read error in standard input'
   stop
end subroutine Read_Input

subroutine Print_Parameters
  use Simulation_Control
  use Particles
  use Statistics
  open(unit=2,file='parameters.dat',status='unknown',action='write')
  !e=1.659
  !LENX=9.356
  !LENZ=62.379
  !SLAB=40
  write(2,*) e
  write(2,*) BoxSize(1)
  write(2,*) BoxSize(DIM)
  write(2,*) SLAB
  close(unit=2)
end subroutine Print_Parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  TIME EVOLUTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     ROUTINES

subroutine Evolve_Sample
!
!  This is the main subroutine, controlling the time evolution of the system.
!
use Particles
use Simulation_Control
use Statistics
implicit none
integer :: step, first, last
double precision::chi = 0.d0
!
!  We need to have the initial temperature ready in case we are going
!  at constant T:
!
call Compute_Temperature
!
!  "Velocity Verlet" integrator (see e.g. Allen and Tildesley book, p. 81).
!  Simple velocity scaling (done on velocities at the previous step)
!  applied when ConstantT is enabled.
!
call refold_positions
first=init + 1 ! This record the iteration step that has finished
last=init + Nsteps
time: do step=first, last

   pos=pos+deltat*vel+0.5d0*(deltat**2)*acc      ! r(t+dt)
   call Compute_Temperature ! T(t)
   if (ConstantT .and. (temperature>0) ) then	  !  veloc rescale for const T
      chi=sqrt( Trequested / temperature )
      vel=chi*vel+0.5d0*deltat*acc                 ! v(t+dt/2)
   else                          !  regular constant E dynamics
      vel=vel+0.5d0*deltat*acc                     ! v(t+dt/2)
   endif
   call Refold_Positions
   call Compute_Forces                                 ! a(t+dt)
   vel=vel+0.5d0*deltat*acc                        ! v(t+dt)
!   call Grouping
!   pos = pos + deltat*vel + 0.5d0*(deltat**2)*acc      ! r(t+dt)
!   call Refold_Positions
!   vel = vel + 0.5d0*deltat*acc                     ! v(t+dt/2)
!   call Compute_Forces                                 ! a(t+dt)
!   vel = vel + 0.5d0*deltat*acc                        ! v(t+dt)
!!!! Every W steps,we do the impose! 
   if(1==mod(step,W)) then
     call Do_The_Impose 
   endif

   if (1==mod(step,Mtime))then
     call Compute_Temperature  ! temperature at t+dt, also ene_kin
     write(*,700)step,'/',last
!!!! accumulate statistics:
     slab_temperature_sum  =slab_temperature_sum  + slab_temperature
     slab_temperature_SQsum=slab_temperature_SQsum+ slab_temperature**2! used to calculate the deviation
     Tsteps=Tsteps+1!count repeated times of measuring temperature
     700 format (1x,i9,a1,i9)
     call Print_Statistics
 endif  
!!!!write_input
   if (0==mod(step,200)) then
     call Write_Sample(step)
   endif
enddo time

end subroutine Evolve_Sample

subroutine Grouping
!
!This subroutine group the particles and count the number of particles in each slab
!
use Particles
use Statistics
integer :: i=1
number_of_atoms=0
do i=1,N
 slab_label(i)=floor((pos(DIM,i)+0.500001d0)*SLAB)
 if (SLAB==slab_label(i)) then 
    slab_label(i)=0
 end if
 number_of_atoms(slab_label(i))=number_of_atoms(slab_label(i))+1
enddo
end subroutine Grouping

subroutine Refold_Positions
!
!  Particles that left the box are refolded back into the box by
!  periodic boundary conditions 
!
use Particles
implicit none
pos= mod(pos,1.d0)
where ( pos >  0.5d0 ) pos = pos - 1.d0
where ( pos < -0.5d0 ) pos = pos + 1.d0
end subroutine Refold_Positions

subroutine Compute_Forces
!
!  Compute forces on atoms from the positions, using the Lennard-Jones 
!  potential.
!  Note double nested loop, giving O(N^2) time: this is a SLOW ROUTINE,
!  unsuitable for large systems.
!
use Particles
use Potential
implicit none
double precision, dimension(DIM) :: Sij,Rij
double precision :: Rsqij,phi,dphi
double precision :: rm2,rm6,rm12
integer::i,j,k,high
integer,parameter ::cut=288*2
!
!  Reset to zero potential energies, forces, virial term
!
acc = 0.d0 ! Reset to zero

!
!  Loop over all pairs of particles
!
do i = 1,N-1                                     ! looping an all pairs
  high=min(i+cut,N) 
  do j = i+1,high
      Sij = pos(:,i) - pos(:,j)                   ! distance vector between i j
      where ( abs(Sij) > 0.5d0 )                  ! (in box scaled units)
         Sij = Sij - sign(1.d0,Sij)               ! periodic boundary conditions -0.5->0.5     
      end where                                   ! applied where needed.
      do k=1,DIM 
       Rij(k)= BoxSize(k)*Sij(k)                  ! go to real space units
      end do
       Rsqij = dot_product(Rij,Rij)               ! compute square distance
      if ( Rsqij <(Rcutoff**2) ) then             ! particles are interacting?
         !  compute Lennard-Jones potenntial
         rm2 = 1.d0/Rsqij                         !  1/r^2
         rm6 = rm2**3                             !  1/r^6
         rm12 = rm6**2                            !  1/r^12
         dphi = 24.d0*rm2*(2.d0*rm12-rm6 )        !  24[2/r^14 - 1/r^8]
         acc(:,i)=acc(:,i)+dphi*Sij               ! accumulate forces
         acc(:,j)=acc(:,j)-dphi*Sij               ! (Fji = -Fij) these forces here are in scaled space
      endif
   enddo
enddo
do i = 1,cut                                     ! looping an all pairs
  do j = N-cut+i,N
      Sij = pos(:,i) - pos(:,j)                   ! distance vector between i j
      where ( abs(Sij) > 0.5d0 )                  ! (in box scaled units)
         Sij = Sij - sign(1.d0,Sij)               ! periodic boundary conditions -0.5->0.5     
      end where                                   ! applied where needed.
      do k=1,DIM 
       Rij(k)= BoxSize(k)*Sij(k)                         ! go to real space units
      end do
       Rsqij = dot_product(Rij,Rij)                ! compute square distance
      if ( Rsqij <(Rcutoff**2) ) then              ! particles are interacting?
         !  compute Lennard-Jones potenntial
         rm2 = 1.d0/Rsqij                         !  1/r^2
         rm6 = rm2**3                             !  1/r^6
         rm12 = rm6**2                            !  1/r^12
         dphi = 24.d0*rm2*( 2.d0*rm12 - rm6 )     !  24[2/r^14 - 1/r^8]
         acc(:,i) = acc(:,i) + dphi*Sij           ! accumulate forces
         acc(:,j) = acc(:,j) - dphi*Sij           !    (Fji = -Fij) these forces here are in scaled space
      endif
   enddo
enddo
end subroutine Compute_Forces

subroutine Compute_Temperature
!
!  Starting from the velocities currently stored in vel, it updates
!  the kinetic energy array ene_kin, and computes and returns the
!  averaged (on particles) kinetic energy ene_kin_aver and the
!  instantaneous temperature.
!
use Particles
use Statistics
implicit none
double precision :: vel_sq_aver
double precision, dimension(DIM) :: real_vel !temp
double precision, dimension(0:SLAB-1) :: slab_vel_sq_aver=0.d0
integer :: i
!!! intialize
temperature=0.d0
slab_temperature=0.d0
!!!!Compute the temperature for each slab 
do i=1,N
   real_vel = BoxSize*vel(:,i)                        ! real space velocity of i
   vel_sq(i) = dot_product(real_vel,real_vel)  ! kin en of each atom 
   !!!!!dealing with the slab model
     slab_vel_sq_aver(slab_label(i))=slab_vel_sq_aver(slab_label(i))+vel_sq(i)
   !!!!!
enddo
do i=0,SLAB-1
 slab_vel_sq_aver(i)=slab_vel_sq_aver(i)/number_of_atoms(i)
enddo
slab_temperature=slab_vel_sq_aver/DIM
!!!!Compute the temperature for each slab
!!!!Compute the temperature for the whole system 
vel_sq_aver = sum( vel_sq ) / real(N)
temperature =vel_sq_aver/DIM

end subroutine Compute_Temperature



subroutine Do_The_Impose
use Particles
use Simulation_Control
use Statistics
implicit none
integer ::step,i,j
real,parameter:: PI=3.1415926535898
!! there are four regions in total:
!! the cold region 1:     0---SLAB/4-2
!! the hot region 1:     SLAB/4---SLAB/2-2
!! the hot region 2:     SLAB/2-1---SLAB*3/4-2
!! the cold region 2:    SLAB*3/4---SLAB-1 
double precision,dimension(0:SLAB/4-2)::max_sqA=0.d0
double precision,dimension(SLAB/4:SLAB/2-2)::min_sqA=0.d0
double precision,dimension(SLAB/2-1:SLAB*3/4-2)::min_sqB=0.d0
double precision,dimension(SLAB*3/4:SLAB-1)::max_sqB=0.d0
double precision :: ener_diff=0.d0,diff=0.d0,X=0.d0,pro=0.d0,hot_sq=0.d0,cold_sq=0.d0,alpha=0.d0
integer,dimension(0:SLAB/4-2)::max_indexA=0
integer,dimension(SLAB/4:SLAB/2-2)::min_indexA=0
integer,dimension(SLAB/2-1:SLAB*3/4-2)::min_indexB=0
integer,dimension(SLAB*3/4:SLAB-1)::max_indexB=0
double precision,dimension(DIM) ::temp,real_vel,real_vel2
!!!search for the maximum and minimum
     !!!First do the initialization   
     !!first initialize all the index to zero to prepare,zero is also a sign that those indices is invalid
     max_indexA=0!!max_indexA are attached to slab 0,1,2,3
     min_indexA=0!!min_indexA are attached to slab 5,6,7,8
     min_indexB=0!!min_indexB are attached to slab 9,10,11,12,13
     max_indexB=0!!max_indexB are attached to slab 15,16,17,18,19
   do i=0,SLAB/4-2
     j=0
     !!!If the index is invalid,look for the first atom that is in certain slab and mark it as the atom that is attached to 	 the maximum index
     do while (0==max_indexA(i))
        j=j+1
        if (i==slab_label(j)) then
          max_indexA(i)=j
          real_vel=Boxsize*vel(:,j)
          max_sqA(i)=dot_product(real_vel,real_vel)!!max_sqA are for comparision in the following step
        endif
     enddo
   enddo
   do i=SLAB/4,SLAB/2-2
     j=0
     !!!If the index is invalid,look for the first atom that is in certain slab and mark it as the atom that is attached to 	 the minimum index
     do while (0==min_indexA(i))
        j=j+1
        if (i==slab_label(j)) then
          min_indexA(i)=j
          real_vel=Boxsize*vel(:,j)
          min_sqA(i)=dot_product(real_vel,real_vel)!!min_sqA are for comparision in the following step
        endif
     enddo
   enddo
   do i=SLAB/2-1,SLAB*3/4-2!!9---13
     j=0
     !!!If the index is invalid,look for the first atom that is in certain slab and mark it as the atom that is attached to 	 the minimum index
     do while (0==min_indexB(i))
        j=j+1
        if (i==slab_label(j)) then
          min_indexB(i)=j
          real_vel=Boxsize*vel(:,j)
          min_sqB(i)=dot_product(real_vel,real_vel)!!min_sqB are for comparision in the following step
        endif
     enddo
   enddo
   do i=SLAB*3/4,SLAB-1!!15---19
     j=0
     !!!If the index is invalid,look for the first atom that is in certain slab and mark it as the atom that is attached to 	 the maximum index
     do while (0==max_indexB(i))
        j=j+1
        if (i==slab_label(j)) then
          max_indexB(i)=j
          real_vel=Boxsize*vel(:,j)
          max_sqB(i)=dot_product(real_vel,real_vel)!!max_sqB are for comparsion in the following step
        endif
     enddo
   enddo
     !!!First do the initialization
!!!comparision
   do i=0,SLAB/4-2!!Find the fastest atom in the slab 0,1,2,3
     do j=max_indexA(i)+1,N!!We only have to consider atoms behind the first atom in the slab in question
       real_vel = BoxSize*vel(:,j)!!retrieve the velocity
       if (max_sqA(i)<dot_product(real_vel,real_vel).and.(i==slab_label(j))) then 
         max_indexA(i)=j
         max_sqA(i)=dot_product(real_vel,real_vel)
       endif
     enddo
   enddo
   do i=SLAB/4,SLAB/2-2!!Find the slowest atom in the slab 5,6,7,8
     do j=min_indexA(i)+1,N!!We only have to consider atoms behind the first atom in the slab in question
       real_vel = BoxSize*vel(:,j)!!retrieve the velocity
       if (min_sqA(i)>dot_product(real_vel,real_vel).and.(i==slab_label(j))) then 
         min_indexA(i)=j
         min_sqA(i)=dot_product(real_vel,real_vel)
       endif
     enddo
   enddo
   do i=SLAB/2-1,SLAB*3/4-2!!Find the slowest atom in the slab 9,10,11,12,13
     do j=min_indexB(i)+1,N!!We only have to consider atoms behind the first atom in the slab in question
       real_vel = BoxSize*vel(:,j)!!retrieve the velocity
       if (min_sqB(i)>dot_product(real_vel,real_vel).and.(i==slab_label(j))) then 
         min_indexB(i)=j
         min_sqB(i)=dot_product(real_vel,real_vel)
       endif
     enddo
   enddo
   do i=SLAB*3/4,SLAB-1!!Find the fastest atom in the slab 15,16,17,18,19
     do j=max_indexB(i)+1,N!!We only have to consider atoms behind the first atom in the slab in question
       real_vel = BoxSize*vel(:,j)!!retrieve the velocity
       if (max_sqB(i)<dot_product(real_vel,real_vel).and.(i==slab_label(j))) then 
         max_indexB(i)=j
         max_sqB(i)=dot_product(real_vel,real_vel)
       endif
     enddo
   enddo
   !!!search for the maximum and minimum
  ener_diff=e*W*deltat!!This is the amplitude of the flux at one imposing process
  !!!real_vel now is attached to the hotest atom, while real_vel2 is attached to the coldest atom 
  
  do i=0,SLAB/4-2!label of cold slabs:0,1,2,3
    diff=ener_diff*cos(2*PI*(i+1)/SLAB)!The energy shift in the slab i
    real_vel=Boxsize*vel(:,max_indexA(i))!real_vel is the velocity of the hotest atom in the cold slab
    real_vel2=Boxsize*vel(:,min_indexA(SLAB/2-2-i))!real-vel2 is the velocity of the coldest atom in the hot slab
    pro=dot_product(real_vel,real_vel2)
    hot_sq=dot_product(real_vel,real_vel)
    cold_sq=dot_product(real_vel2,real_vel2)
    X=(4*diff**2-2*diff*(hot_sq-cold_sq))/(hot_sq*cold_sq-pro**2)
    alpha=1-sqrt(1-X)
    temp=alpha*(real_vel-real_vel2)/2+(2*diff-alpha*(hot_sq-cold_sq)/2) & 
         *(real_vel+real_vel2)/(hot_sq+cold_sq+2*pro)! temp is the the velocity shift
    !!Do the impose here
    temp=temp/Boxsize!scale the velocity shift
    vel(:,max_indexA(i))=vel(:,max_indexA(i))-temp!!This ensure the success of the algorithm
    vel(:,min_indexA(SLAB/2-2-i))=vel(:,min_indexA(SLAB/2-2-i))+temp
  enddo
  do i=SLAB*3/4,SLAB-1!label of cold slabs:15,16,17,18,19
    diff=ener_diff*cos(2*PI*(i+1)/SLAB)!The energy shift in the slab i
    real_vel=Boxsize*vel(:,max_indexB(i))!real_vel is the velocity of the hotest atom in the cold slab
    real_vel2=Boxsize*vel(:,min_indexB(3*SLAB/2-2-i))!real-vel2 is the velocity of the coldest atom in the hot slab
    pro=dot_product(real_vel,real_vel2)
    hot_sq=dot_product(real_vel,real_vel)
    cold_sq=dot_product(real_vel2,real_vel2)
    X=(4*diff**2-2*diff*(hot_sq-cold_sq))/(hot_sq*cold_sq-pro**2)
    alpha=1-sqrt(1-X)
    temp=alpha*(real_vel-real_vel2)/2+(2*diff-alpha*(hot_sq-cold_sq)/2) & 
         *(real_vel+real_vel2)/(hot_sq+cold_sq+2*pro)! temp is the the velocity shift
    !!Do the impose here
    temp=temp/Boxsize!scale the velocity shift
    vel(:,max_indexB(i))=vel(:,max_indexB(i))-temp
    vel(:,min_indexB(SLAB*3/2-2-i))=vel(:,min_indexB(SLAB*3/2-2-i))+temp
  enddo
end subroutine Do_The_Impose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  TERMINATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    ROUTINES

subroutine Terminate
!
!  Termination procedure (called once at the end, after time evolution
!  and before program termination)
!
use Particles
implicit none
!
!  Print a line with averages, etc
!
!
!  Write the final sample (positions, velocities, accelerations) on file
!
!
!  Deallocate all the dynamical arrays to clean memory ('ecological' practice)
!

deallocate( pos )
deallocate( vel )
deallocate( acc )
deallocate( vel_sq )
end subroutine Terminate


subroutine Print_Statistics
!
!  Print the mean value, averaged during the run, of the statistical 
!  quantities which have been accumulated
!  
use Simulation_Control
use Statistics
implicit none
open(unit=4,file='slab.dat',status='unknown',action='write',position='append')
open(unit=7,file='error.dat',status='unknown',action='write',position='append')
write (4,141) slab_temperature_sum(0:SLAB-1)/Tsteps
write (7,141) sqrt(slab_temperature_SQsum(0:SLAB-1)/real(Tsteps-1)/real(Tsteps) &
              -(slab_temperature_sum(0:SLAB-1)**2)/(real(Tsteps)*real(Tsteps)*real(Tsteps-1)))
141 format(12f10.6)
close(unit=4)
close(unit=7)
end subroutine Print_Statistics


subroutine Write_Sample(step)
!
!  Write the final results into file unit 1
!
use Simulation_Control
use Statistics
use Particles
implicit none
double precision, dimension(DIM):: PosAtomReal,VelAtomReal,AccAtomReal!temprary variables
integer::i=0
integer,intent(in)::step
open(unit=1,file='input.dat',status='old',action='write')

write(1,'(1X,L2,I7,3E23.15,2I7)') .TRUE.,N,BoxSize,step,Tsteps
!slab_temperature_sum=0
!slab_temperature_SQsum=0
write(1,'(24E23.15)')slab_temperature_sum,slab_temperature_SQsum
101 format(1X,3E23.15)
do i=1,N
   write(1,101) pos(:,i)*BoxSize
   write(1,101) vel(:,i)*BoxSize
   write(1,101) acc(:,i)*BoxSize
enddo
 close(unit=1)

end subroutine Write_Sample
