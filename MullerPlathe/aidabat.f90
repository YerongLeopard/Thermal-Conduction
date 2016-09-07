! Lennard-Jones potential, 'velocity' Verlet time integration algorithm.
! Computes kinetic and potential energy.

module Particles
!
!  Contains all the data structures containing informations about atom postion,velocity and accelertion
!
   integer, parameter :: DIM=3     !the dimension of thw system
   logical :: VelAcc = .FALSE.     ! velocities and accelerations in file or not
   integer :: N=0!the number of atoms 
   double precision, dimension(DIM) :: BoxSize! real size of the system
!     the following arrays are allocated at run time when the number of atoms N is known.
   double precision, dimension(:,:), allocatable :: pos     ! positions
   double precision, dimension(:,:), allocatable :: vel     ! velocities
   double precision, dimension(:,:), allocatable :: acc     ! accelerations
   integer, dimension(:), allocatable :: slab_label! labeling the slabs 
   double precision, dimension(:), allocatable :: ene_kin   ! kinetic energies
   double precision :: volume, density         		! these are constants 
end module Particles

module Simulation_Control
!
!  Contains the parameters controlling the simulation, which can be changed by the user.
!
   character*80:: title         !  title of this program
   character*80:: SampIn        !  name of file containing input sample
   character*80:: SampOut       !  name of file to contain output sample
   integer:: Nsteps             !  number of time steps to do
   double precision:: deltat         !  time steps (redueced units)
   double precision:: TRequested    !  desired temperature, or <0.d0 if constant E
   logical :: ConstantT     !  made .TRUE. when (TRequested >= 0)
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
   integer,parameter::SLAB=20,W=120,Mtime=7
   integer,dimension(0:SLAB-1):: number_of_atoms=0!!number of atoms in each slabs
   integer::Tsteps=0! steps of measuring the system
   double precision,dimension(0:SLAB-1):: slab_temperature_sum=0.d0,slab_temperature_SQsum=0.d0,slab_temperature=0.d0!!   
   double precision ::flux=0.d0,temperature=0.d0,temperature_sum=0.d0!
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
call Grouping!! gouping the atoms into slabs
call Initial_Printout
end subroutine Initialize


subroutine Read_Sample
!
!  Reads the initial sample from file unit 1
!
use Particles
use Simulation_Control
implicit none
double precision, dimension(DIM):: PosAtomReal,VelAtomReal,AccAtomReal!temporary variables
double precision, dimension(DIM):: Mass_center
integer :: i,k,lisin,lisout
lisin=len_trim(SampIn)
lisout=len_trim(Sampout)
open(unit=1,file=SampIn(1:lisin),status='old',action='read')!!open the file here read the input
open(unit=2,file=SampOut(1:lisout),status='replace',action='write')!!output file
read(1,'(1X,L2,I7,3E23.15)') VelAcc,N,BoxSize!VelAcc judges whether there are velocities and accelerations in file 
if ( N <= 0 ) then
   print*,'Read_Sample: FATAL: N is',N
   stop
endif
!
!  compute volume and density once for all (they do not change in the run)
!
volume  =product(BoxSize)
density = real(N) / volume
!
! allocate space here
!
allocate( pos(DIM,N) )
allocate( vel(DIM,N) )
allocate( acc(DIM,N) )
allocate( ene_kin(N) )!record the kinetic energy here
allocate( slab_label(N) )!used for labeling each atoms
!
!  read the coordinates from the file (one line per atom), normalize them to the box size along each direction and store them.Energies are set initially to zero.
!
do i=1,N
   read(1,*) PosAtomReal
   pos(:,i) = PosAtomReal/BoxSize
   ene_kin(i) = 0.d0
enddo
!
!In the case of requiring certain velocities and accelerations
!
if (VelAcc) then     ! also read velocities and accelerations
   do i=1,N
      read(1,'(1X,3E23.15)') VelAtomReal
      vel(:,i) = VelAtomReal/BoxSize
   enddo
   do i=1,N
      read(1,'(1X,3E23.15)') AccAtomReal
      acc(:,i) = AccAtomReal/BoxSize
   enddo
else                 ! set velocities and accelerations to zero
   vel = 100000000.d5
   acc = 100000000.d5
endif
!
!  translate atoms so that center of mass is at the origin
!
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
title='First Program'!name of the program
SampIn='input.dat'
SampOut='output.dat'
Nsteps=96000
deltat=6.965d-3
TRequested=-1.d0!if TRequested<0. then at constant Energy 
ConstantT = ( TRequested >= 0 )
return   !  successful exit
200 continue
   print*,'Read_Input: FATAL: premature end-of-file in standard input'
   stop
800 continue
   print*,'Read_Input: FATAL: read error in standard input'
   stop
end subroutine Read_Input


subroutine Initial_Printout
!
!  Prints informations on the run parameters on the standard output
!  Leading # are to directly use the output as a gnuplot input.
!
use Particles
use Simulation_Control
implicit none
write(2,'(2a)')'# ',title(1:len_trim(title))
write(2,'(2a)')'# Input sample: ',SampIn(1:len_trim(SampIn)),'# Output sample: ',SampOut(1:len_trim(SampOut))
if (VelAcc) then
   write(2,'(a)')'#  (positions, velocities, accelerations read from file)'
else
   write(2,'(a)')'#  (only positions read from file)'
endif
write(2,'(a,/,a,i8,/,a,es10.4,/,a,es10.4,/,a)')&
   '###########################',&
   'Number of steps:',Nsteps,&
   '# time step:',deltat, &
   '# total time:',Nsteps*deltat,&
   '###########################'
write(2,'(a,i6)')'# Number of atoms:',N
write(2,'(a,3f12.6,/,(a,f15.3),/,(a,f15.3))') &
   '# Box size:',BoxSize,&
   ' Volume:',volume,'Density:',density
if (ConstantT) then
  write(2,'(a,f12.6)')'# Constant T run with T =',TRequested
else
  write(2, '(a)')'# Fixed energy.'
endif
  
!
!  Now print headers of columns
!
print '(a,/,a,/,a)', '#', &
'#  Step       Temperature ', &    
'# ------  -------------------'
end subroutine Initial_Printout

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
integer ::step,i
double precision :: max_sq,min_sq
double precision :: chi
integer::max_index=1,min_index=1
double precision,dimension(DIM) ::temp,real_vel

!
!  We need to have the initial temperature ready in case we are going
!  at constant T:
!
call Compute_Temperature
!
!  "Velocity Verlet" integrator (see e.g. Allen and Tildesley book, p. 81).
!  Simple velocity scaling (done on velocities at the previous step)
!  applied when ConstantT is enabled.
!f
time: do step=1,Nsteps
   call Refold_Positions
   pos = pos + deltat*vel + 0.5d0*(deltat**2)*acc      ! r(t+dt)
   if (ConstantT .and. (temperature > 0.d0) ) then	  !  veloc rescale for const T
      chi = sqrt( Trequested / temperature )
      vel = chi*vel + 0.5d0*deltat*acc                 ! v(t+dt/2)
   else                          !  regular constant E dynamics
      vel = vel + 0.5d0*deltat*acc                     ! v(t+dt/2)
   endif   
   call Compute_Forces                                 ! a(t+dt)
   vel = vel + 0.5d0*deltat*acc                        ! v(t+dt)
   call Refold_Positions
   call Grouping
if(1==mod(step,W)) then
  !!!search for the maximum and minimum
  !!!the maximum in slab 0
do i=1,N
    if (0==slab_label(i)) then
      real_vel=vel(:,i)*BoxSize
      max_sq=dot_product(real_vel,real_vel)!!initialize maximum is slab 0
      max_index=i
      exit
    endif
   enddo
do i=max_index+1,N
     real_vel = BoxSize*vel(:,i)
     if (max_sq<dot_product(real_vel,real_vel) .and. 0==slab_label(i)) then 
        max_index=i
        max_sq=dot_product(real_vel,real_vel)
     endif
enddo
  !!!the maximum in slab 0

  !!!the minimum in slab SLAB/2
do i=1,N
    if (SLAB/2==slab_label(i)) then
      real_vel=vel(:,i)*BoxSize
      min_sq=dot_product(real_vel,real_vel)!!initialize minimum is slab SLAB/2
      min_index=i
      exit
    endif
   enddo
do i=min_index+1,N
     real_vel = BoxSize*vel(:,i)
     if (min_sq>dot_product(real_vel,real_vel) .and. (SLAB/2)==slab_label(i)) then 
        min_index=i
        min_sq=dot_product(real_vel,real_vel)
     endif
enddo
  !!!the minimum is slab SLAB/2
   
  !!!search for the maximum and minimum
   !!!!velocity exchange
   temp=vel(:,max_index)
   vel(:,max_index)=vel(:,min_index)
   vel(:,min_index)=temp
   !!!!velocity exchange
   !!!!calculation of flux
    flux=flux+max_sq-min_sq
    write(*,*) 0.5d0*(max_sq-min_sq)/(W*deltat)
    write(*,*) flux/(dble(step)*deltat*dble(2))
    write(*,*) (dble(step)*deltat*dble(2))
endif


if (0==mod(step,Mtime))then

   call Compute_Temperature  ! temperature at t+dt, also ene_kin
!     accumulate statistics:
   slab_temperature_sum  = slab_temperature_sum  + slab_temperature
   slab_temperature_SQsum= slab_temperature_SQsum+ slab_temperature**2! used to calculate the deviation
   Tsteps=Tsteps+1!count repeated times of measuring temperature
   write(*,700)step,temperature
   700 format (1x,i5,',',f10.5)
endif
enddo time
flux=flux/(dble(Nsteps)*deltat*dble(2))! factor 2 for kninetic energy

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
 slab_label(i)=floor((pos(DIM,i)+0.5d0)*SLAB)
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
integer::i,j,k
!
!  Reset to zero potential energies, forces, virial term
!
acc = 0.d0! Reset to zero

!
!  Loop over all pairs of particles
!
do i = 1,N-1                                     ! looping an all pairs
   do j = i+1,N
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
double precision :: ene_kin_aver
double precision, dimension(DIM) :: real_vel !temp
double precision, dimension(0:SLAB-1) :: slab_ene_kin_aver=0.d0
integer :: i
!!! intialize
temperature=0.d0
slab_temperature=0.d0
!!!!Compute the temperature for each slab 
do i=1,N
   real_vel = BoxSize*vel(:,i)                        ! real space velocity of i
   ene_kin(i) = 0.5d0*dot_product(real_vel,real_vel)  ! kin en of each atom 
   !!!!!dealing with the slab model
     slab_ene_kin_aver(slab_label(i))=slab_ene_kin_aver(slab_label(i))+ene_kin(i)
   !!!!!
enddo
do i=0,SLAB-1
 if (0==number_of_atoms(i)) then 
   slab_ene_kin_aver(i)=0.d0
   continue
 end if
 slab_ene_kin_aver(i)=slab_ene_kin_aver(i)/number_of_atoms(i)
enddo
slab_temperature=2.d0*slab_ene_kin_aver/DIM
!!!!Compute the temperature for each slab

!!!!Compute the temperature for the whole system 
ene_kin_aver = sum( ene_kin ) / real(N)
temperature = 2.d0*ene_kin_aver/DIM
!!!!Compute the temperature for the whole system 
end subroutine Compute_Temperature

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
call Print_Statistics
!
!  Write the final sample (positions, velocities, accelerations) on file
!
call Write_Sample
!
!  Deallocate all the dynamical arrays to clean memory ('ecological' practice)
!
deallocate( pos )
deallocate( vel )
deallocate( acc )
deallocate( ene_kin )
close(unit=2)
end subroutine Terminate


subroutine Print_Statistics
!
!  Print the mean value, averaged during the run, of the statistical 
!  quantities which have been accumulated
!  
use Simulation_Control
use Statistics
implicit none
integer::i=1
open(unit=9,file='slab2.dat',status='replace',action='write')
open(unit=4,file='slab.dat',status='replace',action='write')
open(unit=10,file='error2.dat',status='replace',action='write')
open(unit=7,file='error.dat',status='replace',action='write')
open(unit=8,file='flux.dat',status='replace',action='write')
if ( Nsteps <= 0 ) return               ! protection against '0 steps run'
write (8,*) flux
do i=0,SLAB/2
 write (*,140) i,slab_temperature_sum(i)/Tsteps
 write (4,141) slab_temperature_sum(i)/Tsteps

 write (7,*) sqrt(slab_temperature_SQsum(i)/real(Tsteps-1)/real(Tsteps) &
              -(slab_temperature_sum(i)**2)/(real(Tsteps)*real(Tsteps)*real(Tsteps-1)))
 
enddo
do i=SLAB/2,SLAB-1
  write (9,141) slab_temperature_sum(i)/Tsteps
  write (10,*) sqrt(slab_temperature_SQsum(i)/real(Tsteps-1)/real(Tsteps) &
              -(slab_temperature_sum(i)**2)/(real(Tsteps)*real(Tsteps)*real(Tsteps-1)))
enddo
140 format(i5,',',f8.6)
141 format(f10.6)
close(unit=4)
close(unit=7)
close(unit=8)
close(unit=9)
close(unit=10)
end subroutine Print_Statistics


subroutine Write_Sample
!
!  Write the final results into file unit 1
!
use Particles
implicit none
double precision, dimension(DIM):: PosAtomReal,VelAtomReal,AccAtomReal!temprary variables
integer::i=0
open(unit=1,file='input.dat',status='old',action='write')
open(unit=99,file='distribution.dat',status='replace',action='write')
open(unit=100,file='slab00.dat',status='replace',action='write')
open(unit=101,file='slab01.dat',status='replace',action='write')
open(unit=102,file='slab02.dat',status='replace',action='write')
open(unit=103,file='slab03.dat',status='replace',action='write')
open(unit=104,file='slab04.dat',status='replace',action='write')
open(unit=105,file='slab05.dat',status='replace',action='write')
open(unit=106,file='slab06.dat',status='replace',action='write')
open(unit=107,file='slab07.dat',status='replace',action='write')
open(unit=108,file='slab08.dat',status='replace',action='write')
open(unit=109,file='slab09.dat',status='replace',action='write')
open(unit=110,file='slab10.dat',status='replace',action='write')
open(unit=111,file='slab11.dat',status='replace',action='write')
open(unit=112,file='slab12.dat',status='replace',action='write')
open(unit=113,file='slab13.dat',status='replace',action='write')
open(unit=114,file='slab14.dat',status='replace',action='write')
open(unit=115,file='slab15.dat',status='replace',action='write')
open(unit=116,file='slab16.dat',status='replace',action='write')
open(unit=117,file='slab17.dat',status='replace',action='write')
open(unit=118,file='slab18.dat',status='replace',action='write')
open(unit=119,file='slab19.dat',status='replace',action='write')

write(1,'(1X,L2,I7,3E23.15)') .TRUE.,N,BoxSize
101 format(1X,3E23.15)
do i=1,N
   write(1,101) pos(:,i)*BoxSize
enddo
do i=1,N
   write(1,101) vel(:,i)*BoxSize
   write(99,101) vel(:,i)*BoxSize
   MAXWELL:select case (slab_label(i))
      case (0)
        write(100,101) vel(:,i)*BoxSize 
      case (1)
        write(101,101) vel(:,i)*BoxSize 
      case (2)
        write(102,101) vel(:,i)*BoxSize 
      case (3)
        write(103,101) vel(:,i)*BoxSize 
      case (4)
        write(104,101) vel(:,i)*BoxSize 
      case (5)
        write(105,101) vel(:,i)*BoxSize 
      case (6)
        write(106,101) vel(:,i)*BoxSize 
      case (7)
        write(107,101) vel(:,i)*BoxSize 
      case (8)
        write(108,101) vel(:,i)*BoxSize 
      case (9)
        write(109,101) vel(:,i)*BoxSize 
      case (10)
        write(110,101) vel(:,i)*BoxSize 
      case (11)
        write(111,101) vel(:,i)*BoxSize 
      case (12)
        write(112,101) vel(:,i)*BoxSize 
      case (13)
        write(113,101) vel(:,i)*BoxSize 
      case (14)
        write(114,101) vel(:,i)*BoxSize 
      case (15)
        write(115,101) vel(:,i)*BoxSize 
      case (16)
        write(116,101) vel(:,i)*BoxSize 
      case (17)
        write(117,101) vel(:,i)*BoxSize 
      case (18)
        write(118,101) vel(:,i)*BoxSize 
      case (19)
        write(119,101) vel(:,i)*BoxSize 
   end select MAXWELL
enddo
do i=1,N
   write(1,101) acc(:,i)*BoxSize
enddo
close(unit=1)
close(unit=99)
close(unit=100)
close(unit=101)
close(unit=102)
close(unit=103)
close(unit=104)
close(unit=105)
close(unit=106)
close(unit=107)
close(unit=108)
close(unit=109)
close(unit=110)
close(unit=111)
close(unit=112)
close(unit=113)
close(unit=114)
close(unit=115)
close(unit=116)
close(unit=117)
close(unit=118)
close(unit=119)



end subroutine Write_Sample
