c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mdinitbead  --  initialize a pimd trajectory  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mdinitbead" initializes the velocities and positions
c     for a molecular pimd trajectory, including restarts
c
c
#include "tinker_precision.h"
      subroutine mdinitbead(ibead_glob,dt,restart)
      use atmtyp
      use atoms
      use bath
      use beads
      use bound
      use couple
      use domdec
      use files
      use keys
      use freeze
      use inform
      use iounit
      use langevin
      use math
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use neigh
      use units 
      use uprior
      use usage
      use random_mod
      implicit none
      integer, intent(in) :: ibead_glob
      real(r_p), intent(in) :: dt
      logical, intent(inout) :: restart
      integer i,j,k,idyn,iglob,nh
      integer next
      integer lext,freeunit
      integer ierr
      real(r_p) :: e
      real(r_p) :: maxwell,speed
      real(r_p) :: hmax,hmass
      real(r_p) :: sum,dmass
      real(t_p) :: vec(3)
      real(r_p), allocatable :: speedvec(:)
      logical exist,heavy
      character*7 ext
      character*20 keyword
      character*120 dynfile
      character*120 record
      character*120 string
      character*3 numberbeads
     
c
c     try to restart using prior velocities and accelerations
c
      write(numberbeads, '(i3.3)') ibead_glob
      dynfile = filename(1:leng)//'_beads'//numberbeads//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      restart = exist
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
c
c     Do the domain decomposition
c
         call ddpme3d
         call reinitnl(0)
         call reassignpme(.false.)
         call mechanicstep(0)
         call allocstep
         call nblist(0)
c          call reassignpi(0)
      else
c
c     set velocities and accelerations for cartesian dynamics
c
         a(:,:) = 0.0_re_p
!$acc data present(glob,mass,v,samplevec,nbeads) async
c
#ifdef _OPENACC
         if (.not.host_rand_platform) then
            allocate(speedvec(nloc))
!$acc data create(speedvec) async
            call rand_unitgpu(samplevec(1),nloc)
            call maxwellgpu(mass,nbeads*kelvin,nloc,speedvec)
!$acc parallel loop collapse(2) async
            do i = 1, nloc
               do j = 1, 3
                  iglob = glob(i)
                  if (use(iglob)) then
                     v(j,iglob) = speedvec(i) 
     &                          * real(samplevec(3*(i-1)+j),r_p)
                  else
                     v(j,iglob)    = 0.0_re_p
                  end if
               end do
            end do
!$acc end data
            deallocate(speedvec)
         end if
#endif
         if (host_rand_platform) then
!$acc wait
!$acc update host(glob)
            do i = 1, nloc
               iglob = glob(i)
               if (use(iglob)) then
                  speed = maxwell (mass(iglob),nbeads*kelvin)
                  call ranvec (vec)
                  do j = 1, 3
                     v(j,iglob) = speed * real(vec(j),r_p)
                  end do
               else
                  do j = 1, 3
                     v(j,iglob)    = 0.0_re_p
                  end do
               end if
            end do
!$acc update device(v)
         end if
!$acc end data
      end if

c
c     check for any prior dynamics coordinate sets
c
      i = 0
      exist = .true.
      do while (exist)
         i = i + 1
         lext = 3
         call numeral (i,ext,lext)
         dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                        '.'//ext(1:lext)
         inquire (file=dynfile,exist=exist)
         if (.not.exist .and. i.lt.100) then
            lext = 2
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                        '.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
         if (.not.exist .and. i.lt.10) then
            lext = 1
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                        '.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
      end do
      nprior = i - 1
      return
      end
c
