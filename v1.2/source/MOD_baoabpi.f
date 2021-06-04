c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################################
c     ##                                                                                   ##
c     ##  module baoabopi  --  BAOAB Path Integral Langevin molecular dynamics step    ##
c     ##                                                                                   ##
c     #######################################################################################
c
c
c     "baoabpi" performs a single molecular dynamics time step
c     via the BAOAB recursion formula
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c

      module baoabpi
      implicit none

      contains

c       subroutine baoabpi1 (istep,dt)
c       use atmtyp
c       use atoms
c       use bath
c       use cutoff
c       use domdec
c       use energi
c       use freeze
c       use langevin
c       use mdstuf
c       use moldyn
c       use timestat
c       use units
c       use usage
c       use mpi
c       implicit none
c       integer, intent(in) :: istep
c       real*8, intent(in) :: dt
c       integer i,j,iglob,ibead
c       real*8 dt_x,factor
c       real*8 etot,eksum,epot
c       real*8 temp,pres
c       real*8 part1,part2
c       real*8 a1,a2,normal
c       real*8 ekin(3,3)
c       real*8 stress(3,3)
c       real*8 time0,time1
c c
c c
c c     find quarter step velocities via BAOAB recursion
c c
c       do i = 1, nloc
c          iglob = glob(i)
c          if (use(iglob)) then
c             do j = 1, 3
c                v(j,iglob) = v(j,iglob) + 0.5*dt*a(j,iglob)
c             end do
c          end if
c       end do

c       end subroutine baoabpi1

      subroutine apply_b_pi (istep,dt)
      use atmtyp
      use atoms
      use bath
      use beads
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer, intent(in) :: istep
      real*8, intent(in) :: dt
      integer i,j,iglob,ibead,ierr
      real*8 dt_x,factor
      real*8 pres
      real*8 part1,part2
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)

      allocate (derivs(3,nbloc))
      derivs = 0d0
c
c     get the potential energy and atomic forces
c
      if(contract) then
        call gradfast(eintrapi_loc,derivs)
      else
        call gradient (epotpi_loc,derivs)
      endif
c
c     MPI : get total energy (and allreduce the virial)
c
      if (contract) then
        call reduceen(eintrapi_loc)
        call MPI_BCAST(eintrapi_loc,1,MPI_REAL8,0,COMM_TINKER,ierr)
        call commforcesrespa(derivs,.true.)
      else
        call reduceen(epotpi_loc)
        call MPI_BCAST(epotpi_loc,1,MPI_REAL8,0,COMM_TINKER,ierr)
        call commforces(derivs)
      endif
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
c      write(*,*) 'x 1 = ',x(1),y(1),v(1,1),a(1,1)
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = v(j,iglob) + dt*a(j,iglob)
               !v(j,iglob) = v(j,iglob) + 0.5*dt*a(j,iglob)
            end do
         end if
      end do

      deallocate(derivs)

      end subroutine apply_b_pi

c
c     !subroutine commposbead: communicate all positions to global master to be able
c     to make normal mode transformation and propagation,
c     then communicate back the positions after NM propagation
c
      subroutine apply_aoa_pi(istep,dt,pos_full,vel_full)
      use atoms
      use boxes
      use domdec
      use mpi
      use units
      use beads
      use langevin 
      use bath
      use cutoff
      use energi
      use freeze
      use mdstuf
      use moldyn
      use timestat
      use units
      use usage
      use atmtyp
      use math
      use iounit
      use inform
      use molcul
      implicit none
      real*8, intent(in) :: dt
      integer, intent(in) :: istep
      real*8, intent(inout),allocatable::pos_full(:,:,:),vel_full(:,:,:)
      ! pos_full, vel_full must be filled only for ranktot 0
      integer i,j,k,iproc, ibead,l
      integer :: modstep
      real*8 :: eigx0,eigv0
      real*8 :: dt2, a1,a2, sqrtmass
      real*8 :: mpitime1, mpitime2
      real*8 :: time_tot, time_com

      dt2=0.5*dt

      if(isobaric) extvolold = extvol

      !time_tot=0d0
      !time_com=0d0
      !mpitime1=mpi_wtime()

      ! gather beads info at ranktot 0
      !call gather_polymer(pos_full,vel_full,forces_full)

      !call compute_observables_pi(pos_full,vel_full,forces_full,istep)

      !if(ranktot.eq.0) call mdstatpi(istep,dt)
      
      !mpitime2=mpi_wtime()
      !time_com=mpitime2-mpitime1

c     propagate AOA
      if(ranktot .eq. 0) then

c       transform to normal modes
        DO i=1,n ; DO j=1,3
          vel_full(j,i,:)=matmul(eigmattr,vel_full(j,i,:))
          pos_full(j,i,:)=matmul(eigmattr,pos_full(j,i,:))
        ENDDO ; ENDDO

       if(isobaric) then
c         propagate volume velocity from previous step
          aextvol = 3.d0*nbeads*convert*(
     &        extvol*(presvir-atmsph)/prescon 
     &        +gasconst*kelvin 
     &      )/masspiston
          vextvol = vextvol + dt2*aextvol          
        endif

        call apply_a_pi(pos_full,vel_full,dt2)

        call apply_o_pi(pos_full,vel_full,dt)

        call apply_a_pi(pos_full,vel_full,dt2)

c         update Ekcentroid  and presvir     
        presvir=presvir-prescon*(Ekcentroid/(3*nbeads*volbox))
        Ekcentroid=0
        DO i=1,n ; DO j=1,3          
          Ekcentroid=Ekcentroid+mass(i)*vel_full(j,i,1)**2
        ENDDO ; ENDDO
        Ekcentroid=Ekcentroid/convert        
        presvir=presvir+prescon*(Ekcentroid/(3*nbeads*volbox))

        if(isobaric) then
c         propagate volume velocity
          aextvol = 3.d0*nbeads*convert*(
     &        extvol*(presvir-atmsph)/prescon 
     &        +gasconst*kelvin 
     &      )/masspiston
          vextvol = vextvol + dt2*aextvol
        endif        

c         transform back to coordinates
        DO i=1,n ; DO j=1,3 
          vel_full(j,i,:)=matmul(eigmat,vel_full(j,i,:))
          pos_full(j,i,:)=matmul(eigmat,pos_full(j,i,:))          
        ENDDO ; ENDDO 
      endif

      !mpitime1=mpi_wtime()
      !time_tot=mpitime1-mpitime2

      !call broadcast_polymer(pos_full,vel_full)
      !if(ranktot.eq.0) deallocate(pos_full,vel_full,forces_full)

      if(isobaric) call rescale_box_pi(istep)     
     
      
      !mpitime2=mpi_wtime()
      !time_com=time_com+mpitime2-mpitime1

c      if (ranktot.eq.0) then
c            write(*,*) 'Time spent in com', time_com
c            write(*,*) 'Time spent in propagation', time_tot
c      endif
      end subroutine apply_aoa_pi


       subroutine apply_a_pi(eigpos,eigvel,tau)
      use atoms
      use units
      use beads
      use bath
      implicit none
      real*8, intent(inout), allocatable :: eigpos(:,:,:),eigvel(:,:,:)
      real*8, intent(in) :: tau
      real*8 :: eigx0,eigv0
      real*8 :: acentroid,scale
      integer :: i,j,k

        if(isobaric) then
c         propagate centroid isobaric (half step)
          acentroid = sinh(tau*vextvol)/vextvol
          scale=exp(tau*vextvol)            
          extvol = extvol*exp(3.d0*tau*vextvol)

          DO i=1,n  ; DO j=1,3
            eigpos(j,i,1)=eigpos(j,i,1)*scale
     &        +acentroid*eigvel(j,i,1)
            eigvel(j,i,1)=eigvel(j,i,1)/scale
          ENDDO ; ENDDO 
        else
c         propagate centroid (half step)
          DO i=1,n ; DO j=1,3
            eigpos(j,i,1)=eigpos(j,i,1) + tau*eigvel(j,i,1)
          ENDDO ; ENDDO
        endif

c         propagate springs (half step)
        DO k=2,nbeads
          DO i=1,n ; DO j=1,3           
            eigx0=eigpos(j,i,k)
            eigv0=eigvel(j,i,k)
            eigpos(j,i,k)=eigx0*cos(omkpi(k)*tau)
     $            +eigv0*sin(omkpi(k)*tau)/omkpi(k)
            eigvel(j,i,k)=eigv0*cos(omkpi(k)*tau)
     $           -eigx0*sin(omkpi(k)*tau)*omkpi(k)
          ENDDO ; ENDDO
        ENDDO

      end subroutine apply_a_pi

      subroutine apply_o_pi(eigpos,eigvel,tau)
        use atmtyp
        use atoms
        use units
        use beads
        use bath
        use langevin
        use qtb
        implicit none
        real*8, intent(inout),allocatable :: eigpos(:,:,:),eigvel(:,:,:)
        real*8, intent(in) :: tau
        real*8 :: acentroid,a1p,a2p,gammak,a1,a2
        real*8 :: normal
        integer :: i,j,k

c         propagate langevin (full step, TRPMD)
        DO k=1,nbeads
          gammak=max(gamma,omkpi(k))
          a1 = exp(-gammak*tau)
          a2 = sqrt((1.-a1**2)*nbeads*boltzmann*kelvin)
          DO i=1,n ; DO j=1,3 
            eigvel(j,i,k)=eigvel(j,i,k)*a1 
     &              + a2*normal()/sqrt(mass(i))
          ENDDO; ENDDO
        ENDDO

        if(isobaric) then
c         langevin piston (full step)
          a1p = exp(-gammapiston*tau)
          a2p = sqrt((1-a1p**2)*nbeads*boltzmann*kelvin/masspiston)
          vextvol = vextvol*a1p + a2p*normal()
        endif

      end subroutine apply_o_pi

      subroutine compute_grad_slow()
      use atmtyp
      use atoms
      use bath
      use beads
      use cutoff
      use deriv
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use neigh
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer i,j,iglob,ibead,ierr
      integer iloc,iloc1,iloc2,iloc3,ilocrec1,ilocrec2,ilocrec3
      real*8 dt_x,factor
      real*8 pres
      real*8 part1,part2
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)

      allocate(derivs(3,nbloc))
      derivs=0d0
      
      
      call gradslow(einterpi_loc,derivs)      
      call reduceen(einterpi_loc)
      call MPI_BCAST(einterpi_loc,1,MPI_REAL8,0,COMM_TINKER,ierr)

      call commforcesrespa(derivs,.false.)
      !write(0,*) einterpi_loc

c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
c      write(*,*) 'x 1 = ',x(1),y(1),v(1,1),a(1,1)
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
        do j = 1, 3
            a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
c               v(j,iglob) = v(j,iglob) + dt*a(j,iglob)
c               v(j,iglob) = v(j,iglob) + 0.5*dt*a(j,iglob)
          end do
        end if
      end do
      deallocate(derivs)
      
      end subroutine compute_grad_slow



      subroutine apply_b_slow(polymer,dt)
      use atmtyp
      use atoms
      use bath
      use beads
      use cutoff
      use deriv
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use neigh
      use timestat
      use units
      use usage
      use mpi
      implicit none
      type(POLYMER_COMM_TYPE) :: polymer,polymer_ctr
      real*8, intent(in) :: dt
      integer i,j,k,iglob,ibead,ierr
      integer iloc,iloc1,iloc2,iloc3,ilocrec1,ilocrec2,ilocrec3
      real*8 dt_x,factor
      real*8 pres
      real*8 part1,part2
      real*8 time0,time1

      if(ranktot.eq.0) then
        DO i=1,n;  DO j=1,3
          polymer%vel(j,i,:)=polymer%vel(j,i,:)+dt*convert
     &                            *polymer%forces_slow(j,i,:)/mass(i)
          enddo;enddo
      endif

      end subroutine apply_b_slow

      end module baoabpi

