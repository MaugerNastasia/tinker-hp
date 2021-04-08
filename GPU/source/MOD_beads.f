c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module beads   --  pimd variables                              ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     nbeads  number of global replicas used for PIMD simulations
c     nbeadsloc  number of local replicas used for PIMD simulations
c     nproctot number of process for beads parallelism
c     nproc number of process for gradient parallelism
c     ncomm number of communicators for beads parallelism
c     locbead : array to switch from global to local beads
c     ibeadsloc : number of the current local beads
c     
c  dynamical variables replicated for PIMD simulations
c
c     pospi : array of positions
c     velpi : array of velocities
c     api : array of accelerations
c     
c     array used for domain decomposition within a bead:
c      glob
c      loc
c      repart
c      repartrec
c      domlen
c      domlenrec number of reciprocal atoms in the reciprocal domains
c      domlenpole
c      domlenpolerec
c      globrec
c      locrec
c      globrec1
c      locrec1
c      bufbegrec
c      bufbegpole
c      bufbeg
c     buflen1,buflen2,buf1,buf2,bufbeg1,bufbeg2 explicit direct-reciprocal atomic correspondance, 
c      nloc
c      nbloc
c      nlocrec
c      nblocrec local + neighbors reciprocal number of atoms
c      nlocnl local nl number of atoms
c      nblocrecdir local + neighbors direct+reciprocal number of atoms
c
c      molecule 
c      nmolelocpi,molculeglobpi
c
c      VDW :
c        nvdwblocpi,vdwglobpi,nvdwlocnlpi,vdwglobnlpi
c        nvlstpi,vlstpi
c
c      BONDS:
c        nbondlocpi,bndglobpi
c
c      STRETCH-BEND:
c        nstrbndlocpi,strbndglobpi
c
c      ANGLE-ANGLE:
c        nanganglocpi,angangglobpi
c
c      OP-BENDING:
c        nopbendlocpi,opbendglobpi
c
c      OP-DIST:
c        nopdistlocpi,opdistglobpi
c
c      IMPROP:
c        niproplocpi,impropglobpi
c
c      IMPTOR:
c        nitorslocpi,imptorglobpi 
c
c      TORSION:
c        ntorslocpi,torsglobpi 
c
c      PITORSION:
c        npitorslocpi,pitorsglobpi 
c
c      STRETCH-TORSION:
c        nstrtorlocpi,strtorglobpi 
c
c      TORSION-TORSION:
c        ntortorlocpi,tortorglobpi 
c
c      ANGLE:
c        nanglelocpi,angleglobpi 
c
c      CHARGE:
c        nionreclocpi,chgrecglobpi 
c        nionlocpi,chgglobpi
c        nionlocnlpi,chgglobnlpi
c        nelstpi,elstpi
c
c      MULTIPOLE:
c        npolereclocpi,polerecglobpi 
c        npolelocpi,poleglobpi,polelocpi
c        npolelocnlpi,poleglobnlpi
c        
c      POLARIZATION:
c        udaltpi,upaltpi,nualtpi
c        uindpi,uinppi
c        
c
c      STATS:
c        etot_sumpi,etot2_sumpi,eint_sumpi,eint2_sumpi
c        epot_sumpi,epot2_sumpi,ekin_sumpi,ekin2_sumpi
c        temp_sumpi,temp2_sumpi,pres_sumpi,pres2_sumpi
c        dens_sumpi,dens2_sumpi
c
c      TIME:
c        timestep
c
#include "tinker_precision.h"
      module beads
      implicit none
      integer :: ibead_loaded_loc,ibead_loaded_glob
      integer :: nbeads,nbeadsloc,nprocbeads   
!$acc declare create(nbeads)


      integer :: rank_beadloc, ncomm  

      logical :: contract
      integer :: nbeads_ctr, nbeadsloc_ctr  

      real(r_p) :: temppi,temppi_cl
      real(r_p) :: epotpi_loc,etotpi_loc
      real(r_p) :: eksumpi_loc,ekinpi_loc(3,3)
      real(r_p) :: ekprim,ekvir,presvir
      real(r_p) :: ekcentroid
      real(r_p) :: ekprim_ave, epotpi_ave, temppi_ave      
      real(r_p) :: eintrapi_loc
      real(r_p) :: einterpi_loc
      
      real(8), allocatable :: eigmat(:,:)
      real(8), allocatable :: eigmattr(:,:)
      real(8), allocatable ::  omkpi(:)
      
      real(8), allocatable :: contractor_mat(:,:)
      real(8), allocatable :: uncontractor_mat(:,:)

      TYPE POLYMER_COMM_TYPE
        integer :: nbeads
        integer :: nbeadsloc_max
        integer, allocatable :: nbeadscomm(:),nloccomm(:,:) 
        integer, allocatable :: globbeadcomm(:,:)
        integer, allocatable :: repart_dof_beads(:,:,:)

        real(8), allocatable :: pos(:,:,:),vel(:,:,:)
        real(8), allocatable :: forces(:,:,:)
        real(8), allocatable :: forces_slow(:,:,:)
      END TYPE     

      TYPE BEAD_TYPE
        real(r_p), allocatable :: x(:), y(:),z(:)
        real(r_p), allocatable :: v(:,:), a(:,:)
        real(r_p) :: dedv
        real(r_p) :: epot, etot, eksum, ekin(3,3)
        real(r_p) :: eintra
        real(r_p) :: einter

        integer :: ibead_loc
        integer :: ibead_glob
        integer :: contraction_level

        !PARALLELISM
        integer :: nloc,nbloc,nlocrec
        integer :: nblocrec,nlocnl
        integer :: nblocrecdir
        integer, allocatable :: glob(:),loc(:),repart(:)
        integer, allocatable :: ineignl(:) !,locnl(:)
        integer, allocatable :: repartrec(:), domlen(:)
        integer, allocatable :: domlenrec(:),domlenpole(:)
        integer, allocatable :: domlenpolerec(:),globrec(:)
        integer, allocatable :: locrec(:),globrec1(:)
        integer, allocatable :: locrec1(:)
        integer, allocatable :: bufbegrec(:)
        integer, allocatable :: bufbegpole(:),bufbeg(:)
        integer, allocatable :: buflen1(:),buflen2(:)
        integer, allocatable :: buf1(:),buf2(:)
        integer, allocatable :: bufbeg1(:),bufbeg2(:)
        integer :: nmoleloc
        integer, allocatable :: molculeglob(:)

        !
        !     BOND-STRETCHING
        !
        integer :: nbondloc
        integer, allocatable :: bndglob(:)
        !
        !     STRETCH-BENDING
        !
        integer :: nstrbndloc
        integer, allocatable :: strbndglob(:)
        !
        !     UREY-BRADLEY
        !
        integer :: nureyloc
        integer, allocatable :: ureyglob(:)
        !
        !     ANGLE-ANGLE
        !
        integer :: nangangloc
        integer, allocatable :: angangglob(:)
        !
        !     OP-BENDING
        !
        integer :: nopbendloc
        integer, allocatable :: opbendglob(:)
        !
        !     OP-DIST
        !
        integer :: nopdistloc
        integer, allocatable :: opdistglob(:)
        !
        !     IMPROP
        !
        integer :: niproploc
        integer, allocatable :: impropglob(:)
        !
        !     IMPTOR
        !
        integer :: nitorsloc
        integer, allocatable :: imptorglob(:)
        !
        !     TORSION
        !
        integer :: ntorsloc
        integer, allocatable ::torsglob(:)
        !
        !     PITORSION
        !
        integer :: npitorsloc
        integer, allocatable :: pitorsglob(:)
        !
        !     STRETCH-TORSION
        !
        integer :: nstrtorloc
        integer, allocatable ::strtorglob(:)
        !
        !     TORSION-TORSION
        !
        integer :: ntortorloc
        integer, allocatable :: tortorglob(:)
        !
        !     ANGLE
        !
        integer :: nangleloc
        integer, allocatable :: angleglob(:)
        !
        !     CHARGE
        !
        integer :: nionrecloc, nionloc, nionlocnl
        integer, allocatable :: chgrecglob(:)
        integer, allocatable ::chgglob(:)
        integer, allocatable :: chgglobnl(:)
        integer, allocatable :: nelst(:),elst(:,:)
        !
        !     MULTIPOLE
        !
        integer :: npolerecloc,npoleloc,npolebloc, npolelocnl
        integer, allocatable :: polerecglob(:)
        integer, allocatable :: poleglob(:)
        integer, allocatable :: poleloc(:)
        integer, allocatable :: poleglobnl(:)
        !
        !      POLARIZATION
        !
        integer :: nualt
        real(t_p), allocatable :: udalt(:,:,:),upalt(:,:,:)
        real(t_p), allocatable :: uind(:,:),uinp(:,:)
        !
        !     VDW
        !
        integer :: nvdwbloc,nvdwlocnl
        integer, allocatable :: vdwglob(:)
        integer, allocatable :: vdwglobnl(:)
        integer, allocatable :: nvlst(:),vlst(:,:)
        !
        !     STAT
        !
        real(t_p) :: etot_sum,etot2_sum
        real(t_p) :: eint_sum,eint2_sum
        real(t_p) :: epot_sum,epot2_sum
        real(t_p) :: ekin_sum,ekin2_sum
        real(t_p) :: temp_sum,temp2_sum
        real(t_p) :: pres_sum,pres2_sum
        real(t_p) :: dens_sum,dens2_sum
          
        !
        !     TIME
        !
        real(t_p) :: timestep
      END TYPE

      TYPE(BEAD_TYPE), allocatable :: beadsloc(:)
      TYPE(BEAD_TYPE), allocatable :: beadsloc_ctr(:)


      save

      contains

      subroutine deallocate_polymer(polymer)
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer

        if(allocated(polymer%pos)) deallocate(polymer%pos)
        if(allocated(polymer%vel)) deallocate(polymer%vel)
        if(allocated(polymer%forces)) deallocate(polymer%forces)

        if(allocated(polymer%nbeadscomm)) deallocate(polymer%nbeadscomm)
        if(allocated(polymer%nloccomm)) deallocate(polymer%nloccomm)
        if(allocated(polymer%globbeadcomm)) 
     &                      deallocate(polymer%globbeadcomm)
        if(allocated(polymer%repart_dof_beads)) 
     &                      deallocate(polymer%repart_dof_beads)

      end subroutine

      subroutine allocpi()
      use angle
      use atoms
      use bath
      use bitor
      use domdec
      use molcul
      use neigh
      use sizes
      use pitors
      use potent
      use tors
      use uprior
      use units
      use mdstuf
      implicit none
      integer ibead,i,ierr,j,k
      real(8), allocatable :: WORK(:)
      real(8), allocatable :: WORK_ctr(:),eigMat_ctr(:,:),omkpi_ctr(:)

      ibead_loaded_loc = -1
      ibead_loaded_glob = -1

      if(ranktot.eq.0) then

        Ekcentroid=nfree*nbeads*kelvin*gasconst      
      
        if (allocated(eigmat)) deallocate(eigmat)
        allocate(eigmat(nbeads,nbeads))
        if (allocated(eigmattr)) deallocate(eigmattr)
        allocate(eigmattr(nbeads,nbeads))
        if (allocated(omkpi)) deallocate(omkpi)
        allocate(omkpi(nbeads))        
        allocate(WORK(3*nbeads))

        eigmat=0
        DO i=1,nbeads-1
          eigmat(i,i)=2
          eigmat(i+1,i)=-1
          eigmat(i,i+1)=-1
        ENDDO
        eigmat(1,nbeads)=-1
        eigmat(nbeads,1)=-1
        eigmat(nbeads,nbeads)=2
        call DSYEV('V','U',nbeads,eigMat,nbeads, 
     $       omkpi,WORK,3*nbeads,ierr)
        omkpi(1)=0
        omkpi(:)=sqrt(omkpi)*(nbeads*boltzmann*kelvin/hbar)
        eigmattr=transpose(eigmat)
        deallocate(WORK)


        if (contract) then
          allocate(eigmat_ctr(nbeads_ctr,nbeads_ctr))
          allocate(omkpi_ctr(nbeads_ctr))
          allocate(WORK_ctr(3*nbeads_ctr))
          eigmat_ctr=0
          do i=1,nbeads_ctr-1
            eigmat_ctr(i,i)=2
            eigmat_ctr(i+1,i)=-1
            eigmat_ctr(i,i+1)=1
          enddo
          eigmat_ctr(1,nbeads_ctr)=-1
          eigmat_ctr(nbeads_ctr,1)=-1
          eigmat_ctr(nbeads_ctr,nbeads_ctr)=2
          call DSYEV('V','U',nbeads_ctr,eigMat_ctr,nbeads_ctr, 
     $        omkpi_ctr,WORK_ctr,3*nbeads_ctr,ierr)

          if (allocated(contractor_mat)) deallocate(contractor_mat)
          allocate(contractor_mat(nbeads_ctr,nbeads))
          contractor_mat(:,:) = 0._8
          do i=1,nbeads_ctr ; do j=1,nbeads            
            do k=1,nbeads_ctr
              contractor_mat(i,j) = contractor_mat(i,j)
     &          + eigmat_ctr(i,k)*eigmat(j,k)
            enddo
          enddo; enddo
          contractor_mat=contractor_mat*sqrt(nbeads_ctr*1._8/nbeads)

          if (allocated(uncontractor_mat)) deallocate(uncontractor_mat)
          allocate(uncontractor_mat(nbeads,nbeads_ctr))
          uncontractor_mat = nbeads*transpose(contractor_mat)/nbeads_ctr
          deallocate(WORK_ctr,eigMat_ctr,omkpi_ctr)
        endif
      endif

      if(allocated(beadsloc)) deallocate(beadsloc)
      allocate(beadsloc(nbeadsloc))
      do i=1,nbeadsloc
        call allocbead(beadsloc(i))
        beadsloc(i)%ibead_loc = i
        beadsloc(i)%ibead_glob = rank_beadloc*int(nbeads/ncomm) + i
        ! OK if the last rank_beadloc takes the remaining beads
        !write(0,*) ranktot,"loc=",i,"glob=",beadsloc(i)%ibead_glob
      enddo

      if(contract) then
        if(allocated(beadsloc_ctr)) deallocate(beadsloc_ctr)
        allocate(beadsloc_ctr(nbeadsloc_ctr))
        do i=1,nbeadsloc_ctr
          call allocbead(beadsloc_ctr(i))
          beadsloc_ctr(i)%ibead_loc = i
          beadsloc_ctr(i)%ibead_glob = 
     &       rank_beadloc*int(nbeads_ctr/ncomm) + i
          ! OK if the last rank_beadloc takes the remaining beads
        enddo
      endif

      end subroutine allocpi

      subroutine allocbead(bead)
      use angle
      use atoms
      use bath
      use bitor
      use domdec
      use molcul
      use neigh
      use sizes
      use pitors
      use potent
      use tors
      use uprior
      use units
      use mdstuf
      implicit none
      TYPE(BEAD_TYPE), intent(inout) :: bead

      call deallocate_bead(bead)

      allocate(bead%x(n))
      allocate(bead%y(n))
      allocate(bead%z(n))
      allocate(bead%v(3,n))
      allocate(bead%a(3,n))
    
    
      allocate(bead%glob(n))
      allocate(bead%loc(n))
      allocate(bead%repart(n))
      allocate(bead%repartrec(n))
      allocate(bead%domlen(nproc))
      allocate(bead%domlenrec(nproc))
      allocate(bead%domlenpole(nproc))
      allocate(bead%domlenpolerec(nproc))
      allocate(bead%globrec(n))
      allocate(bead%locrec(n))
      allocate(bead%globrec1(n))
      allocate(bead%locrec1(n))
      allocate(bead%bufbegrec(nproc))
      allocate(bead%bufbegpole(nproc))
      allocate(bead%bufbeg(nproc))
      allocate(bead%buflen1(nproc))
      allocate(bead%buflen2(nproc))
      allocate(bead%buf1(n))
      allocate(bead%buf2(n))
      allocate(bead%bufbeg1(nproc))
      allocate(bead%bufbeg2(nproc))
      allocate(bead%molculeglob(nmol))
      allocate(bead%ineignl(n))
      !allocate(bead%locnl(n))      


      if (use_vdw) then
        allocate(bead%vdwglob(n))
        allocate(bead%vdwglobnl(n))
      end if

      if (use_bond) then
        allocate(bead%bndglob(8*n))
      end if

      if (use_strbnd) then
        allocate(bead%strbndglob(nangle))
      end if
      
      if (use_urey) then
        allocate(bead%ureyglob(nangle))
      end if

      if (use_angang) then
        allocate(bead%angangglob(ntors))
      end if

      if (use_opbend) then
        allocate(bead%opbendglob(nangle))
      end if

      if (use_opdist) then
        allocate(bead%opdistglob(n))
      end if

      if (use_improp) then
        allocate(bead%impropglob(6*n))
      end if

      if (use_imptor) then
        allocate(bead%imptorglob(6*n))
      end if

      if (use_tors) then
        allocate(bead%torsglob(6*n))
      end if

      if (use_pitors) then 
        allocate(bead%pitorsglob(ntors))
      end if

      if (use_strtor) then
        allocate(bead%strtorglob(ntors))
      end if

      if (use_tortor) then
        allocate(bead%tortorglob(nbitor))
      end if

      if (use_angle) then
        allocate(bead%angleglob(4*n))
      end if

      if (use_charge) then
        allocate(bead%chgrecglob(n))
        allocate(bead%chgglob(n))
        allocate(bead%chgglobnl(n))
      end if

      if (use_polar) then
        allocate(bead%udalt(maxualt,3,n))
        allocate(bead%upalt(maxualt,3,n))
        allocate(bead%uind(3,n))
        allocate(bead%uinp(3,n))
      end if

      if (use_mpole) then
        allocate(bead%polerecglob(n))
        allocate(bead%poleglob(n))
        allocate(bead%poleloc(n))
        allocate(bead%poleglobnl(n))
      end if


      end subroutine

      subroutine deallocate_bead(bead)
        implicit none
        TYPE(BEAD_TYPE), intent(inout) :: bead

        if (allocated(bead%x)) deallocate(bead%x)
        if (allocated(bead%y)) deallocate(bead%y)
        if (allocated(bead%z)) deallocate(bead%z)
        if (allocated(bead%v)) deallocate(bead%v)
        if (allocated(bead%a)) deallocate(bead%a)   
        if (allocated(bead%glob)) deallocate(bead%glob)
        if (allocated(bead%loc)) deallocate(bead%loc)
        if (allocated(bead%repart)) deallocate(bead%repart)
        if (allocated(bead%repartrec)) deallocate(bead%repartrec)
        if (allocated(bead%domlen)) deallocate(bead%domlen)
        if (allocated(bead%domlenrec)) deallocate(bead%domlenrec)
        if (allocated(bead%domlenpole)) deallocate(bead%domlenpole)
        if (allocated(bead%domlenpolerec))deallocate(bead%domlenpolerec)
        if (allocated(bead%globrec)) deallocate(bead%globrec)
        if (allocated(bead%locrec)) deallocate(bead%locrec)
        if (allocated(bead%globrec1)) deallocate(bead%globrec1)
        if (allocated(bead%locrec1)) deallocate(bead%locrec1)
        if (allocated(bead%bufbegrec)) deallocate(bead%bufbegrec)
        if (allocated(bead%bufbegpole)) deallocate(bead%bufbegpole)
        if (allocated(bead%bufbeg)) deallocate(bead%bufbeg)
        if (allocated(bead%buflen1)) deallocate(bead%buflen1)
        if (allocated(bead%buflen2)) deallocate(bead%buflen2)
        if (allocated(bead%bufbeg1)) deallocate(bead%bufbeg1)
        if (allocated(bead%buf1)) deallocate(bead%buf1)
        if (allocated(bead%buf2)) deallocate(bead%buf2)
        if (allocated(bead%bufbeg2)) deallocate(bead%bufbeg2)
        if (allocated(bead%molculeglob)) deallocate(bead%molculeglob)
        if (allocated(bead%ineignl)) deallocate(bead%ineignl)
        !if (allocated(bead%locnl)) deallocate(bead%locnl)
        if (allocated(bead%vdwglob)) deallocate(bead%vdwglob)
        if (allocated(bead%vdwglobnl)) deallocate(bead%vdwglobnl)
        if (allocated(bead%bndglob)) deallocate(bead%bndglob)
        if (allocated(bead%strbndglob)) deallocate(bead%strbndglob)
        if (allocated(bead%ureyglob)) deallocate(bead%ureyglob)
        if (allocated(bead%angangglob)) deallocate(bead%angangglob)
        if (allocated(bead%opbendglob)) deallocate(bead%opbendglob)
        if (allocated(bead%opdistglob)) deallocate(bead%opdistglob)
        if (allocated(bead%impropglob)) deallocate(bead%impropglob)
        if (allocated(bead%imptorglob)) deallocate(bead%imptorglob)
        if (allocated(bead%torsglob)) deallocate(bead%torsglob)
        if (allocated(bead%pitorsglob)) deallocate(bead%pitorsglob)
        if (allocated(bead%strtorglob)) deallocate(bead%strtorglob)
        if (allocated(bead%tortorglob)) deallocate(bead%tortorglob)
        if (allocated(bead%angleglob)) deallocate(bead%angleglob)
        if (allocated(bead%chgrecglob)) deallocate(bead%chgrecglob)
        if (allocated(bead%chgglob)) deallocate(bead%chgglob)
        if (allocated(bead%chgglobnl)) deallocate(bead%chgglobnl)
        if (allocated(bead%udalt)) deallocate(bead%udalt)
        if (allocated(bead%upalt)) deallocate(bead%upalt)
        if (allocated(bead%uind)) deallocate(bead%uind)
        if (allocated(bead%uinp)) deallocate(bead%uinp)
        if (allocated(bead%polerecglob)) deallocate(bead%polerecglob)
        if (allocated(bead%poleglob)) deallocate(bead%poleglob)
        if (allocated(bead%poleloc)) deallocate(bead%poleloc)
        if (allocated(bead%poleglobnl)) deallocate(bead%poleglobnl)

      end subroutine deallocate_bead
      

      subroutine initbead(istep,bead,skip_parameters,updateGPU)
      use angang
      use angle
      use atoms
      use atmlst
      use atmtyp
      use bond
      use charge
      use domdec
      use energi
      use improp
      use imptor
      use moldyn
      use mpole
      use neigh
      use opbend, only: nopbendloc
      use opdist
      use pitors
      use polar
      use potent
      use sizes
      use stat
      use strbnd
      use strtor
      use timestat
      use tors
      use tortor
      use uprior
      use urey
      use vdw
      use virial
      implicit none
      integer, intent(in) :: istep
      type(BEAD_TYPE), intent(inout) :: bead 
      logical,intent(in) :: skip_parameters,updateGPU
      integer ibead,i,modnl

c
c     positions, speed, mass
c
      if(updateGPU) then
!$acc update host(x(:),y(:),z(:),v(:,:),a(:,:))
      endif

      bead%x = x
      bead%y = y
      bead%z = z
      bead%v = v
      bead%a = a

!$acc update host(eksumpi_loc,ekinpi_loc) async
      if(contract) then
!$acc update host(eintrapi_loc,einterpi_loc) async
!$acc wait
        bead%eintra=eintrapi_loc
        bead%einter=einterpi_loc
        bead%epot=eintrapi_loc+einterpi_loc
        bead%etot=eintrapi_loc+einterpi_loc+eksumpi_loc
      else
!$acc update host(epotpi_loc) async
!$acc wait
        bead%epot=epotpi_loc
        bead%etot=epotpi_loc+eksumpi_loc
      endif

      bead%eksum=eksumpi_loc
      bead%ekin=ekinpi_loc
      bead%dedv = dedv

      if(updateGPU) then
!$acc update host(glob(:))
      endif
      bead%nloc = nloc
      bead%glob = glob
c
c     STAT
c
      bead%etot_sum = etot_sum
      bead%etot2_sum= etot2_sum
      bead%eint_sum = eint_sum
      bead%eint2_sum= eint2_sum
      bead%epot_sum = epot_sum
      bead%epot2_sum= epot2_sum
      bead%ekin_sum = ekin_sum
      bead%ekin2_sum= ekin2_sum
      bead%temp_sum = temp_sum
      bead%temp2_sum= temp2_sum
      bead%pres_sum = pres_sum
      bead%pres2_sum= pres2_sum
      bead%dens_sum = dens_sum
      bead%dens2_sum= dens2_sum

      if(skip_parameters) return

      modnl = mod(istep,ineigup)
c
c     parallelism
c
      bead%nbloc = nbloc
      bead%nlocrec = nlocrec
      bead%nblocrec = nblocrec
      bead%nlocnl = nlocnl
      bead%nblocrecdir = nblocrecdir

      if(updateGPU) then
!$acc update host(loc(:),ineignl(:)) async
!$acc update host(repart(:),repartrec(:)) async
!$acc update host(domlenrec(:),domlenpole(:),domlenpolerec(:)) async
!$acc update host(globrec(:),locrec(:),globrec1(:),locrec1(:)) async
!$acc update host(bufbegpole(:)) async
!$acc update host(buflen2(:),buf1(:),buf2(:)) async
!$acc update host(bufbeg1(:),bufbeg2(:)) async
!$acc wait
      endif
      bead%loc = loc
      bead%ineignl = ineignl
      !bead%locnl = locnl
      bead%repart = repart
      bead%repartrec = repartrec
      bead%domlen = domlen
      bead%domlenrec = domlenrec
      bead%domlenpole = domlenpole
      bead%domlenpolerec = domlenpolerec
      bead%globrec = globrec
      bead%locrec = locrec
      bead%globrec1 = globrec1
      bead%locrec1 = locrec1
      bead%bufbegrec = bufbegrec
      bead%bufbegpole = bufbegpole
      bead%bufbeg = bufbeg
      bead%buflen1 = buflen1
      bead%buf1(1:nblocrecdir) = buf1
      bead%buflen2 = buflen2
      bead%buf2(1:nblocrecdir) = buf2
      bead%bufbeg1 =  bufbeg1
      bead%bufbeg2 = bufbeg2

c
c     VDW 
c
      if (use_vdw) then
        if(updateGPU) then
!$acc update host(nvdwbloc,nvdwlocnl,vdwglob(:),vdwglobnl(:))
        endif
        bead%nvdwbloc = nvdwbloc
        bead%vdwglob = vdwglob
        bead%vdwglobnl = vdwglobnl
        bead%nvdwlocnl = nvdwlocnl
      end if
c
c     BONDS
c
      if (use_bond) then
        if(updateGPU) then
!$acc update host(nbondloc,bndglob(:))
        endif
        bead%bndglob = bndglob
        bead%nbondloc = nbondloc
      end if
c
c     STRETCH-BEND
c
      if (use_strbnd) then
        if(updateGPU) then
!$acc update host(nstrbndloc,strbndglob(:))
        endif
        bead%strbndglob = strbndglob
        bead%nstrbndloc = nstrbndloc
      end if
c
c     UREY-BRADLEY
c
      if (use_urey) then
        if(updateGPU) then
!$acc update host(nureyloc,ureyglob(:))
        endif
        bead%ureyglob = ureyglob
        bead%nureyloc = nureyloc
      end if
c
c     ANGlE-ANGLE
c
      if (use_angang) then
        if(updateGPU) then
!$acc update host(nangangloc,angangglob(:))
        endif
        bead%angangglob = angangglob
        bead%nangangloc = nangangloc
      end if
c
c     OP-BENDING
c
      if (use_opbend) then
        if(updateGPU) then
!$acc update host(nopbendloc,opbendglob(:))
        endif
        bead%opbendglob = opbendglob
        bead%nopbendloc = nopbendloc
      end if
c
c     OP-DIST
c
      if (use_opdist) then
        if(updateGPU) then
!$acc update host(nopdistloc,opdistglob(:))
        endif
        bead%opdistglob = opdistglob
        bead%nopdistloc = nopdistloc
      end if
c
c     IMPROP
c
      if (use_improp) then
        if(updateGPU) then
!$acc update host(niproploc,impropglob(:))
        endif
        bead%impropglob = impropglob
        bead%niproploc = niproploc
      end if
c
c     IMPTOR
c
      if (use_imptor) then
        if(updateGPU) then
!$acc update host(nitorsloc,imptorglob(:))
        endif
        bead%imptorglob = imptorglob
        bead%nitorsloc = nitorsloc
      end if
c
c     TORSION
c
      if (use_tors) then
        if(updateGPU) then
!$acc update host(ntorsloc,torsglob(:))
        endif
        bead%torsglob = torsglob
        bead%ntorsloc = ntorsloc
      end if
c
c     PITORSION
c
      if (use_pitors) then
        if(updateGPU) then
!$acc update host(npitorsloc,pitorsglob(:))
        endif
        bead%pitorsglob = pitorsglob
        bead%npitorsloc = npitorsloc
      end if
c
c     STRETCH-TORSION
c
      if (use_strtor) then
        if(updateGPU) then
!$acc update host(nstrtorloc,strtorglob(:))
        endif
        bead%strtorglob = strtorglob
        bead%nstrtorloc = nstrtorloc
      end if
c
c     TORSION-TORSION
c
      if (use_tortor) then
        if(updateGPU) then
!$acc update host(ntortorloc,tortorglob(:))
        endif
        bead%tortorglob = tortorglob
        bead%ntortorloc = ntortorloc
      end if
c
c     ANGLE
c
      if (use_angle) then
        if(updateGPU) then
!$acc update host(nangleloc,angleglob(:))
        endif
        bead%angleglob = angleglob
        bead%nangleloc = nangleloc
      end if
c
c     CHARGE
c
      if (use_charge) then
        if(updateGPU) then
!$acc update host(nionloc,nionlocnl,nionrecloc) async
!$acc update host(chgrecglob(:),chgglob(:),chgglobnl(:)) async
!$acc wait
        endif
        bead%chgrecglob = chgrecglob
        bead%nionrecloc = nionrecloc
        bead%chgglob = chgglob
        bead%nionloc = nionloc
        bead%chgglobnl = chgglobnl
        bead%nionlocnl = nionlocnl
      end if
c      if ((use_charge.or.use_mpole).and.(modnl.eq.0)) then
c        nelstpi(1:nlocnl,ibead) = nelst
c        elstpi(:,1:nlocnl,ibead) = elst
c      end if
c
c     MULTIPOLE
c
      if (use_mpole) then
        if(updateGPU) then
!$acc update host(npoleloc,npolebloc,poleglobnl) async
!$acc update host(npolelocnl,npolerecloc) async
!$acc update host(polerecglob(:),poleglob(:),poleloc(:)) async
!$acc wait
        endif
        bead%polerecglob = polerecglob
        bead%npolerecloc = npolerecloc
        bead%poleglob = poleglob
        bead%poleloc = poleloc
        bead%npoleloc = npoleloc
        bead%npolebloc = npolebloc
        bead%poleglobnl = poleglobnl
        bead%npolelocnl = npolelocnl
      end if
c
c     POLARIZATION
c
      if (use_polar) then
        if(updateGPU) then
!$acc update host(udalt(:,:,:),upalt(:,:,:)) async
!$acc update host(uind(:,:),uinp(:,:)) async
!$acc wait
        endif
        bead%nualt = nualt
        bead%udalt = udalt
        bead%upalt = upalt
        bead%uind = uind
        bead%uinp = uinp
      end if

c
c     TIME
c
      bead%timestep = timestep

      end subroutine initbead

      subroutine resize_nl_arrays_bead(istep,bead)
      use angle
      use atoms
      use bath
      use bitor
      use domdec
      use molcul
      use neigh
      use sizes
      use pitors
      use potent
      use tors
      use uprior
      use units
      use mdstuf
      implicit none
      integer, intent(in) :: istep
      TYPE(BEAD_TYPE), intent(inout) :: bead
      integer nblocrecdirmax,modnl
  
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return    

      !nlocnlmax = maxval(nlocnlpi)
      if (use_vdw) then
        if (allocated(bead%nvlst)) deallocate(bead%nvlst)
        allocate(bead%nvlst(bead%nlocnl))
        if (allocated(bead%vlst)) deallocate(bead%vlst)
        allocate(bead%vlst(maxvlst,bead%nlocnl))
      end if
      if (use_mpole.or.use_charge) then
        if (allocated(bead%nelst)) deallocate(bead%nelst)
        allocate(bead%nelst(bead%nlocnl))
        if (allocated(bead%elst)) deallocate(bead%elst)
        allocate(bead%elst(maxelst,bead%nlocnl))
      end if

      end subroutine resize_nl_arrays_bead
c
      subroutine savebeadnl(istep,bead,updateGPU)
      use domdec
      use neigh
      use potent
      implicit none
      TYPE(BEAD_TYPE), intent(inout) :: bead
      integer, intent(in) ::  istep
      LOGICAL, intent(in) :: updateGPU
      integer modnl

      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return

      if (use_vdw) then
        if(updateGPU) then
!$acc update host(nvlst(:),vlst(:,:))
        endif
        bead%nvlst = nvlst
        bead%vlst = vlst
      end if
      
      if ((use_charge).or.(use_mpole)) then
        if(updateGPU) then
!$acc update host(nelst(:),elst(:,:))
        endif
        bead%nelst = nelst
        bead%elst = elst
      end if
      return
      end
      

      subroutine pushbead(istep,bead,skip_parameters,updateGPU)
      use angang
      use angle
      use atoms
      use atmlst
      use atmtyp
      use bond
      use charge
      use domdec
      use energi
      use improp
      use imptor
      use moldyn
      use mpole
      use neigh
      use opbend, only: nopbendloc
      use opdist
      use pitors
      use polar
      use potent
      use sizes
      use stat
      use strbnd
      use strtor
      use timestat
      use tors
      use tortor
      use uprior
      use urey
      use vdw
      use virial
      implicit none
      integer, intent(in) :: istep
      type(BEAD_TYPE), intent(inout) :: bead
      LOGICAL, intent(in) :: skip_parameters,updateGPU
      integer modnl

      ibead_loaded_glob = bead%ibead_glob
      ibead_loaded_loc = bead%ibead_loc

      !write(0,*) "push energies"

      if(contract) then
        eintrapi_loc=bead%eintra
        einterpi_loc=bead%einter
        epotpi_loc=bead%eintra+bead%einter
        eksumpi_loc=bead%eksum
        ekinpi_loc=bead%ekin
        dedv=bead%dedv
!$acc update device(eintrapi_loc,einterpi_loc) async
      else
        epotpi_loc=bead%epot
        eksumpi_loc=bead%eksum
        ekinpi_loc=bead%ekin
        dedv=bead%dedv
      endif
!$acc update device(eksumpi_loc,ekinpi_loc,epotpi_loc) async

      !write(0,*) "push position"      
      x = bead%x
      y = bead%y
      z = bead%z
      v = bead%v
      a= bead%a
      if(updateGPU) then
!$acc update device(x(:),y(:),z(:),v(:,:),a(:,:)) async
      endif

      nloc = bead%nloc
      glob = bead%glob
      if(updateGPU) then
!$acc update device(glob(:)) async
      endif
c
c     STAT
c
      ! write(0,*) "push stat"
      etot_sum  = bead%etot_sum
      etot2_sum = bead%etot2_sum
      eint_sum  = bead%eint_sum
      eint2_sum = bead%eint2_sum
      epot_sum  = bead%epot_sum
      epot2_sum = bead%epot2_sum
      ekin_sum  = bead%ekin_sum
      ekin2_sum = bead%ekin2_sum
      temp_sum  = bead%temp_sum
      temp2_sum = bead%temp2_sum
      pres_sum  = bead%pres_sum
      pres2_sum = bead%pres2_sum
      dens_sum  = bead%dens_sum
      dens2_sum = bead%dens2_sum
c
      if (skip_parameters) return

      modnl = mod(istep,ineigup)
c
c     parallelism
c
      ! write(0,*) "push parallelism"
      nbloc = bead%nbloc
      nlocrec = bead%nlocrec
      nblocrec = bead%nblocrec
      nlocnl = bead%nlocnl
      nblocrecdir = bead%nblocrecdir

      loc = bead%loc
      ineignl = bead%ineignl
      ! locnl = bead%locnl
      repart = bead%repart
      repartrec = bead%repartrec
      domlen = bead%domlen
      domlenrec = bead%domlenrec
      domlenpole = bead%domlenpole
      domlenpolerec = bead%domlenpolerec
      globrec = bead%globrec
      locrec = bead%locrec
      globrec1 = bead%globrec1
      locrec1 = bead%locrec1
      bufbegrec = bead%bufbegrec
      bufbegpole = bead%bufbegpole
      bufbeg = bead%bufbeg
      buflen1 = bead%buflen1
      buf1 = bead%buf1
      buflen2 = bead%buflen2
      buf2 = bead%buf2
      bufbeg1 = bead%bufbeg1
      bufbeg2 = bead%bufbeg2

      if(updateGPU) then
!$acc update device(loc(:),ineignl(:)) async
!$acc update device(repart(:),repartrec(:)) async
!$acc update device(domlenrec(:),domlenpole(:),domlenpolerec(:)) async
!$acc update device(globrec(:),locrec(:),globrec1(:),locrec1(:)) async
!$acc update device(bufbegpole(:)) async
!$acc update device(buflen2(:),buf1(:),buf2(:)) async
!$acc update device(bufbeg1(:),bufbeg2(:)) async
      endif

c
c     VDW
c
      if (use_vdw) then
        ! write(0,*) "push VDW"
        nvdwbloc = bead%nvdwbloc
        vdwglob = bead%vdwglob
        vdwglobnl = bead%vdwglobnl
        nvdwlocnl = bead%nvdwlocnl
        nvlst = bead%nvlst
        vlst = bead%vlst
        if(updateGPU) then
!$acc update device(nvdwbloc,nvdwlocnl,vdwglob(:),vdwglobnl(:)) async
!$acc update device(nvlst(:),vlst(:,:)) async
        endif
      endif
c
c     BOND
c
      if (use_bond) then
        !write(0,*) "push BOND"
        nbondloc = bead%nbondloc
        bndglob = bead%bndglob
        if(updateGPU) then
!$acc update device(nbondloc,bndglob(:)) async
        endif
      endif
c
c     STRETCH-BEND
c
      if (use_strbnd) then
        !write(0,*) "push STRETCH-BEND"
        nstrbndloc = bead%nstrbndloc
        strbndglob = bead%strbndglob
        if(updateGPU) then
!$acc update device(nstrbndloc,strbndglob(:)) async
        endif
      endif
c
c     UREY-BRADLEY
c
      if (use_urey) then
        !write(0,*) "push UREY-BRADLEY"
        nureyloc = bead%nureyloc
        ureyglob = bead%ureyglob
        if(updateGPU) then
!$acc update device(nureyloc,ureyglob(:)) async
        endif
      endif
c
c     ANGLE-ANGLE
c
      if (use_angang) then
        !write(0,*) "push ANGLE-ANGLE"
        nangangloc = bead%nangangloc
        angangglob = bead%angangglob
        if(updateGPU) then
!$acc update device(nangangloc,angangglob(:)) async
        endif
      endif
c
c     OP-BENDING
c
      if (use_opbend) then
        !write(0,*) "push OP-BENDING"
        nopbendloc = bead%nopbendloc
        opbendglob = bead%opbendglob
        if(updateGPU) then
!$acc update device(nopbendloc,opbendglob(:)) async
        endif
      endif
c
c     OP-DIST
c
      if (use_opdist) then
        !write(0,*) "push  OP-DIST"
        nopdistloc = bead%nopdistloc
        opdistglob = bead%opdistglob
        if(updateGPU) then
!$acc update device(nopdistloc,opdistglob(:)) async
        endif
      endif
c
c     IMPROP
c
      if (use_improp) then
       ! write(0,*) "push IMPROP"
        niproploc = bead%niproploc
        impropglob = bead%impropglob
        if(updateGPU) then
!$acc update device(niproploc,impropglob(:)) async
        endif
      endif
c
c     IMPTOR
c
      if (use_imptor) then
       ! write(0,*) "push IMPTOR"
        nitorsloc = bead%nitorsloc
        imptorglob = bead%imptorglob
        if(updateGPU) then
!$acc update device(nitorsloc,imptorglob(:)) async
        endif
      endif
c
c     TORSION
c
      if (use_tors) then
        !write(0,*) "push TORSION"
        ntorsloc = bead%ntorsloc
        torsglob = bead%torsglob
        if(updateGPU) then
!$acc update device(ntorsloc,torsglob(:)) async
        endif
      endif
c
c     PITORSION
c
      if (use_pitors) then
        !write(0,*) "push PITORSION"
        npitorsloc = bead%npitorsloc
        pitorsglob = bead%pitorsglob
        if(updateGPU) then
!$acc update device(npitorsloc,pitorsglob(:)) async
        endif
      endif
c
c     STRETCH-TORSION
c
      if (use_strtor) then
        !write(0,*) "push STRETCH-TORSION"
        nstrtorloc = bead%nstrtorloc
        strtorglob = bead%strtorglob
        if(updateGPU) then
!$acc update device(nstrtorloc,strtorglob(:)) async
        endif
      endif
c
c     TORSION-TORSION
c
      if (use_tortor) then
        !write(0,*) "push TORSION-TORSION"
        ntortorloc = bead%ntortorloc
        tortorglob = bead%tortorglob
        if(updateGPU) then
!$acc update device(ntortorloc,tortorglob(:)) async
        endif
      endif
c
c     ANGLE
c
      if (use_angle) then
        !write(0,*) "push ANGLE"
        nangleloc = bead%nangleloc
        angleglob = bead%angleglob
        if(updateGPU) then
!$acc update device(nangleloc,angleglob(:)) async
        endif
      endif
c
c     CHARGE
c
      if (use_charge) then
        !write(0,*) "push CHARGE"
        nionrecloc = bead%nionrecloc
        chgrecglob = bead%chgrecglob
        nionloc = bead%nionloc
        chgglob = bead%chgglob
        chgglobnl = bead%chgglobnl
        nionlocnl = bead%nionlocnl
        if(updateGPU) then
!$acc update device(nionloc,nionlocnl,nionrecloc) async
!$acc update device(chgrecglob(:),chgglob(:),chgglobnl(:)) async
        endif
      endif
      if (use_charge.or.use_mpole) then
        nelst = bead%nelst
        elst = bead%elst
        if(updateGPU) then
!$acc update device(nelst(:),elst(:,:)) async
        endif
      endif
c
c     MULTIPOLE
c
      if (use_mpole) then
        !write(0,*) "push MULTIPOLE"
        npolerecloc = bead%npolerecloc
        polerecglob = bead%polerecglob
        npoleloc = bead%npoleloc
        npolebloc = bead%npolebloc
        poleglob = bead%poleglob
        poleloc = bead%poleloc
        poleglobnl = bead%poleglobnl
        npolelocnl = bead%npolelocnl
        if(updateGPU) then
!$acc update device(npoleloc,npolebloc,poleglobnl) async
!$acc update device(npolelocnl,npolerecloc) async
!$acc update device(polerecglob(:),poleglob(:),poleloc(:)) async
        endif
      endif
c
c     POLARIZATION
c
      if (use_polar) then
        !write(0,*) "push POLARIZATION"
        nualt = bead%nualt
        udalt = bead%udalt
        upalt = bead%upalt
        uind = bead%uind
        uinp = bead%uinp
        if(updateGPU) then
!$acc update device(udalt(:,:,:),upalt(:,:,:)) async
!$acc update device(uind(:,:),uinp(:,:)) async
!$acc wait
        endif
      endif
c
c     TIME
c 
      timestep = bead%timestep

      ! write(0,*) "pushbead done"

      !SYNCHRONIZE GPU with CPU
      if(updateGPU) then 
!$acc wait
      endif

      end subroutine pushbead


      subroutine compute_observables_pi(pos,vel,forces,istep,dt)
      use atoms
      use units
      use atmtyp
      use math
      use boxes
      use mdstuf
      use bath
      use mpi
      use domdec
      !use qtb
      IMPLICIT NONE
      real(8), intent(in),allocatable :: pos(:,:,:),vel(:,:,:)
     &                                   ,forces(:,:,:)
      real(r_p), intent(in) :: dt
      integer, intent(in) :: istep
      real(8), allocatable :: centroid(:,:),vel_centroid(:,:)
      real(8) :: omp,omp2,dedv_mean
      integer :: ibead,i,j,k,ierr
      real(8) :: buffer_energy(4)

c
c     reduce potential and kinetic energies
c
      buffer_energy=0
      DO ibead=1,nbeadsloc        
        buffer_energy(1)=buffer_energy(1) + beadsloc(ibead)%eksum
        buffer_energy(2)=buffer_energy(2) + beadsloc(ibead)%dedv
        if(contract) then
          buffer_energy(3)=buffer_energy(3) + beadsloc(ibead)%eintra
        else
          buffer_energy(3)=buffer_energy(3) + beadsloc(ibead)%epot
        endif
      ENDDO

      if(contract) then
        DO ibead=1,nbeadsloc_ctr
          buffer_energy(4)=buffer_energy(4) + beadsloc_ctr(ibead)%einter
        ENDDO
      endif

      if ((ranktot.eq.0).and.(contract)) then
         call MPI_REDUCE(MPI_IN_PLACE,buffer_energy,4,MPI_REAL8
     $      ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         eksumpi_loc=buffer_energy(1)/(nbeads*nproc)
         dedv_mean=buffer_energy(2)/(nbeads*nproc)
         eintrapi_loc=buffer_energy(3)/(nbeads*nproc)
         einterpi_loc=buffer_energy(4)/(nbeads_ctr*nproc)
         epotpi_loc =  eintrapi_loc +  einterpi_loc
      else if ((ranktot.eq.0).and.(contract.eqv..false.)) then
         call MPI_REDUCE(MPI_IN_PLACE,buffer_energy,3,MPI_REAL8
     $      ,MPI_SUM,0,MPI_COMM_WORLD,ierr)         
         eksumpi_loc=buffer_energy(1)/(nbeads*nproc)
         dedv_mean=buffer_energy(2)/(nbeads*nproc)
         epotpi_loc=buffer_energy(3)/(nbeads*nproc)
      endif

      if((ranktot.ne.0).and.(contract)) then
          call MPI_REDUCE(buffer_energy,buffer_energy,4,MPI_REAL8
     $     ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      else if ((ranktot.ne.0).and.(contract.eqv..false.)) then
        call MPI_REDUCE(buffer_energy,buffer_energy,3,MPI_REAL8
     $     ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
      

      if(ranktot.eq.0) then
        allocate(centroid(3,n),vel_centroid(3,n))
        centroid(:,:)=0.d0
        vel_centroid(:,:)=0.d0
        DO ibead=1,nbeads
          DO i=1,n
            centroid(:,i)=centroid(j,i)+pos(:,i,ibead)
            vel_centroid(:,i)=vel_centroid(:,i)+vel(:,i,ibead)
          ENDDO
        ENDDO  
        centroid(:,:)=centroid(:,:)/REAL(nbeads)
        vel_centroid(:,:)=vel_centroid(:,:)/REAL(nbeads) 

        Ekcentroid=0.d0
        DO i=1,n ; DO j=1,3          
          Ekcentroid=Ekcentroid+mass(i)*vel_centroid(j,i)**2
        ENDDO ; ENDDO
        Ekcentroid=Ekcentroid*nbeads/convert

        !if (ir) then
        !  k = mod(istep-1,nseg)+1
        !  vad(:,:,k)=vel_centroid(:,:)
        !  if ((mod(istep,nseg).eq.0)) then
        !      call irspectra_pimd
        !      compteur=compteur+1
        ! endif
        !endif

        deallocate(vel_centroid)

        omp=nbeads*boltzmann*kelvin/hbar        
        omp2=omp*omp

c       COMPUTE PRIMITIVE KINETIC ENERGY
        ekprim=0.d0
        DO ibead=1,nbeads-1
          DO i=1,n
        !    if (atomic(i).eq.0) cycle
            DO j=1,3
              ekprim = ekprim - 0.5*mass(i)*omp2
     &          *(pos(j,i,ibead+1)-pos(j,i,ibead))**2
            ENDDO
          ENDDO
        ENDDO  
        DO i=1,n
        ! if (atomic(i).eq.0) cycle
          DO j=1,3
            ekprim = ekprim - 0.5*mass(i)*omp2
     &          *(pos(j,i,nbeads)-pos(j,i,1))**2
          ENDDO
        ENDDO  
        ekprim = (ekprim/nbeads
     &          + 0.5*nbeads*nfree*boltzmann*kelvin)/convert

c       COMPUTE VIRIAL KINETIC ENERGY
        ekvir=0.d0
        DO ibead=1,nbeads
          DO i=1,n
        !  if (atomic(i).eq.0) cycle
            DO j=1,3
              ekvir=ekvir+(pos(j,i,ibead)-centroid(j,i))
     &                      *forces(j,i,ibead)
            ENDDO
          ENDDO
        ENDDO

        presvir = prescon*( -dedv_mean + (Ekcentroid
     &               - ekvir)/(3*nbeads*volbox) )

        ekvir=0.5d0*(nfree*boltzmann*kelvin/convert-ekvir/nbeads)
        temppi = 2.0d0 * ekvir / (nfree * gasconst)
        temppi_cl = 2.0d0 * eksumpi_loc / (nfree * gasconst)

      endif

      end subroutine compute_observables_pi
      

      subroutine rescale_box_pi(istep)
c      rescale the simulation box according to extvol computed by ranktot 0
        use mpi
        use domdec
        use bath
        use boxes, only: volbox,xbox,ybox,zbox
        implicit none
        integer, intent(in) :: istep
        integer :: ierr
        real(8) :: third,scale

c       broadcast new volume  
        call MPI_BCAST(extvol,1,MPI_TPREC,0,MPI_COMM_WORLD,ierr)
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c       rescale box
        third=1.d0/3.d0
        scale=(extvol/extvolold)**third
        xbox = (extvol)**third
        ybox = (extvol)**third
        zbox = (extvol)**third  
        
!$acc update device(xbox,ybox,zbox)
c
c     propagate the new box dimensions to other lattice values
c 
        call lattice   
        call ddpme3dnpt(scale,istep)
      end subroutine rescale_box_pi

c ---------------------------------------------------------------------
c    COMMUNICATION ROUTINES

      subroutine prepare_loaded_bead (istep)
      use atmtyp
      use atomsMirror
      use cutoff
      use domdec
      use energi
      use freeze
      use mdstuf
      use moldyn
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer, intent(in) :: istep
      real*8 time0,time1
      real*8 oterm,hterm
      integer :: i,iglob

!$acc data  present(x,y,z,xold,yold,zold,v,a,mass,glob,use)

c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
      call reassign
c
c     -> reciprocal space
c
      call reassignpme(.false.)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
c
c     communicate positions
c
      time0 = mpi_wtime()
      call commpos
      call commposrec
      time1 = mpi_wtime()
      timecommstep = timecommstep + time1 - time0

      call reCast_position
c
c
      call reinitnl(istep)
c
      time0 = mpi_wtime()
      call mechanicstep(istep)
      time1 = mpi_wtime()
c
      timeparam = timeparam + time1 - time0
c
      time0 = mpi_wtime()
      call allocstep
      time1 = mpi_wtime()
      timeclear = timeclear  + time1 - time0
c      write(*,*) 'x 0 = ',x(1),y(1),v(1,1),a(1,1)
c      write(*,*) 'istep 1 = ',istep

      !rebuild the neighbor lists
      if (use_list) call nblist(istep)

!$acc end data
          
      end subroutine prepare_loaded_bead



      subroutine get_polymer_info(polymer,beads
     &   ,get_nloccomm,get_globbeadcomm)
      use domdec
      use mpi
      use atoms
      use atmtyp
      use units
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      type(BEAD_TYPE), intent(in) :: beads(:)
      LOGICAL, intent(in) :: get_nloccomm,get_globbeadcomm 
      integer i,l,iproc,ibead,ierr,iglob
      integer status(MPI_STATUS_SIZE),tagmpi
      integer, allocatable :: reqrec(:),reqsend(:),buffer(:,:)
      integer :: nsend,isend

      allocate(reqsend(nproctot))
      allocate(reqrec(nproctot))

      if(allocated(polymer%nbeadscomm)) deallocate(polymer%nbeadscomm)
      nsend=0
      if(get_nloccomm) then        
        if(allocated(polymer%nloccomm)) deallocate(polymer%nloccomm)
        nsend=nsend+1
      endif
      if(get_globbeadcomm) then      
        if(allocated(polymer%globbeadcomm)) then 
          deallocate(polymer%globbeadcomm)
        endif
        nsend=nsend+1
      endif      

c     first get the number of beads per process
c   
      if (ranktot.eq.0) then
        allocate(polymer%nbeadscomm(nproctot))  
        polymer%nbeadscomm(1)=size(beads)
        do i = 1, nproctot-1
          tagmpi = nproctot*ranktot + i + 1
          call MPI_IRECV(polymer%nbeadscomm(i+1),1,
     $     MPI_INT,i,tagmpi,MPI_COMM_WORLD,reqrec(i),ierr)
          call MPI_WAIT(reqrec(i),status,ierr)
c          write(*,*) 'nbeads of ',i,' = ',nbeadscomm(i+1)
        end do
       
      else
        tagmpi = ranktot + 1
        call MPI_ISEND(size(beads),1,MPI_INT,0,tagmpi,MPI_COMM_WORLD,
     $   reqsend(1),ierr)
        call MPI_WAIT(reqsend(1),status,ierr)
      end if

c
c     get the number of atoms per process
c   
      if(nsend==0) return

      if (ranktot.eq.0) then        
        polymer%nbeadsloc_max = maxval(polymer%nbeadscomm)
        if(get_nloccomm) then
          allocate(polymer%nloccomm(polymer%nbeadsloc_max,nproctot))
          do i=1,size(beads) 
            polymer%nloccomm(i,1)=beads(i)%nloc
          enddo
        endif
        if(get_globbeadcomm) then
          allocate(polymer%globbeadcomm(polymer%nbeadsloc_max,nproctot))
          do i=1,size(beads) 
            polymer%globbeadcomm(i,1)=beads(i)%ibead_glob
          enddo
        endif        
        do i = 1, nproctot-1
          allocate(buffer(polymer%nbeadscomm(i+1),nsend))
          tagmpi = i + 1
          call MPI_IRECV(buffer,nsend*polymer%nbeadscomm(i+1),
     $     MPI_INT,i,tagmpi,MPI_COMM_WORLD,reqrec(i),ierr)
          call MPI_WAIT(reqrec(i),status,ierr)
c          write(*,*) 'nloc of ',i,' = ',nloccomm(i+1)
          isend=1
          if(get_nloccomm) then
            polymer%nloccomm(1:polymer%nbeadscomm(i+1),i+1)
     &         =buffer(:,isend)
            isend=isend+1
          endif
          if(get_globbeadcomm) then
            polymer%globbeadcomm(1:polymer%nbeadscomm(i+1),i+1)
     &         =buffer(:,isend)
            isend=isend+1
          endif
          deallocate(buffer)
        end do
      else
        tagmpi = ranktot + 1
        allocate(buffer(size(beads),nsend))
        isend=1
        if(get_nloccomm) then
          do i=1,size(beads) 
            buffer(i,isend)=beads(i)%nloc
          enddo
          isend=isend+1
        endif
        if(get_globbeadcomm) then
          do i=1,size(beads) 
            buffer(i,isend)=beads(i)%ibead_glob
          enddo
          isend=isend+1
        endif
        call MPI_ISEND(buffer,nsend*size(beads),MPI_INT
     $   ,0,tagmpi,MPI_COMM_WORLD,reqsend(1),ierr)
        call MPI_WAIT(reqsend(1),status,ierr)
      end if

      end subroutine get_polymer_info


      subroutine gather_polymer(polymer,beads, 
     &                send_pos,send_vel,send_forces)
      use domdec
      use mpi
      use atoms
      use atmtyp
      use units
      implicit none
      type(BEAD_TYPE), intent(in) :: beads(:)
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      ! pos_full,vel_full,forces_full will be allocated for ranktot 0 (and deallocated for all the other procs)
      ! get_repart_dof_beads must be called before to fill out the polymer info concerning parallelization
      LOGICAL, intent(in) :: send_pos,send_vel,send_forces
      integer i,j,k,l,iproc,ibead,ierr,iglob
      integer status(MPI_STATUS_SIZE),tagmpi
      real(r_p), allocatable :: buffer(:,:),indexposcomm(:,:,:,:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer :: nsend,isend,nloc_max

      allocate(reqsend(nproctot))
      allocate(reqrec(nproctot))

      nsend=1 ! always send repart_dof_beads
      if(send_pos) nsend = nsend + 3
      if(send_vel) nsend = nsend + 3
      if(send_forces) nsend = nsend + 3

      if(allocated(polymer%repart_dof_beads))then
        deallocate(polymer%repart_dof_beads)
      endif
      if(send_pos) then
        if(allocated(polymer%pos)) deallocate(polymer%pos)
      endif
      if(send_vel) then
        if(allocated(polymer%vel)) deallocate(polymer%vel)
      endif
      if(send_forces) then
        if(allocated(polymer%forces)) deallocate(polymer%forces)
      endif

c     gather number of beads and degrees of freedom
      call get_polymer_info(polymer,beads,.TRUE.,.TRUE.)

      if (ranktot.eq.0) then
        ! number of beads in the polymer is the sum of the communicated beads 
        ! divided by the spatial number of procs
        polymer%nbeads = SUM(polymer%nbeadscomm)/nproc        
        allocate(polymer%repart_dof_beads(n,polymer%nbeads,nproctot))
        allocate(indexposcomm(nsend,n,polymer%nbeads,nproctot))        
      end if     
      
c
c     get their indexes and positions, velocities and forces
c   
      if (ranktot.eq.0) then
        do i = 1, nproctot-1
          tagmpi = i + 1
          do k=1,polymer%nbeadscomm(i+1)
            call MPI_IRECV(indexposcomm(1,1,k,i+1),
     $     nsend*polymer%nloccomm(k,i+1),
     $     MPI_RPREC,i,tagmpi,MPI_COMM_WORLD,reqrec(i),ierr)
            call MPI_WAIT(reqrec(i),status,ierr)
          enddo
        end do
      else
        tagmpi = ranktot + 1        
        do k = 1, size(beads)
         allocate(buffer(nsend,beads(k)%nloc))
         buffer = 0.d0
         do i = 1, beads(k)%nloc
          iglob = beads(k)%glob(i)
          buffer(1,i) = real(iglob,r_p)
          isend=2
          if(send_pos) then
            buffer(isend,i) = beads(k)%x(iglob)
            buffer(isend+1,i) = beads(k)%y(iglob)
            buffer(isend+2,i) = beads(k)%z(iglob)
            isend = isend + 3
          endif

          if(send_vel) then
            buffer(isend,i) = beads(k)%v(1,iglob)
            buffer(isend+1,i) = beads(k)%v(2,iglob)
            buffer(isend+2,i) = beads(k)%v(3,iglob)
            isend = isend + 3
          endif

          if(send_forces) then
            buffer(isend,i)=beads(k)%a(1,iglob)*mass(iglob)/convert
            buffer(isend+1,i)=beads(k)%a(2,iglob)*mass(iglob)/convert
            buffer(isend+2,i)=beads(k)%a(3,iglob)*mass(iglob)/convert
            isend = isend + 3
          endif
         end do
         call MPI_ISEND(buffer,nsend*beads(k)%nloc,MPI_RPREC
     $   ,0,tagmpi,MPI_COMM_WORLD,reqsend(1),ierr)
         call MPI_WAIT(reqsend(1),status,ierr)

         deallocate(buffer)
        enddo
      end if
      
      if(ranktot .eq. 0) then
c       ORDER POSITIONS AND VELOCITES IN pos_full AND vel_full
        if(send_pos) allocate(polymer%pos(3,n,polymer%nbeads))
        if(send_vel) allocate(polymer%vel(3,n,polymer%nbeads))
        if(send_forces) allocate(polymer%forces(3,n,polymer%nbeads))
          
        DO k=1,size(beads)
          DO i=1,beads(k)%nloc
            iglob=beads(k)%glob(i)
            if(send_pos) then
              polymer%pos(1,iglob,k) =beads(k)%x(iglob)
              polymer%pos(2,iglob,k) =beads(k)%y(iglob)
              polymer%pos(3,iglob,k) =beads(k)%z(iglob)
            endif
            if(send_vel) polymer%vel(:,iglob,k) =beads(k)%v(:,iglob)
            if(send_forces) then
              polymer%forces(:,iglob,k)= 
     &             beads(k)%a(:,iglob)*mass(iglob)/convert
            endif
          ENDDO
        ENDDO

        do iproc = 1, nproctot-1
          do k=1,polymer%nbeadscomm(iproc+1)
            ibead = polymer%globbeadcomm(k,iproc+1)
            DO i=1,polymer%nloccomm(k,iproc+1)      
              isend=1   
              
              iglob=nint(indexposcomm(1,i,k,iproc+1))
              polymer%repart_dof_beads(i,k,iproc+1)=iglob
              isend=1
              if(send_pos) then
                DO j=1,3     
                 polymer%pos(j,iglob,ibead)=
     &               indexposcomm(j+isend,i,k,iproc+1)  
                ENDDO
                isend=isend+3
              endif
              if(send_vel) then
                DO j=1,3     
                  polymer%vel(j,iglob,ibead)=
     &               indexposcomm(j+isend,i,k,iproc+1)  
                ENDDO
                isend=isend+3
              endif
              if(send_forces) then
                DO j=1,3
                  polymer%forces(j,iglob,ibead)=
     &               indexposcomm(j+isend,i,k,iproc+1)  
                ENDDO
                isend=isend+3
              endif
            ENDDO 
          enddo
        enddo

        deallocate(indexposcomm)        
      endif

      deallocate(reqsend,reqrec)

      end subroutine gather_polymer

      subroutine broadcast_polymer(polymer,beads, 
     &                send_pos,send_vel,send_forces)
      use domdec
      use mpi
      use atoms
      use atmtyp
      use units
      implicit none
      type(BEAD_TYPE), intent(inout) :: beads(:)
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      LOGICAL, intent(in) :: send_pos,send_vel,send_forces
      integer i,j,k,l,iproc,ibead,ierr,iglob
      integer status(MPI_STATUS_SIZE),tagmpi
      real(r_p), allocatable :: buffer(:,:),indexposcomm(:,:,:,:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer :: nsend,isend,nloc_max

      allocate(reqsend(nproctot))
      allocate(reqrec(nproctot))

      nsend=0
      if(send_pos) nsend = nsend + 3
      if(send_vel) nsend = nsend + 3
      if(send_forces) nsend = nsend + 3
      if(nsend==0 .and. ranktot==0) then
        write(0,*) "Error: broadcast_polymer ",
     &     "was called without requesting any send"
        call fatal
      endif
      
      if(ranktot.eq.0) then
        if(.not.allocated(polymer%repart_dof_beads)) then
          write(0,*) "Error: repart_dof_beads not allocated"
     &     ," in broadcast_polymer"
          call fatal
        endif
c       PUT BACK POSITIONS AND VELOCITES AT CORRECT INDICES
        allocate(indexposcomm(nsend,n,polymer%nbeads,nproctot))
        DO k=1,size(beads)
          DO i=1,beads(k)%nloc 
            iglob=beads(k)%glob(i)
            if(send_pos) then
              beads(k)%x(iglob)=polymer%pos(1,iglob,k)
              beads(k)%y(iglob)=polymer%pos(2,iglob,k)
              beads(k)%z(iglob)=polymer%pos(3,iglob,k)
            endif
            if(send_vel) then
              beads(k)%v(:,iglob)=polymer%vel(:,iglob,k)
            endif
            if(send_forces) then
              beads(k)%a(:,iglob) = 
     &           polymer%forces(:,iglob,k)/mass(iglob)*convert
            endif
          ENDDO
        ENDDO

        do iproc = 1, nproctot-1
          do k=1,polymer%nbeadscomm(iproc+1)
            ibead = polymer%globbeadcomm(k,iproc+1)
            DO i=1,polymer%nloccomm(k,iproc+1)
              iglob=polymer%repart_dof_beads(i,k,iproc+1)
              isend=0
              if(send_pos) then
                DO j=1,3
                  indexposcomm(j+isend,i,k,iproc+1)=
     &               polymer%pos(j,iglob,ibead) 
                ENDDO
                isend=isend+3
              endif
              if(send_vel) then
                DO j=1,3
                  indexposcomm(j+isend,i,k,iproc+1)=
     &               polymer%vel(j,iglob,ibead) 
                ENDDO
                isend=isend+3
              endif
              if(send_forces) then
                DO j=1,3
                  indexposcomm(j+isend,i,k,iproc+1)=
     &               polymer%forces(j,iglob,ibead) 
                ENDDO
                isend=isend+3
              endif
            ENDDO
          enddo     
        enddo   
      endif
c
c     communicate back positions
c
      if (ranktot.eq.0) then
        do i = 1, nproctot-1
          tagmpi = i + 1
          do k=1,polymer%nbeadscomm(i+1)
            call MPI_ISEND(indexposcomm(1,1,k,i+1)
     $     ,nsend*polymer%nloccomm(k,i+1)
     $     ,MPI_RPREC,i,tagmpi,MPI_COMM_WORLD,reqsend(i),ierr)
            call MPI_WAIT(reqsend(i),status,ierr)
          enddo
        end do

        deallocate(indexposcomm)
      else        
        
        tagmpi = ranktot + 1        
        do k = 1, size(beads)
          allocate(buffer(nsend,beads(k)%nloc))
          call MPI_IRECV(buffer,nsend*beads(k)%nloc,MPI_RPREC
     $   ,0,tagmpi,MPI_COMM_WORLD,reqrec(1),ierr)
          call MPI_WAIT(reqrec(1),status,ierr)
          do i = 1, beads(k)%nloc
            iglob = beads(k)%glob(i)
            isend=0
            if(send_pos) then
              beads(k)%x(iglob) = buffer(isend+1,i)
              beads(k)%y(iglob) = buffer(isend+2,i)
              beads(k)%z(iglob) = buffer(isend+3,i)
              isend=isend+3
            endif
            if(send_vel) then
              DO j=1,3
                beads(k)%v(j,iglob) = buffer(j+isend,i)
              ENDDO
              isend=isend+3
            endif
            if(send_forces) then
              DO j=1,3
                beads(k)%a(j,iglob) = buffer(j+isend,i)
     &            /mass(iglob)*convert
              ENDDO
              isend=isend+3
            endif
           ! write(*,*) 'j = ',glob(j),'k = ',k,'x = ',pospi(1,glob(j),k)
          end do
          deallocate(buffer)
        end do
      end if

      deallocate(reqsend,reqrec)

      end subroutine broadcast_polymer

c     -----------------------------------------
c      CONTRACTION SUBROUTINES
      

      subroutine contract_polymer(polymer,polymer_ctr)
        use atoms
        implicit none
        type(POLYMER_COMM_TYPE), intent(in) :: polymer
        type(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
        integer :: i,j,k

        if(allocated(polymer_ctr%pos)) 
     &               deallocate(polymer_ctr%pos)
        allocate(polymer_ctr%pos(3,n,nbeads_ctr))

        DO i=1,n ; DO j=1,3
          polymer_ctr%pos(j,i,:)=matmul(contractor_mat
     &                    ,polymer%pos(j,i,:))
        ENDDO ; ENDDO        

        !do k=1,nbeads_ctr
        !  do i=1,n
        !    write(17+k,*) i,polymer_ctr%pos(:,i,k)
        !  enddo
        !  FLUSH(17+k)
        !enddo
        !call fatal

      end subroutine contract_polymer

      subroutine project_forces_contract(polymer,polymer_ctr)
        use atoms
        implicit none
        type(POLYMER_COMM_TYPE),intent(inout) :: polymer
        type(POLYMER_COMM_TYPE),intent(in) :: polymer_ctr
        integer :: i,j

        if(allocated(polymer%forces_slow)) then
          deallocate(polymer%forces_slow)
        endif
        allocate(polymer%forces_slow(3,n,nbeads))
      
        DO i=1,n ; DO j=1,3
          polymer%forces_slow(j,i,:)=matmul(uncontractor_mat
     &         ,polymer_ctr%forces(j,i,:))
        ENDDO ; ENDDO

        if(allocated(polymer%forces)) then
          polymer%forces=polymer%forces+polymer%forces_slow
        endif

       
      end subroutine project_forces_contract

      end module
