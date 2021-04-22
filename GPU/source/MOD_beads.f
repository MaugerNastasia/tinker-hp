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

        !NBLIST
        integer, allocatable :: ineignl(:)

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
        integer nionlocloop,nionlocnlloop
        integer nionlocnlb,nionlocnlb_pair,nionlocnlb2_pair
     &       ,nshortionlocnlb2_pair
        integer, allocatable :: chgrecglob(:)
        integer, allocatable ::chgglob(:)
        integer, allocatable :: chgglobnl(:)
        integer, allocatable :: nelst(:),elst(:,:),shortelst(:,:)
        integer, allocatable :: nelstc(:),nshortelst(:),nshortelstc(:)
        integer, allocatable :: eblst(:),ieblst(:)
        integer, allocatable :: shorteblst(:),ishorteblst(:)
        integer,allocatable :: celle_key(:),celle_glob(:)
     &              ,celle_pole(:),celle_loc(:),celle_ploc(:)
     &              ,celle_plocnl(:)
        integer,allocatable :: celle_chg(:)
        real(t_p),allocatable:: celle_x(:),celle_y(:),celle_z(:)
        !
        !     MULTIPOLE
        !
        integer :: npolerecloc,npoleloc,npolebloc, npolelocnl
        integer :: npolelocnlb,npolelocnlb_pair
        integer :: npolelocnlb2_pair,nshortpolelocnlb2_pair
        integer, allocatable :: polerecglob(:)
        integer, allocatable :: poleglob(:)
        integer, allocatable :: poleloc(:),polelocnl(:),polerecloc(:)
        integer, allocatable :: poleglobnl(:)
        integer npolerecloc_old,npolereclocloop
        integer npolelocloop, npolelocnlloop,npoleblocloop
        !
        !      POLARIZATION
        !
        integer :: nualt
        real(t_p), allocatable :: udalt(:,:,:),upalt(:,:,:)
        real(t_p), allocatable :: uind(:,:),uinp(:,:)
        !
        !     VDW
        !
        integer :: nvdwloc,nvdwbloc,nvdwlocnl,nvdwlocnlb,nvdwlocnlb_pair
        integer :: nvdwblocloop
        integer :: nvdwlocnlb2_pair,nshortvdwlocnlb2_pair
        integer, allocatable :: vdwglob(:)
        integer, allocatable :: vdwglobnl(:),vdwlocnl(:)
        integer, allocatable :: nvlst(:),vlst(:,:)
        integer, allocatable :: nshortvlst(:),shortvlst(:,:)
        integer, allocatable :: vblst(:),ivblst(:)
        integer, allocatable :: shortvblst(:),ishortvblst(:)
        integer,allocatable :: cellv_key(:),cellv_glob(:)
     &                        ,cellv_loc(:),cellv_jvdw(:)


        !
        !     SCALING FACTORS
        !
        integer :: n_vscale,n_mscale,n_cscale
        integer :: n_uscale,n_dpscale,n_dpuscale
        integer,allocatable :: vcorrect_ik(:,:)
        real(t_p), allocatable :: vcorrect_scale(:)
        integer,allocatable :: mcorrect_ik(:,:)
        real(t_p), allocatable :: mcorrect_scale(:)
        integer,allocatable :: ccorrect_ik(:,:)
        real(t_p), allocatable :: ccorrect_scale(:)
        integer,allocatable :: ucorrect_ik(:)
        real(t_p), allocatable :: ucorrect_scale(:)
        integer,allocatable :: dpcorrect_ik(:)
        real(t_p), allocatable :: dpcorrect_scale(:)
        integer,allocatable :: dpucorrect_ik(:)
        real(t_p), allocatable :: dpucorrect_scale(:)
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
!$acc enter data copyin(eigmat,eigmattr,omkpi)
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
!        write(0,*) "allocbea=",i,"/",nbeadsloc
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

      if(allocated(ineignl)) 
     &   allocate(bead%ineignl(n))
         


      if (use_vdw) then
        allocate(bead%vdwglob(n))
        allocate(bead%vdwglobnl(n))
        allocate(bead%vdwlocnl(n))
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
        allocate(bead%udalt(3,n,maxualt))
        allocate(bead%upalt(3,n,maxualt))
        allocate(bead%uind(3,n))
        allocate(bead%uinp(3,n))
      end if

      if (use_mpole) then
        allocate(bead%polerecglob(n))
        allocate(bead%poleglob(n))
        allocate(bead%poleloc(n))
        allocate(bead%poleglobnl(n))
        allocate(bead%polelocnl(n))
        allocate(bead%polerecloc(n))
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
        if (allocated(bead%vdwlocnl)) deallocate(bead%vdwlocnl)
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
        if (allocated(bead%polerecloc)) deallocate(bead%polerecloc)
        if (allocated(bead%polelocnl)) deallocate(bead%polelocnl)
        if (allocated(bead%poleglobnl)) deallocate(bead%poleglobnl)

        if (allocated(bead%nvlst)) deallocate(bead%nvlst)
        if (allocated(bead%vlst)) deallocate(bead%vlst)
        if (allocated(bead%nelst)) deallocate(bead%nelst)
        if (allocated(bead%elst)) deallocate(bead%elst)

        if (allocated(bead%nelstc)) deallocate(bead%nelstc)
        if (allocated(bead%shortelst)) deallocate(bead%shortelst)
        if (allocated(bead%nshortelst)) deallocate(bead%nshortelst)
        if (allocated(bead%shortelst)) deallocate(bead%shortelst)
        if (allocated(bead%nshortelstc)) deallocate(bead%nshortelstc)

        if (allocated(bead%eblst)) deallocate(bead%eblst)
        if (allocated(bead%ieblst)) deallocate(bead%ieblst)
        if (allocated(bead%shorteblst)) deallocate(bead%shorteblst)
        if (allocated(bead%ishorteblst)) deallocate(bead%ishorteblst)

        if (allocated(bead%celle_key)) deallocate(bead%celle_key)
        if (allocated(bead%celle_glob)) deallocate(bead%celle_glob)
        if (allocated(bead%celle_pole)) deallocate(bead%celle_pole)
        if (allocated(bead%celle_loc)) deallocate(bead%celle_loc)
        if (allocated(bead%celle_ploc)) deallocate(bead%celle_ploc)
        if (allocated(bead%celle_plocnl)) deallocate(bead%celle_plocnl)
        if (allocated(bead%celle_chg)) deallocate(bead%celle_chg)
        if (allocated(bead%celle_x)) deallocate(bead%celle_x)
        if (allocated(bead%celle_y)) deallocate(bead%celle_y)
        if (allocated(bead%celle_z)) deallocate(bead%celle_z)

        if (allocated(bead%nshortvlst)) deallocate(bead%nshortvlst)
        if (allocated(bead%shortvlst)) deallocate(bead%shortvlst)
        if (allocated(bead%vblst)) deallocate(bead%vblst)
        if (allocated(bead%ivblst)) deallocate(bead%ivblst)
        if (allocated(bead%shortvblst)) deallocate(bead%shortvblst)
        if (allocated(bead%ishortvblst)) deallocate(bead%ishortvblst)
        if (allocated(bead%cellv_key)) deallocate(bead%cellv_key)
        if (allocated(bead%cellv_glob)) deallocate(bead%cellv_glob)
        if (allocated(bead%cellv_loc)) deallocate(bead%cellv_loc)
        if (allocated(bead%cellv_jvdw)) deallocate(bead%cellv_jvdw)


      end subroutine deallocate_bead
      

      subroutine initbead(istep,bead,skip_parameters)
      use angang
      use angle
      use atomsMirror
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
      logical,intent(in) :: skip_parameters
      integer ibead,i,modnl

      !write(0,*) "initbead pos"
!$acc wait
c
c     positions, speed, mass
c
!$acc update host(x(:),y(:),z(:),v(:,:),a(:,:))

      bead%x = x
      bead%y = y
      bead%z = z
      bead%v = v
      bead%a = a

      !write(0,*) "initbead ener"

!$acc update host(eksumpi_loc,ekinpi_loc)
      if(contract) then
!$acc update host(eintrapi_loc,einterpi_loc)
        bead%eintra=eintrapi_loc
        bead%einter=einterpi_loc
        bead%epot=eintrapi_loc+einterpi_loc
        bead%etot=eintrapi_loc+einterpi_loc+eksumpi_loc
      else
!$acc update host(epotpi_loc)
        bead%epot=epotpi_loc
        bead%etot=epotpi_loc+eksumpi_loc
      endif

      bead%eksum=eksumpi_loc
      bead%ekin=ekinpi_loc
      bead%dedv = dedv

!$acc update host(glob(:))
      !write(0,*) "initbead glob",allocated(glob),allocated(bead%glob)

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
      !write(0,*) "initbead bloc",nblocrecdir

      bead%nbloc = nbloc
      bead%nlocrec = nlocrec
      bead%nblocrec = nblocrec
      bead%nlocnl = nlocnl
      bead%nblocrecdir = nblocrecdir

!$acc update host(loc(:),ineignl(:)) async
!$acc update host(repart(:),repartrec(:)) async
!$acc update host(domlenrec(:),domlenpole(:),domlenpolerec(:)) async
!$acc update host(globrec(:),locrec(:),globrec1(:),locrec1(:)) async
!$acc update host(bufbegpole(:)) async
!$acc update host(buflen2(:),buf1(:),buf2(:)) async
!$acc update host(bufbeg1(:),bufbeg2(:)) async
!$acc wait
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
      !write(0,*) "initbead vdw"

      if (use_vdw) then
!$acc update host(vdwglob(:),vdwglobnl(:),vdwlocnl(:))
        bead%nvdwbloc = nvdwbloc
        bead%vdwglob = vdwglob
        bead%vdwglobnl = vdwglobnl
        bead%vdwlocnl = vdwlocnl
        bead%nvdwlocnl = nvdwlocnl
        bead%nvdwlocnlb = nvdwlocnlb
        bead%nvdwlocnlb_pair = nvdwlocnlb_pair
        bead%nvdwlocnlb2_pair = nvdwlocnlb2_pair
        bead%nshortvdwlocnlb2_pair = nshortvdwlocnlb2_pair
        bead%nvdwblocloop = nvdwblocloop
      end if
c
c     BONDS
c
      if (use_bond) then
!$acc update host(bndglob(:))
        bead%bndglob = bndglob
        bead%nbondloc = nbondloc
      end if
c
c     STRETCH-BEND
c
      if (use_strbnd) then
!$acc update host(strbndglob(:))
        bead%strbndglob = strbndglob
        bead%nstrbndloc = nstrbndloc
      end if
c
c     UREY-BRADLEY
c
      if (use_urey) then
!$acc update host(ureyglob(:))
        bead%ureyglob = ureyglob
        bead%nureyloc = nureyloc
      end if
c
c     ANGlE-ANGLE
c
      if (use_angang) then
!$acc update host(angangglob(:))
        bead%angangglob = angangglob
        bead%nangangloc = nangangloc
      end if
c
c     OP-BENDING
c
      if (use_opbend) then
!$acc update host(opbendglob(:))
        bead%opbendglob = opbendglob
        bead%nopbendloc = nopbendloc
      end if
c
c     OP-DIST
c
      if (use_opdist) then
!$acc update host(opdistglob(:))
        bead%opdistglob = opdistglob
        bead%nopdistloc = nopdistloc
      end if
c
c     IMPROP
c
      if (use_improp) then
!$acc update host(impropglob(:))
        bead%impropglob = impropglob
        bead%niproploc = niproploc
      end if
c
c     IMPTOR
c
      if (use_imptor) then
!$acc update host(imptorglob(:))
        bead%imptorglob = imptorglob
        bead%nitorsloc = nitorsloc
      end if
c
c     TORSION
c
      if (use_tors) then
!$acc update host(torsglob(:))
        bead%torsglob = torsglob
        bead%ntorsloc = ntorsloc
      end if
c
c     PITORSION
c
      if (use_pitors) then
!$acc update host(pitorsglob(:))
        bead%pitorsglob = pitorsglob
        bead%npitorsloc = npitorsloc
      end if
c
c     STRETCH-TORSION
c
      if (use_strtor) then
!$acc update host(strtorglob(:))
        bead%strtorglob = strtorglob
        bead%nstrtorloc = nstrtorloc
      end if
c
c     TORSION-TORSION
c
      if (use_tortor) then
!$acc update host(tortorglob(:))
        bead%tortorglob = tortorglob
        bead%ntortorloc = ntortorloc
      end if
c
c     ANGLE
c
      if (use_angle) then
!$acc update host(angleglob(:))
        bead%angleglob = angleglob
        bead%nangleloc = nangleloc
      end if
c
c     CHARGE
c
      !write(0,*) "initbead charge"
      if (use_charge) then
!$acc update host(chgrecglob(:),chgglob(:),chgglobnl(:))
        bead%chgrecglob = chgrecglob
        bead%nionrecloc = nionrecloc
        bead%chgglob = chgglob
        bead%nionloc = nionloc
        bead%chgglobnl = chgglobnl
        bead%nionlocnl = nionlocnl
        bead%nionlocloop = nionlocloop
        bead%nionlocnlloop = nionlocnlloop
        bead%nionlocnlb = nionlocnlb
        bead%nionlocnlb_pair = nionlocnlb_pair
        bead%nionlocnlb2_pair = nionlocnlb2_pair
        bead%nshortionlocnlb2_pair = nshortionlocnlb2_pair
      end if
c      if ((use_charge.or.use_mpole).and.(modnl.eq.0)) then
c        nelstpi(1:nlocnl,ibead) = nelst
c        elstpi(:,1:nlocnl,ibead) = elst
c      end if
c
c     MULTIPOLE
c
      !,*) "initbead mpole"
      if (use_mpole) then
!$acc update host(polerecglob(:),poleglob(:),poleloc(:)) async
!$acc update host(polelocnl(:),polerecloc(:)) async
!$acc wait
        if(nlocrec>0) then
          bead%polerecglob(1:nlocrec) = polerecglob(:)
        endif
        bead%npolerecloc = npolerecloc
        bead%poleglob(1:nbloc) = poleglob(:)
        bead%poleloc = poleloc
        bead%polerecloc = polerecloc
        bead%npoleloc = npoleloc
        bead%npolebloc = npolebloc
        if(nlocnl>0) then
          bead%poleglobnl(1:nlocnl) = poleglobnl
        endif
        bead%npolelocnl = npolelocnl
        bead%polelocnl = polelocnl
        bead%npolelocnlb = npolelocnlb
        bead%npolelocnlb_pair = npolelocnlb_pair
        bead%npolelocnlb2_pair = npolelocnlb2_pair
        bead%nshortpolelocnlb2_pair = nshortpolelocnlb2_pair
        bead%npolerecloc_old = npolerecloc_old
        bead%npolelocloop = npolelocloop
        bead%npolelocnlloop = npolelocnlloop
        bead%npoleblocloop = npoleblocloop
        bead%npolereclocloop = npolereclocloop
        
      end if
c
c     POLARIZATION
c
      if (use_polar) then
!$acc update host(udalt(:,:,:),upalt(:,:,:)) async
!$acc update host(uind(:,:),uinp(:,:)) async
!$acc wait
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

      subroutine resize_nl_arrays_bead(bead)
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
      use polpot
      use chgpot
      use mplpot
      use vdwpot
      implicit none
      TYPE(BEAD_TYPE), intent(inout) :: bead
      integer nblocrecdirmax,modnl

!$acc wait

      if (allocated(nvlst)) then
        if (allocated(bead%nvlst)) deallocate(bead%nvlst)
        allocate(bead%nvlst(size(nvlst)))
      endif

      if(allocated(vlst)) then
        if (allocated(bead%vlst)) deallocate(bead%vlst)
        allocate(bead%vlst(size(vlst,1),size(vlst,2)))
      endif

      if(allocated(nelst)) then
        if (allocated(bead%nelst)) deallocate(bead%nelst)
        allocate(bead%nelst(size(nelst)))
      endif

      if (allocated(elst)) then        
        if (allocated(bead%elst)) deallocate(bead%elst)
        allocate(bead%elst(size(elst,1),size(elst,2)))
      end if

      if (allocated(nelstc)) then        
        if (allocated(bead%nelstc)) deallocate(bead%nelstc)
        allocate(bead%nelstc(size(nelstc)))
      end if

      if (allocated(shortelst)) then        
        if (allocated(bead%shortelst)) deallocate(bead%shortelst)
        allocate(bead%shortelst(size(shortelst,1),size(shortelst,2)))
      end if

      if (allocated(nshortelst)) then        
        if (allocated(bead%nshortelst)) deallocate(bead%nshortelst)
        allocate(bead%nshortelst(size(nshortelst)))
      end if

      if (allocated(nshortelstc)) then        
        if (allocated(bead%nshortelstc)) deallocate(bead%nshortelstc)
        allocate(bead%nshortelstc(size(nshortelstc)))
      end if

      if (allocated(eblst)) then        
        if (allocated(bead%eblst)) deallocate(bead%eblst)
        allocate(bead%eblst(size(eblst)))
      end if

       if (allocated(ieblst)) then        
        if (allocated(bead%ieblst)) deallocate(bead%ieblst)
        allocate(bead%ieblst(size(ieblst)))
      end if

      if (allocated(shorteblst)) then        
        if (allocated(bead%shorteblst)) deallocate(bead%shorteblst)
        allocate(bead%shorteblst(size(shorteblst)))
      end if

      if (allocated(ishorteblst)) then        
        if (allocated(bead%ishorteblst)) deallocate(bead%ishorteblst)
        allocate(bead%ishorteblst(size(ishorteblst)))
      end if

      if (allocated(nshortvlst)) then        
        if (allocated(bead%nshortvlst)) deallocate(bead%nshortvlst)
        allocate(bead%nshortvlst(size(nshortvlst)))
      end if

      if (allocated(shortvlst)) then        
        if (allocated(bead%shortvlst)) deallocate(bead%shortvlst)
        allocate(bead%shortvlst(size(shortvlst,1),size(shortvlst,2)))
      end if

      if (allocated(vblst)) then        
        if (allocated(bead%vblst)) deallocate(bead%vblst)
        allocate(bead%vblst(size(vblst)))
      end if

      if (allocated(ivblst)) then        
        if (allocated(bead%ivblst)) deallocate(bead%ivblst)
        allocate(bead%ivblst(size(ivblst)))
      end if

      if (allocated(shortvblst)) then        
        if (allocated(bead%shortvblst)) deallocate(bead%shortvblst)
        allocate(bead%shortvblst(size(shortvblst)))
      end if

      if (allocated(ishortvblst)) then        
        if (allocated(bead%ishortvblst)) deallocate(bead%ishortvblst)
        allocate(bead%ishortvblst(size(ishortvblst)))
      end if

      if(allocated(celle_glob))  then
        if (allocated(bead%celle_glob)) deallocate(bead%celle_glob)
        allocate(bead%celle_glob(size(celle_glob)))
      endif

      if(allocated(celle_pole))  then
        if (allocated(bead%celle_pole)) deallocate(bead%celle_pole)
        allocate(bead%celle_pole(size(celle_pole)))
      endif

      if(allocated(celle_plocnl))  then
        if (allocated(bead%celle_plocnl)) deallocate(bead%celle_plocnl)
        allocate(bead%celle_plocnl(size(celle_plocnl)))
      endif

      if(allocated(celle_key))  then
        if (allocated(bead%celle_key)) deallocate(bead%celle_key)
        allocate(bead%celle_key(size(celle_key)))
      endif

      if(allocated(celle_chg))  then
        if (allocated(bead%celle_chg)) deallocate(bead%celle_chg)
        allocate(bead%celle_chg(size(celle_chg)))
      endif

      if(allocated(celle_loc))  then
        if (allocated(bead%celle_loc)) deallocate(bead%celle_loc)
        allocate(bead%celle_loc(size(celle_loc)))
      endif

      if(allocated(celle_ploc))  then
        if (allocated(bead%celle_ploc)) deallocate(bead%celle_ploc)
        allocate(bead%celle_ploc(size(celle_ploc)))
      endif

      if(allocated(celle_x))  then
        if (allocated(bead%celle_x)) deallocate(bead%celle_x)
        allocate(bead%celle_x(size(celle_x)))
      endif

      if(allocated(celle_y))  then
        if (allocated(bead%celle_y)) deallocate(bead%celle_y)
        allocate(bead%celle_y(size(celle_y)))
      endif

      if(allocated(celle_z))  then
        if (allocated(bead%celle_z)) deallocate(bead%celle_z)
        allocate(bead%celle_z(size(celle_z)))
      endif

      if(allocated(cellv_key))  then
        if (allocated(bead%cellv_key)) deallocate(bead%cellv_key)
        allocate(bead%cellv_key(size(cellv_key)))
      endif

      if(allocated(cellv_glob))  then
        if (allocated(bead%cellv_glob)) deallocate(bead%cellv_glob)
        allocate(bead%cellv_glob(size(cellv_glob)))
      endif

      if(allocated(cellv_loc))  then
        if (allocated(bead%cellv_loc)) deallocate(bead%cellv_loc)
        allocate(bead%cellv_loc(size(cellv_loc)))
      endif

      if(allocated(cellv_jvdw))  then
        if (allocated(bead%cellv_jvdw)) deallocate(bead%cellv_jvdw)
        allocate(bead%cellv_jvdw(size(cellv_jvdw)))
      endif

      if(allocated(vcorrect_ik))  then
        if (allocated(bead%vcorrect_ik)) deallocate(bead%vcorrect_ik)
        allocate(bead%vcorrect_ik(size(vcorrect_ik,1)
     &        ,size(vcorrect_ik,2)))
      endif

      if(allocated(vcorrect_scale))  then
        if (allocated(bead%vcorrect_scale)) 
     &        deallocate(bead%vcorrect_scale)
        allocate(bead%vcorrect_scale(size(vcorrect_scale)))
      endif

      if(allocated(mcorrect_ik))  then
        if (allocated(bead%mcorrect_ik)) deallocate(bead%mcorrect_ik)
        allocate(bead%mcorrect_ik(size(mcorrect_ik,1)
     &        ,size(mcorrect_ik,2)))
      endif

      if(allocated(mcorrect_scale))  then
        if (allocated(bead%mcorrect_scale)) 
     &        deallocate(bead%mcorrect_scale)
        allocate(bead%mcorrect_scale(size(mcorrect_scale)))
      endif

      if(allocated(ccorrect_ik))  then
        if (allocated(bead%ccorrect_ik)) deallocate(bead%ccorrect_ik)
        allocate(bead%ccorrect_ik(size(ccorrect_ik,1)
     &        ,size(ccorrect_ik,2)))
      endif

      if(allocated(ccorrect_scale))  then
        if (allocated(bead%ccorrect_scale)) 
     &        deallocate(bead%ccorrect_scale)
        allocate(bead%ccorrect_scale(size(ccorrect_scale)))
      endif

      if(allocated(ucorrect_ik))  then
        if (allocated(bead%ucorrect_ik)) deallocate(bead%ucorrect_ik)
        allocate(bead%ucorrect_ik(size(ucorrect_ik)))
      endif

       if(allocated(ucorrect_scale))  then
        if (allocated(bead%ucorrect_scale)) 
     &        deallocate(bead%ucorrect_scale)
        allocate(bead%ucorrect_scale(size(ucorrect_scale)))
      endif

      if(allocated(dpcorrect_ik))  then
        if (allocated(bead%dpcorrect_ik)) deallocate(bead%dpcorrect_ik)
        allocate(bead%dpcorrect_ik(size(dpcorrect_ik)))
      endif

      if(allocated(dpcorrect_scale))  then
        if (allocated(bead%dpcorrect_scale)) 
     &        deallocate(bead%dpcorrect_scale)
        allocate(bead%dpcorrect_scale(size(dpcorrect_scale)))
      endif

      if(allocated(dpucorrect_ik))  then
        if (allocated(bead%dpucorrect_ik)) 
     &        deallocate(bead%dpucorrect_ik)
        allocate(bead%dpucorrect_ik(size(dpucorrect_ik)))
      endif

       if(allocated(dpucorrect_scale))  then
        if (allocated(bead%dpucorrect_scale)) 
     &        deallocate(bead%dpucorrect_scale)
        allocate(bead%dpucorrect_scale(size(dpucorrect_scale)))
      endif


      end subroutine resize_nl_arrays_bead
c
      subroutine savebeadnl(istep,bead)
      use domdec
      use neigh
      use potent
      use polpot
      use chgpot
      use mplpot
      use vdwpot
      implicit none
      TYPE(BEAD_TYPE), intent(inout) :: bead
      integer, intent(in) ::  istep
      integer modnl

!$acc wait

      modnl = mod(istep,ineigup)
      !write(0,*) "modnl=",modnl,istep,ineigup

      if(modnl==0) call resize_nl_arrays_bead(bead)
c
c      ! COPY ARRAYS THAT CHANGE EACH STEP
c
      if (allocated(celle_loc)) then
!$acc update host(celle_loc)
        bead%celle_loc = celle_loc
      endif

      if (allocated(celle_ploc)) then
!$acc update host(celle_ploc)
        bead%celle_ploc = celle_ploc
      endif

      if (allocated(celle_x)) then
!$acc update host(celle_x)
        bead%celle_x = celle_x
      endif

      if (allocated(celle_y)) then
!$acc update host(celle_y)
        bead%celle_y = celle_y
      endif

      if (allocated(celle_z)) then
!$acc update host(celle_z)
        bead%celle_z = celle_z
      endif


      if (modnl.ne.0) return
c
c      ! COPY ARRAYS THAT CHANGE ONLY WHEN WE RECOMPUTE THE NBLIST
c

      if (allocated(nvlst)) then
!$acc update host(nvlst)
        bead%nvlst = nvlst
      endif

      if(allocated(vlst)) then
!$acc update host(vlst)
        bead%vlst = vlst
      endif

      if(allocated(nelst)) then
!$acc update host(nelst)
        bead%nelst = nelst
      endif

      if (allocated(elst)) then        
!$acc update host(elst)
        bead%elst = elst
      end if

      if (allocated(nelstc)) then        
!$acc update host(nelstc)
        bead%nelstc = nelstc
      end if

      if (allocated(shortelst)) then        
!$acc update host(shortelst)
        bead%shortelst = shortelst
      end if

      if (allocated(nshortelst)) then        
!$acc update host(nshortelst)
        bead%nshortelst = nshortelst
      end if

      if (allocated(nshortelstc)) then        
!$acc update host(nshortelstc)
        bead%nshortelstc = nshortelstc
      end if

      if (allocated(eblst)) then        
!$acc update host(eblst)
        bead%eblst = eblst
      end if

       if (allocated(ieblst)) then        
!$acc update host(ieblst)
        bead%ieblst = ieblst
      end if

      if (allocated(shorteblst)) then        
!$acc update host(shorteblst)
        bead%shorteblst = shorteblst
      end if

      if (allocated(ishorteblst)) then        
!$acc update host(ishorteblst)
        bead%ishorteblst = ishorteblst
      end if
      

      if (allocated(nshortvlst)) then        
!$acc update host(nshortvlst)
        bead%nshortvlst = nshortvlst
      end if

      if (allocated(shortvlst)) then        
!$acc update host(shortvlst)
        bead%shortvlst = shortvlst
      end if

      if (allocated(vblst)) then        
!$acc update host(vblst)
        bead%vblst = vblst
      end if

      if (allocated(ivblst)) then        
!$acc update host(ivblst)
        bead%ivblst = ivblst
      end if

      if (allocated(shortvblst)) then        
!$acc update host(shortvblst)
        bead%shortvblst = shortvblst
      end if

      if (allocated(ishortvblst)) then        
!$acc update host(ishortvblst)
        bead%ishortvblst = ishortvblst
      end if

      if(allocated(celle_glob))  then
!$acc update host(celle_glob)
        bead%celle_glob = celle_glob
      endif

      if(allocated(celle_pole))  then
!$acc update host(celle_pole)
        bead%celle_pole = celle_pole
      endif

      if(allocated(celle_plocnl))  then
!$acc update host(celle_plocnl)
        bead%celle_plocnl = celle_plocnl
      endif
      
      if(allocated(celle_key))  then
!$acc update host(celle_key)
        bead%celle_key = celle_key
      endif

      if(allocated(celle_chg))  then
!$acc update host(celle_chg)
        bead%celle_chg = celle_chg
      endif

      if(allocated(cellv_key))  then
!$acc update host(cellv_key)
        bead%cellv_key = cellv_key
      endif

      if(allocated(cellv_glob))  then
!$acc update host(cellv_glob)
        bead%cellv_glob = cellv_glob
      endif

      if(allocated(cellv_loc))  then
!$acc update host(cellv_loc)
        bead%cellv_loc = cellv_loc
      endif

       if(allocated(cellv_jvdw))  then
!$acc update host(cellv_jvdw)
        bead%cellv_jvdw = cellv_jvdw
      endif

      ! SCALING FACTORS
      bead%n_vscale=n_vscale
      bead%n_mscale=n_mscale
      bead%n_cscale=n_cscale
      bead%n_uscale=n_uscale
      bead%n_dpscale=n_dpscale
      bead%n_dpuscale=n_dpuscale
      
      if(allocated(vcorrect_ik) .and.size(vcorrect_ik)>0) then
      !  write(0,*) "vcorrect_ik",size(vcorrect_ik)
!$acc update host(vcorrect_ik)
        bead%vcorrect_ik = vcorrect_ik
      endif

      if(allocated(vcorrect_scale).and.size(vcorrect_scale)>0) then
      !  write(0,*) "vcorrect_scale",size(vcorrect_scale)
!$acc update host(vcorrect_scale)
        bead%vcorrect_scale = vcorrect_scale
      endif

      if(allocated(mcorrect_ik).and.size(mcorrect_ik)>0) then
      !  write(0,*) "mcorrect_ik",size(mcorrect_ik)
!$acc update host(mcorrect_ik)
        bead%mcorrect_ik = mcorrect_ik
      endif

      if(allocated(mcorrect_scale).and.size(mcorrect_scale)>0) then
      !  write(0,*) "mcorrect_scale",size(mcorrect_scale)
!$acc update host(mcorrect_scale)
        bead%mcorrect_scale = mcorrect_scale
      endif

      if(allocated(ccorrect_ik).and.size(ccorrect_ik)>0) then
      !  write(0,*) "ccorrect_ik",size(ccorrect_ik)
!$acc update host(ccorrect_ik)
        bead%ccorrect_ik = ccorrect_ik
      endif

      if(allocated(ccorrect_scale).and.size(ccorrect_scale)>0) then
      !   write(0,*) "ccorrect_scale",size(ccorrect_scale)
!$acc update host(ccorrect_scale)
        bead%ccorrect_scale = ccorrect_scale
      endif

      if(allocated(ucorrect_ik).and. size(ucorrect_ik)>0) then
       !   write(0,*) "ucorrect_ik",size(ucorrect_ik)
!$acc update host(ucorrect_ik)
        bead%ucorrect_ik = ucorrect_ik
      endif

       if(allocated(ucorrect_scale).and.size(ucorrect_scale)>0) then
       !  write(0,*) "ucorrect_scale",size(ucorrect_scale)
!$acc update host(ucorrect_scale)
        bead%ucorrect_scale = ucorrect_scale
      endif

      if(allocated(dpcorrect_ik).and.size(dpcorrect_ik)>0)  then
      !  write(0,*) "dpcorrect_ik",size(dpcorrect_ik)
!$acc update host(dpcorrect_ik)
        bead%dpcorrect_ik = dpcorrect_ik
      endif

      if(allocated(dpcorrect_scale).and.size(dpcorrect_scale)>0) then
       !   write(0,*) "dpcorrect_scale",size(dpcorrect_scale)
!$acc update host(dpcorrect_scale)
        bead%dpcorrect_scale = dpcorrect_scale
      endif

      if(allocated(dpucorrect_ik).and.size(dpucorrect_ik)>0) then
      !  write(0,*) "dpucorrect_ik",size(dpucorrect_ik)
!$acc update host(dpucorrect_ik)
        bead%dpucorrect_ik = dpucorrect_ik
      endif

       if(allocated(dpucorrect_scale).and.size(dpucorrect_scale)>0) then
       ! write(0,*) "dpucorrect_scale",size(dpucorrect_scale)
!$acc update host(dpucorrect_scale)
        bead%dpucorrect_scale = dpucorrect_scale
      endif

!$acc wait
      end
      

      subroutine pushbead(istep,bead,skip_parameters)
      use angang
      use angle
      use atomsMirror
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
      use polpot
      use chgpot
      use mplpot
      use vdwpot
      use utilgpu,only:prmem_request
      implicit none
      integer, intent(in) :: istep
      type(BEAD_TYPE), intent(inout) :: bead
      LOGICAL, intent(in) :: skip_parameters
      integer modnl,n1,n2

!$acc wait
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
!$acc update device(x(:),y(:),z(:),v(:,:),a(:,:)) async

      nloc = bead%nloc
      glob = bead%glob
!$acc update device(glob(:)) async
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

!$acc update device(loc(:),ineignl(:)) async
!$acc update device(repart(:),repartrec(:)) async
!$acc update device(domlenrec(:),domlenpole(:),domlenpolerec(:)) async
!$acc update device(globrec(:),locrec(:),globrec1(:),locrec1(:)) async
!$acc update device(bufbegpole(:)) async
!$acc update device(buflen2(:),buf1(:),buf2(:)) async
!$acc update device(bufbeg1(:),bufbeg2(:)) async

c
c     VDW
c
      if (use_vdw) then
        ! write(0,*) "push VDW"
        nvdwbloc = bead%nvdwbloc
        vdwglob = bead%vdwglob
        vdwglobnl = bead%vdwglobnl
        vdwlocnl = bead%vdwlocnl
        nvdwlocnl = bead%nvdwlocnl
        nvdwblocloop = bead%nvdwblocloop
        nvdwlocnlb = bead%nvdwlocnlb
        nvdwlocnlb_pair = bead%nvdwlocnlb_pair
        nvdwlocnlb2_pair = bead%nvdwlocnlb2_pair
        nshortvdwlocnlb2_pair = bead%nshortvdwlocnlb2_pair
!$acc update device(vdwglob(:),vdwglobnl(:),vdwlocnl(:)) async
      endif
c
c     BOND
c
      if (use_bond) then
        !write(0,*) "push BOND"
        nbondloc = bead%nbondloc
        bndglob = bead%bndglob
!$acc update device(bndglob(:)) async
      endif
c
c     STRETCH-BEND
c
      if (use_strbnd) then
        !write(0,*) "push STRETCH-BEND"
        nstrbndloc = bead%nstrbndloc
        strbndglob = bead%strbndglob
!$acc update device(strbndglob(:)) async
      endif
c
c     UREY-BRADLEY
c
      if (use_urey) then
        !write(0,*) "push UREY-BRADLEY"
        nureyloc = bead%nureyloc
        ureyglob = bead%ureyglob
!$acc update device(ureyglob(:)) async
      endif
c
c     ANGLE-ANGLE
c
      if (use_angang) then
        !write(0,*) "push ANGLE-ANGLE"
        nangangloc = bead%nangangloc
        angangglob = bead%angangglob
!$acc update device(angangglob(:)) async
      endif
c
c     OP-BENDING
c
      if (use_opbend) then
        !write(0,*) "push OP-BENDING"
        nopbendloc = bead%nopbendloc
        opbendglob = bead%opbendglob
!$acc update device(opbendglob(:)) async
      endif
c
c     OP-DIST
c
      if (use_opdist) then
        !write(0,*) "push  OP-DIST"
        nopdistloc = bead%nopdistloc
        opdistglob = bead%opdistglob
!$acc update device(opdistglob(:)) async
      endif
c
c     IMPROP
c
      if (use_improp) then
       ! write(0,*) "push IMPROP"
        niproploc = bead%niproploc
        impropglob = bead%impropglob
!$acc update device(impropglob(:)) async
      endif
c
c     IMPTOR
c
      if (use_imptor) then
       ! write(0,*) "push IMPTOR"
        nitorsloc = bead%nitorsloc
        imptorglob = bead%imptorglob
!$acc update device(imptorglob(:)) async
      endif
c
c     TORSION
c
      if (use_tors) then
        !write(0,*) "push TORSION"
        ntorsloc = bead%ntorsloc
        torsglob = bead%torsglob
!$acc update device(torsglob(:)) async
      endif
c
c     PITORSION
c
      if (use_pitors) then
        !write(0,*) "push PITORSION"
        npitorsloc = bead%npitorsloc
        pitorsglob = bead%pitorsglob
!$acc update device(pitorsglob(:)) async
      endif
c
c     STRETCH-TORSION
c
      if (use_strtor) then
        !write(0,*) "push STRETCH-TORSION"
        nstrtorloc = bead%nstrtorloc
        strtorglob = bead%strtorglob
!$acc update device(strtorglob(:)) async
      endif
c
c     TORSION-TORSION
c
      if (use_tortor) then
        !write(0,*) "push TORSION-TORSION"
        ntortorloc = bead%ntortorloc
        tortorglob = bead%tortorglob
!$acc update device(tortorglob(:)) async
      endif
c
c     ANGLE
c
      if (use_angle) then
        !write(0,*) "push ANGLE"
        nangleloc = bead%nangleloc
        angleglob = bead%angleglob
!$acc update device(angleglob(:)) async
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
        nionlocloop = bead%nionlocloop
        nionlocnlloop=bead%nionlocnlloop
        nionlocnlb = bead%nionlocnlb
        nionlocnlb_pair = bead%nionlocnlb_pair
        nionlocnlb2_pair = bead%nionlocnlb2_pair
        nshortionlocnlb2_pair = bead%nshortionlocnlb2_pair
      
!$acc update device(chgrecglob(:),chgglob(:),chgglobnl(:)) async
      endif

c
c     MULTIPOLE
c
      if (use_mpole) then
        !write(0,*) "push MULTIPOLE"
        npolerecloc = bead%npolerecloc
        if(nlocrec>0) then
          polerecglob = bead%polerecglob(1:nlocrec)
!$acc update device(polerecglob(:)) async
        endif
        npoleloc = bead%npoleloc
        npolebloc = bead%npolebloc
        poleglob = bead%poleglob(1:nbloc)
        poleloc = bead%poleloc
        polelocnl = bead%polelocnl
        if(nlocnl>0) then
          poleglobnl = bead%poleglobnl(1:nlocnl)
!$acc update device(poleglobnl(:)) async
        endif
        npolelocnl = bead%npolelocnl
        npolelocnlb = bead%npolelocnlb
        npolelocnlb_pair = bead%npolelocnlb_pair
        npolelocnlb2_pair = bead%npolelocnlb2_pair
        nshortpolelocnlb2_pair = bead%nshortpolelocnlb2_pair
        npolerecloc_old = bead%npolerecloc_old
        npolelocloop = bead%npolelocloop
        npolelocnlloop = bead%npolelocnlloop
        npoleblocloop = bead%npoleblocloop
        npolereclocloop = bead%npolereclocloop

!$acc update device(poleglob(:),poleloc(:),polelocnl(:)) async
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
!$acc update device(udalt(:,:,:),upalt(:,:,:)) async
!$acc update device(uind(:,:),uinp(:,:)) async
      endif

      !NEIGHBORLIST

      if (allocated(celle_loc)) then
        n1=size(bead%celle_loc)
        call prmem_request(celle_loc,n1
     &        ,async=.true.)
        celle_loc(1:n1) = bead%celle_loc(1:n1)
!$acc update device(celle_loc) async
      endif

      if (allocated(celle_ploc)) then
        n1=size(bead%celle_ploc)
        call prmem_request(celle_ploc,n1
     &        ,async=.true.)
        celle_ploc(1:n1) = bead%celle_ploc(1:n1)
!$acc update device(celle_ploc) async
      endif

      if (allocated(celle_x)) then
        n1=size(bead%celle_x)
        call prmem_request(celle_x,n1
     &        ,async=.true.)
        celle_x(1:n1) = bead%celle_x(1:n1)
!$acc update device(celle_x) async
      endif

      if (allocated(celle_y)) then
        n1=size(bead%celle_y)
        call prmem_request(celle_y,n1
     &        ,async=.true.)
        celle_y(1:n1) = bead%celle_y(1:n1)
!$acc update device(celle_y) async
      endif

      if (allocated(celle_z)) then
        n1=size(bead%celle_z)
        call prmem_request(celle_z,n1
     &        ,async=.true.)
        celle_z(1:n1) = bead%celle_z(1:n1)
!$acc update device(celle_z) async
      endif

      if (allocated(nvlst)) then
        n1=size(bead%nvlst)
        call prmem_request(nvlst,n1
     &        ,async=.true.)
        nvlst(1:n1) =  bead%nvlst(1:n1)
!$acc update device(nvlst) async
      endif

      if(allocated(vlst)) then
        n1=size(bead%vlst,1)
        n2=size(bead%vlst,2)
        call prmem_request(vlst,n1,n2
     &        ,async=.true.)
        vlst(1:n1,1:n2) = bead%vlst(1:n1,1:n2) 
!$acc update device(vlst) async
      endif

      if(allocated(nelst)) then
        n1=size(bead%nelst)
        call prmem_request(nelst,n1
     &        ,async=.true.)
        bead%nelst(1:n1) = bead%nelst(1:n1)
!$acc update device(nelst) async
      endif

      if (allocated(elst)) then  
        n1=size(bead%elst,1)
        n2=size(bead%elst,2)
        call prmem_request(elst,n1,n2
     &        ,async=.true.)
        elst(1:n1,1:n2) = bead%elst(1:n1,1:n2) 
!$acc update device(elst) async
      end if

      if (allocated(nelstc)) then   
        n1=size(bead%nelstc)
        call prmem_request(nelstc,n1
     &        ,async=.true.)
        nelstc(1:n1) = bead%nelstc(1:n1)
!$acc update device(nelstc) async
      end if

      if (allocated(shortelst)) then 
        n1=size(bead%shortelst,1)
        n2=size(bead%shortelst,2)
        call prmem_request(shortelst,n1,n2
     &        ,async=.true.)
        shortelst(1:n1,1:n2) = bead%shortelst(1:n1,1:n2) 
!$acc update device(shortelst) async
      end if

      if (allocated(nshortelst)) then      
        n1=size(bead%nshortelst)
        call prmem_request(nshortelst,n1
     &        ,async=.true.)
        nshortelst(1:n1) = bead%nshortelst(1:n1)
!$acc update device(nshortelst) async
      end if

      if (allocated(nshortelstc)) then   
        n1=size(bead%nshortelstc)
        call prmem_request(nshortelstc,n1
     &        ,async=.true.)
        nshortelstc(1:n1) = bead%nshortelstc(1:n1)
!$acc update device(nshortelstc) async
      end if

      if (allocated(eblst)) then   
        n1=size(bead%eblst)
        call prmem_request(eblst,n1
     &        ,async=.true.)
        eblst(1:n1) = bead%eblst(1:n1)
!$acc update device(eblst) async
      end if

      if (allocated(ieblst)) then  
        n1=size(bead%ieblst)
        call prmem_request(ieblst,n1
     &        ,async=.true.)
        ieblst(1:n1) = bead%ieblst(1:n1)
!$acc update device(ieblst) async
      end if

      if (allocated(shorteblst)) then 
        n1=size(bead%shorteblst)
        call prmem_request(shorteblst,n1
     &        ,async=.true.)
        shorteblst(1:n1) = bead%shorteblst(1:n1)
!$acc update device(shorteblst) async
      end if

      if (allocated(ishorteblst)) then   
        n1=size(bead%ishorteblst)
        call prmem_request(ishorteblst,n1
     &        ,async=.true.)
        ishorteblst(1:n1) = bead%ishorteblst(1:n1)
!$acc update device(ishorteblst) async
      end if
      

      if (allocated(nshortvlst)) then 
        n1=size(bead%nshortvlst)
        call prmem_request(nshortvlst,n1
     &        ,async=.true.)
        nshortvlst(1:n1) = bead%nshortvlst(1:n1)
!$acc update device(nshortvlst) async
      end if

      if (allocated(shortvlst)) then   
        n1=size(bead%shortvlst,1)
        n2=size(bead%shortvlst,2)
        call prmem_request(shortvlst,n1,n2
     &        ,async=.true.)
        shortvlst(1:n1,1:n2) = bead%shortvlst(1:n1,1:n2) 
!$acc update device(shortvlst) async
      end if

      if (allocated(vblst)) then   
        n1=size(bead%vblst)
        call prmem_request(vblst,n1
     &        ,async=.true.)
        vblst(1:n1) = bead%vblst(1:n1)
!$acc update device(vblst) async
      end if

      if (allocated(ivblst)) then  
        n1=size(bead%ivblst)
        call prmem_request(ivblst,n1
     &        ,async=.true.)
        ivblst(1:n1) = bead%ivblst(1:n1)
!$acc update device(ivblst) async
      end if

      if (allocated(shortvblst)) then  
        n1=size(bead%shortvblst)
        call prmem_request(shortvblst,n1
     &        ,async=.true.)
        shortvblst(1:n1) = bead%shortvblst(1:n1)
!$acc update device(shortvblst) async
      end if

      if (allocated(ishortvblst)) then
        n1=size(bead%ishortvblst)
        call prmem_request(ishortvblst,n1
     &        ,async=.true.)
        ishortvblst(1:n1) = bead%ishortvblst(1:n1)
!$acc update device(ishortvblst) async
      end if

      if(allocated(celle_glob))  then
        n1=size(bead%celle_glob)
        call prmem_request(celle_glob,n1
     &        ,async=.true.)
        celle_glob(1:n1) = bead%celle_glob(1:n1)
!$acc update device(celle_glob) async
      endif

      if(allocated(celle_pole))  then
        n1=size(bead%celle_pole)
        call prmem_request(celle_pole,n1
     &        ,async=.true.)
        celle_pole(1:n1) = bead%celle_pole(1:n1)
!$acc update device(celle_pole) async
      endif

      if(allocated(celle_plocnl))  then
        n1=size(bead%celle_plocnl)
        call prmem_request(celle_plocnl,n1
     &        ,async=.true.)
        celle_plocnl(1:n1) = bead%celle_plocnl(1:n1)
!$acc update device(celle_plocnl) async
      endif
      
      if(allocated(celle_key))  then
        n1=size(bead%celle_key)
        call prmem_request(celle_key,n1
     &        ,async=.true.)
        celle_key(1:n1) = bead%celle_key(1:n1)
!$acc update device(celle_key) async
      endif

      if(allocated(celle_chg))  then
        n1=size(bead%celle_chg)
        call prmem_request(celle_chg,n1
     &        ,async=.true.)
        celle_chg(1:n1) = bead%celle_chg(1:n1)
!$acc update device(celle_chg) async
      endif

      if(allocated(cellv_key))  then
        n1=size(bead%cellv_key)
        call prmem_request(cellv_key,n1
     &        ,async=.true.)
        cellv_key(1:n1) = bead%cellv_key(1:n1) 
!$acc update device(cellv_key) async
      endif

      if(allocated(cellv_glob))  then
        n1=size(bead%cellv_glob)
        call prmem_request(cellv_glob,n1
     &        ,async=.true.)
        cellv_glob(1:n1) = bead%cellv_glob(1:n1) 
!$acc update device(cellv_glob) async
      endif

      if(allocated(cellv_loc))  then
        n1=size(bead%cellv_loc)
        call prmem_request(cellv_loc,n1
     &        ,async=.true.)
        cellv_loc(1:n1) = bead%cellv_loc(1:n1) 
!$acc update device(cellv_loc) async
      endif

       if(allocated(cellv_jvdw))  then
        n1=size(bead%cellv_jvdw)
        call prmem_request(cellv_jvdw,n1
     &        ,async=.true.)
        cellv_jvdw(1:n1) = bead%cellv_jvdw(1:n1)  
!$acc update device(cellv_jvdw) async
      endif

      !SCALING FACTORS
      n_vscale=bead%n_vscale
      n_mscale=bead%n_mscale
      n_cscale=bead%n_cscale
      n_uscale=bead%n_uscale
      n_dpscale=bead%n_dpscale
      n_dpuscale=bead%n_dpuscale

      if(allocated(vcorrect_ik).and.size(vcorrect_ik)>0) then
        n1=size(bead%vcorrect_ik,1)
        n2=size(bead%vcorrect_ik,2)
        call prmem_request(vcorrect_ik,n1,n2
     &        ,async=.true.)
        vcorrect_ik(1:n1,1:n2) = bead%vcorrect_ik(1:n1,1:n2) 
!$acc update device(vcorrect_ik) async
      endif

      if(allocated(vcorrect_scale).and.size(vcorrect_scale)>0) then
        n1=size(bead%vcorrect_scale)
        call prmem_request(vcorrect_scale,n1
     &        ,async=.true.)
        vcorrect_scale(1:n1) = bead%vcorrect_scale(1:n1)  
!$acc update device(vcorrect_scale) async
      endif

      if(allocated(mcorrect_ik).and.size(mcorrect_ik)>0) then
        n1=size(bead%mcorrect_ik,1)
        n2=size(bead%mcorrect_ik,2)
        call prmem_request(mcorrect_ik,n1,n2
     &        ,async=.true.)
        mcorrect_ik(1:n1,1:n2) = bead%mcorrect_ik(1:n1,1:n2) 
!$acc update device(mcorrect_ik) async
      endif

      if(allocated(mcorrect_scale).and.size(mcorrect_scale)>0) then
        n1=size(bead%mcorrect_scale)
        call prmem_request(mcorrect_scale,n1
     &        ,async=.true.)
        mcorrect_scale(1:n1) = bead%mcorrect_scale(1:n1)  
!$acc update device(mcorrect_scale) async
      endif

      if(allocated(ccorrect_ik).and.size(ccorrect_ik)>0) then
        n1=size(bead%ccorrect_ik,1)
        n2=size(bead%ccorrect_ik,2)
        call prmem_request(ccorrect_ik,n1,n2
     &        ,async=.true.)
        ccorrect_ik(1:n1,1:n2) = bead%ccorrect_ik(1:n1,1:n2) 
!$acc update device(ccorrect_ik) async
      endif

      if(allocated(ccorrect_scale).and.size(ccorrect_scale)>0) then
        n1=size(bead%ccorrect_scale)
        call prmem_request(ccorrect_scale,n1
     &        ,async=.true.)
        ccorrect_scale(1:n1) = bead%ccorrect_scale(1:n1)
!$acc update device(ccorrect_scale) async
      endif

      if(allocated(ucorrect_ik).and.size(ucorrect_ik)>0) then
        n1=size(bead%ucorrect_ik)
        call prmem_request(ucorrect_ik,n1
     &        ,async=.true.)
        ucorrect_ik(1:n1) = bead%ucorrect_ik(1:n1)
!$acc update device(ucorrect_ik) async
      endif

      if(allocated(ucorrect_scale).and.size(ucorrect_scale)>0) then
        n1=size(bead%ucorrect_scale)
        call prmem_request(ucorrect_scale,n1
     &        ,async=.true.)
        ucorrect_scale(1:n1) = bead%ucorrect_scale(1:n1)
!$acc update device(ucorrect_scale) async
      endif

      if(allocated(dpcorrect_ik).and.size(dpcorrect_ik)>0) then
        n1=size(bead%dpcorrect_ik)
        call prmem_request(dpcorrect_ik,n1
     &        ,async=.true.)
        dpcorrect_ik(1:n1) = bead%dpcorrect_ik(1:n1)
!$acc update device(dpcorrect_ik) async
      endif

      if(allocated(dpcorrect_scale).and.size(dpcorrect_scale)>0) then
        n1=size(bead%dpcorrect_scale)
        call prmem_request(dpcorrect_scale,n1
     &        ,async=.true.)
        dpcorrect_scale(1:n1) = bead%dpcorrect_scale(1:n1)
!$acc update device(dpcorrect_scale) async
      endif

      if(allocated(dpucorrect_ik).and.size(dpucorrect_ik)>0) then
        n1=size(bead%dpucorrect_ik)
        call prmem_request(dpucorrect_ik,n1
     &        ,async=.true.)
        dpucorrect_ik(1:n1) = bead%dpucorrect_ik(1:n1)
!$acc update device(dpucorrect_ik) async
      endif

      if(allocated(dpucorrect_scale).and.size(dpucorrect_scale)>0) then
        n1=size(bead%dpucorrect_scale)
        call prmem_request(dpucorrect_scale,n1
     &        ,async=.true.)
        dpucorrect_scale(1:n1) = bead%dpucorrect_scale(1:n1)
!$acc update device(dpucorrect_scale) async
      endif


c
c     TIME
c 
      timestep = bead%timestep

      ! write(0,*) "pushbead done"

      !SYNCHRONIZE GPU with CPU
!$acc wait

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

!$acc wait
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
            centroid(:,i)=centroid(:,i)+pos(:,i,ibead)
            vel_centroid(:,i)=vel_centroid(:,i)+vel(:,i,ibead)
          ENDDO
        ENDDO  
        centroid(:,:)=centroid(:,:)/REAL(nbeads)
        vel_centroid(:,:)=vel_centroid(:,:)/REAL(nbeads) 

        Ekcentroid=0.d0
        DO i=1,n ; DO j=1,3          
          Ekcentroid=Ekcentroid+mass(i)*vel_centroid(j,i)**2
        ENDDO ; ENDDO
        Ekcentroid=0.5d0*Ekcentroid/convert

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
        !ekprim = ekprim/nbeads/convert + eksumpi_loc

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
        ekvir=0.5d0*ekvir/nbeads

        presvir = prescon*( -dedv_mean + 2.d0*(Ekcentroid
     &               - ekvir)/(3.d0*volbox) )

        ekvir=0.5d0*nfree*boltzmann*kelvin/convert-ekvir
        !ekvir= Ekcentroid - ekvir
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

      !rebuild the neighbor lists
      if (use_list) call nblist(istep)

!$acc wait
          
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
