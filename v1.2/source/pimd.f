c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  program pimd  --  run pimd molecular or stochastic dynamics  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "pimd" computes a molecular dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program pimd
      use mpi
      implicit none
      integer ierr,nthreadsupport
c      call MPI_INIT(ierr)
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
      call pimd_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end

      subroutine pimd_bis
      use atoms
      use bath
      use beads
      use bound
      use boxes
      use domdec
      use keys
      use inform
      use iounit
      use mdstuf
      use moldyn
      use neigh
      use timestat
      use mpi
      implicit none
      integer i,istep,nstep,ierr,k,ibead,ibeadglob!,nstep_therm
      integer iglob
      integer mode,next,nseg
      real*8 dt,dtdump,time0,time1
      real*8 timesteppimd
      logical exist,query,restart
      character*20 keyword
      character*120 record
      real*8 pres_tmp
      character*120 string
      real*8 maxwell
      logical isobaric_save
      integer nbeads_para
c
c
 1000 Format(' Time for 100 Steps: ',f15.4,/,
     $  ' Ave. Time per step: ',f15.4)
 1010 Format(' ns per day: ',f15.4)
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
c
c     check for keywords containing any altered parameters
c
      nbeads = 1 
      nbeads_ctr = 0
      nbeadsloc_ctr = 0
      nbeadsloc = 1 
      nproc = 1
!      nstep_therm=0
      do i = 1, nkey
        next = 1
        record = keyline(i)
        call gettext (record,keyword,next)
        call upcase (keyword)
        string = record(next:120)
        if (keyword(1:7) .eq. 'NBEADS ') then
          read (string,*,err=5,end=5) nbeads
        elseif (keyword(1:11) .eq. 'NBEADS_CTR ') then
            read (string,*,err=5,end=5) nbeads_ctr
        end if
   5  continue
      end do

      if(nbeads_ctr > 0 .AND. nbeads_ctr<nbeads) then
        contract = .TRUE. ; nbeads_para = nbeads_ctr
      elseif(nbeads_ctr>nbeads) then
        write(0,*) "nbeads_ctr must be lower than nbeads"
      else
        contract = .FALSE.; nbeads_para = nbeads
      endif

      if (nproctot.lt.nbeads_para) then
        nproc = 1
      else if (mod(nproctot,nbeads_para).ne.0) then
        if (ranktot.eq.0) then
          write(iout,*) 'inconsistent number of process for parallelism'
          write(iout,*) 'the total number of processors should be lower
     $     or a multiple of nbeads'
          call fatal
        end if
      else
        nproc = nproctot/nbeads_para
      end if

c
      call initmpi
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      nbeadsloc = int(nbeads/ncomm)
      if (rank_beadloc.eq.(ncomm-1)) nbeadsloc = 
     $    nbeads-(ncomm-1)*int(nbeads/ncomm)

      if(contract) then
        if (ranktot.eq.0 .AND. (mod(nbeads_ctr,2)==0)) then
          write(iout,*) 'nbeads_ctr must be an odd integer!'
          call fatal
        end if
        nbeadsloc_ctr = int(nbeads_ctr/ncomm)
        if (rank_beadloc.eq.(ncomm-1)) nbeadsloc_ctr = 
     $      nbeads_ctr-(ncomm-1)*int(nbeads_ctr/ncomm)
      endif
c      write(*,*) 'np = ',nproc,'nbeadsloc = ',nbeadsloc,'ranktot = ',
c     $ ranktot,ncomm,int(nbeads/ncomm)
c      if (ranktot.eq.0) write(*,*) 'bead_rank = ',bead_rank

      call cutoffs
      call unitcell
      call lattice
c
c     setup for MPI
c
c
c     get nprocforce to do correct allocation of MPI arrays
c
      call drivermpi
      call reinitnl(0)
c
c     allocate some arrays
c
      if (allocated(v)) deallocate (v)
      allocate (v(3,n))
      if (allocated(a)) deallocate (a)
      allocate (a(3,n))
      if (allocated(aalt)) deallocate (aalt)
      allocate (aalt(3,n))
c
      a = 0d0
      v = 0d0
      aalt = 0d0
c
      call mechanic
      call nblist(0)
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'BAOAB'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
c         else if (keyword(1:7) .eq. 'NBEADS ') then
c            read (string,*,err=5,end=5) nbeads
c         else if (keyword(1:11) .eq. 'NPROCBEADS ') then
c            read (string,*,err=5,end=5) nprocbeads
         end if
      end do
c      write(*,*) 'nbeads = ',nbeads
c      write(*,*) 'nbeadsloc = ',nbeadsloc
c      if (nprocbeads.gt.nbeads) then
c        write(iout,*) 'too many processes for beads parallelism'
c      else
c         nbeadsloctemp = nbeads/nprocbeads
c      end if
c
c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if

c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0d0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   70    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  110    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=120,end=120)  mode
  120    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,130)
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  200          continue
            end do
         end if
        if (mode.eq.1 .or. mode.eq.3) then
            if(ranktot.eq.0) then
               write(iout,*) 'NVE, NPH and NPT not implemented yet!'
               call fatal
            endif
        endif
        if(mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  240          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  mode
  250    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,260)
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  320          continue
            end do
         end if
      end if
c
c     setup dynamics
c
      call mdinit(dt)      
      call allocpi()      

      do ibead = 1, nbeadsloc
        call mdinitbead(beadsloc(ibead)%ibead_glob,dt,restart) 
        call dedvcalc()
        call initbead(0,beadsloc(ibead),.FALSE.)
        call resize_nl_arrays_bead(0,beadsloc(ibead))           
        call savebeadnl(0,beadsloc(ibead))          
      end do

      if(contract) then
        write(0,*) "entering initialize_pimd_contractions"
        call initialize_pimd_contractions()
      endif


      if(isobaric) then
        extvol = volbox
        extvolold = volbox
        masspiston=masspiston*nbeads
        if(ranktot.eq.0) then
          vextvol = maxwell(masspiston,kelvin)
          aextvol = 0d0
        endif
      endif      
      
      if (restart .and. (ranktot.eq.0)) then
        write(*,*) "Restarting from previous dynamics."
      endif
c
c     print out a header line for the dynamics computation
c

       if ((ranktot.eq.0).and.(contract.eqv..false.)) then
        write (iout,'(A,A,A,i3,A)') ' Path Integral Molecular Dynamics',
     &       ' Trajectory via BAOAB Algorithm',
     &       ' with ',nbeads,' beads'
        else if ((ranktot.eq.0).and.(contract.eqv..true.)) then
        write (iout,'(A,A,A,i3,A,i3,A)')'Contracted Path Integral',
     &       ' Molecular Dynamics Trajectory via BAOAB Algorithm',
     &       ' with ',nbeads,' beads and ',nbeads_ctr,'', 
     &       ' contracted beads'
        endif

c
c     integrate equations of motion to take a time step
c

      do istep = 1, nstep
        ! perform a step of PIMD integration
        call integrate_pimd(istep,dt)
      end do

c     perform any final tasks before program exit
c
      call final

      end subroutine

      subroutine integrate_pimd(istep,dt)
        use atoms
        use bath
        use beads
        use bound
        use boxes
        use domdec
        use keys
        use inform
        use iounit
        use mdstuf
        use moldyn
        use timestat
        use mpi
        use baoabpi
        implicit none
        integer, intent(in) :: istep
        real*8, intent(in) :: dt
        integer :: ibead,i,j,iglob
        real*8 time0,time1,time00,time01
        real*8 timesteppimd, timebeads,timepush,timereass,timeinit
        real*8 timededvcalc,timebaoab2,timeaoa,time_com
        real*8 pres_tmp
        real*8, allocatable :: derivs(:,:)
        TYPE(POLYMER_COMM_TYPE) :: polymer,polymer_ctr
        logical :: skip_parameters_copy

        skip_parameters_copy=(nbeadsloc+nbeadsloc_ctr).eq.1

        timesteppimd=0d0
        timepush=0d0
        timereass=0d0
        timeinit=0d0
        time_com=0.d0

        time1 = mpi_wtime()
c       GATHER polymer info at ranktot 0 (arrays are allocated inside the subroutine)
        !write(0,*) "ranktot",ranktot,"entering gather_polymer"
        call gather_polymer(polymer,beadsloc,.TRUE.,.TRUE.,.TRUE.)

        !!!!!!!!CONTRACTIONS!!!!!!!!!!!!!!!!!!!!!!!
        if(contract) then
          call gather_polymer(polymer_ctr,beadsloc_ctr
     &                            ,.FALSE.,.FALSE.,.FALSE.)
          if(ranktot.eq.0) then
            call contract_polymer(polymer,polymer_ctr) 
          endif
          call broadcast_polymer(polymer_ctr,beadsloc_ctr
     &                            ,.TRUE.,.FALSE.,.FALSE.)
          do ibead = 1, nbeadsloc_ctr
            call pushbead(istep,beadsloc_ctr(ibead)
     &                            ,.false.)
            call prepare_loaded_bead(istep)
            call resize_nl_arrays_bead(istep,beadsloc_ctr(ibead))
            call compute_grad_slow()
            call savebeadnl(istep,beadsloc_ctr(ibead))
            call initbead(istep,beadsloc_ctr(ibead)
     &                           ,.false.)
          enddo
          call gather_polymer(polymer_ctr,beadsloc_ctr
     &                          ,.FALSE.,.FALSE.,.TRUE.)
          if(ranktot.eq.0) then
            call project_forces_contract(polymer,polymer_ctr) 
            call apply_b_slow(polymer,dt) 
          endif

          call deallocate_polymer(polymer_ctr)
        endif

        !write(0,*) "ranktot",ranktot,"entering compute_observables"
c       gather and compute OBSERVABLES
        call compute_observables_pi(polymer%pos
     &                    ,polymer%vel,polymer%forces,istep,dt)
        if(allocated(polymer%forces)) deallocate(polymer%forces)   
        if(allocated(polymer%forces_slow))
     &      deallocate(polymer%forces_slow)       

c       aggregate statistical AVERAGES and PRINT system info 
        if(ranktot.eq.0) call mdstatpi(istep,dt)
    
        time0 = mpi_wtime() 
        time_com=time_com+time0-time1

        
c        PROPAGATE AOA
        !write(0,*) "ranktot",ranktot,"entering aoapi"
        call apply_aoa_pi(istep,dt,polymer%pos,polymer%vel)
        time1 = mpi_wtime()
        timeaoa=time1-time0

c       BROADCAST polymer info from ranktot 0 to everyone else
        call broadcast_polymer(polymer,beadsloc,.TRUE.,.TRUE.,.FALSE.)

        if(allocated(polymer%pos)) deallocate(polymer%pos)
        if(allocated(polymer%vel)) deallocate(polymer%vel)
        
        time0=mpi_wtime()
        time_com=time_com+time0-time1

        !write(0,*) "ranktot",ranktot,"entering baoabpi2"
        do ibead = 1, nbeadsloc
          time00=mpi_wtime()
          ! LOAD current bead
          call pushbead(istep,beadsloc(ibead),skip_parameters_copy)
          time01=mpi_wtime()
          timepush=timepush+time01-time00

          ! REASSIGN all positions in the loaded bead
          ! and reconstruction of neighborlist
          call prepare_loaded_bead(istep)
          call resize_nl_arrays_bead(istep,beadsloc(ibead))          

          ! APPLY B (compute gradient)
          time00=mpi_wtime()
          call apply_b_pi(istep,dt)
          time01=mpi_wtime()
          timebaoab2=timebaoab2+time01-time00
          time00=mpi_wtime()

          ! COMPUTE kinetic energy of loaded bead
          call kinetic (eksumpi_loc,ekinpi_loc,temppi)
          etotpi_loc = eksumpi_loc + epotpi_loc

          ! SAVE bead trajectory
          call mdsavebeads (istep,dt)

          ! COMPUTE VIRIAL dE/dV
          !write(0,*) "ranktot",ranktot,"entering dedvcalc"
          call dedvcalc()
          time01=mpi_wtime()
          timededvcalc=timededvcalc+time01-time00
          
                   
          time00=mpi_wtime()
          ! STORE new neighborlist if necessary
          !write(0,*) "ranktot",ranktot,"entering savebeadnl" 
          if(.not.skip_parameters_copy) then
            call savebeadnl(istep,beadsloc(ibead))
          endif

          ! SAVE current bead
          call initbead(istep,beadsloc(ibead),skip_parameters_copy)
          time01=mpi_wtime()
          timeinit=timeinit+time01-time00

        end do

        call deallocate_polymer(polymer)

      end subroutine
c
      subroutine initialize_pimd_contractions()
        use atoms
        use bath
        use beads
        use bound
        use boxes
        use domdec
        use keys
        use inform
        use iounit
        use mdstuf
        use moldyn
        use neigh
        use timestat
        use mpi
        implicit none
        type(POLYMER_COMM_TYPE) :: polymer,polymer_ctr
        integer :: i,j,k,ibead,iproc,nbeadscomm,tagmpi
        integer :: ierr
        integer status(MPI_STATUS_SIZE)
        real*8, allocatable :: buffer(:,:,:)
        integer, allocatable :: reqrec(:),reqsend(:)

        allocate(reqsend(nproctot))
        allocate(reqrec(nproctot))

        call gather_polymer(polymer,beadsloc,.TRUE.,.FALSE.,.FALSE.)
        call get_polymer_info(polymer_ctr,beadsloc_ctr,.FALSE.,.TRUE.)
        if(ranktot.eq.0) then          
          call contract_polymer(polymer,polymer_ctr)
          call deallocate_polymer(polymer)
          do k=1,nbeadsloc_ctr
            ibead = beadsloc_ctr(k)%ibead_glob
            beadsloc_ctr(k)%x = polymer_ctr%pos(1,:,ibead)
            beadsloc_ctr(k)%y = polymer_ctr%pos(2,:,ibead)
            beadsloc_ctr(k)%z = polymer_ctr%pos(3,:,ibead)           
          enddo

          do iproc=1,nproctot-1
            ! SEND contracted beads
            tagmpi = iproc+1
            nbeadscomm = polymer_ctr%nbeadscomm(iproc+1)
            write(0,*) nbeadscomm
            allocate(buffer(3,n,nbeadscomm))
            do k=1,nbeadscomm
              ibead = polymer_ctr%globbeadcomm(k,iproc+1)
              buffer(:,:,k) = polymer_ctr%pos(:,:,ibead)
            enddo
            call MPI_ISEND(buffer,3*n*nbeadscomm,MPI_REAL8
     $        ,iproc,tagmpi,MPI_COMM_WORLD,reqsend(iproc),ierr)
            call MPI_WAIT(reqsend(iproc),status,ierr)
            deallocate(buffer)
          enddo

          call deallocate_polymer(polymer_ctr)
        else
          !RECEIVE contracted beads
          tagmpi=ranktot+1
          allocate(buffer(3,n,nbeadsloc_ctr))
          call MPI_IRECV(buffer,3*n*nbeadsloc_ctr,MPI_REAL8
     $      ,0,tagmpi,MPI_COMM_WORLD,reqrec(1),ierr)
          call MPI_WAIT(reqrec(1),status,ierr)
          do k=1,nbeadsloc_ctr
            beadsloc_ctr(k)%x=buffer(1,:,k)
            beadsloc_ctr(k)%y=buffer(2,:,k)
            beadsloc_ctr(k)%z=buffer(3,:,k)
          enddo
        endif

        
        do k=1,nbeadsloc_ctr
          call resize_nl_arrays_bead(0,beadsloc_ctr(k))   
          call pushbead(0,beadsloc_ctr(k),.FALSE.)
          call ddpme3d
          call reinitnl(0)
          call reassignpme(.false.)
          call mechanicstep(0)
          call nblist(0)
          call initbead(0,beadsloc_ctr(k),.FALSE.)
          call resize_nl_arrays_bead(0,beadsloc_ctr(k))   
          call savebeadnl(0,beadsloc_ctr(k))
        enddo        
        
        
        if(ranktot.eq.0) then
          write(*,*) "initialize_pimd_contractions done."
        endif

      end subroutine initialize_pimd_contractions





      subroutine contract_polymer(polymer,polymer_ctr)
        use atoms
        use atmtyp
        use bath
        use beads
        use bound
        use boxes
        use domdec
        use keys
        use inform
        use iounit
        use mdstuf
        use moldyn
        use neigh
        use timestat
        use mpi
        implicit none
        type(POLYMER_COMM_TYPE), intent(inout) :: polymer,polymer_ctr
        integer :: i,j,k,ibead,iproc,nbeadscomm,tagmpi

        if(allocated(polymer_ctr%pos)) 
     &               deallocate(polymer_ctr%pos)
        allocate(polymer_ctr%pos(3,n,nbeads_ctr))

        DO i=1,n ; DO j=1,3
          polymer_ctr%pos(j,i,:)=matmul(contractor_mat
     &                    ,polymer%pos(j,i,:))
        ENDDO ; ENDDO        

      end subroutine contract_polymer

      subroutine project_forces_contract(polymer,polymer_ctr)
        use atoms
        use atmtyp
        use bath
        use beads
        use bound
        use boxes
        use domdec
        use keys
        use inform
        use iounit
        use mdstuf
        use moldyn
        use neigh
        use timestat
        use mpi
        implicit none
        type(POLYMER_COMM_TYPE) :: polymer,polymer_ctr
        real*8, allocatable :: forces_full(:,:,:)
        real*8, allocatable :: eigmat_small(:,:)
        integer :: i,j,k,ibead,iproc,nbeadscomm,tagmpi

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
        


      
