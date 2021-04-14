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
#include "tinker_precision.h"
      program pimd
      use mpi
#ifdef _OPENACC
      use utilgpu,only: bind_gpu
#endif
      implicit none
      integer ierr,nthreadsupport
#ifdef _OPENACC
      call bind_gpu
#endif
      call MPI_INIT(ierr)
c      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
      call pimd_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end
c
      subroutine pimd_bis
      use atoms
      use bath
      use beads
      use bound
      use boxes, only: volbox
      use domdec
      use keys
      use inform
      use iounit
      use mdstuf
      use moldyn
      use neigh
      use timestat
      use mpi
      use utilgpu ,only: rec_queue
      use tinMemory
      use utils, only: set_to_zero2m
      implicit none
      integer i,istep,nstep,ierr,k,ibead,ibeadglob!,nstep_therm
      integer iglob
      integer mode,next,nseg
      real(r_p) dt,dtdump,time0,time1
      logical exist,query,restart
      character*20 keyword
      character*240 record
      character*240 string
      real(r_p) maxwell
      integer nbeads_para
c
c
 1000 Format(' Time for ',I4,' Steps: ',f15.4,/,
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
        string = record(next:240)
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

!$acc update device(nbeads)

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
      call prmem_requestm(    v,3,n)
      call prmem_requestm(    a,3,n)
      call prmem_requestm( aalt,3,n)
      call prmem_requestm(aalt2,3,n)
c
      call set_to_zero2m(   a,    v,3*n,rec_queue)
      call set_to_zero2m(aalt,aalt2,3*n,rec_queue)
c
      call mechanic
      call nblist(0)
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0
      atmsph = 0.0
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
         string = record(next:240)
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
      dt = -1.0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0)  dt = 1.0
   70    continue
      end do
      dt = 0.001 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0)  dtdump = 0.1
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
            kelvin = -1.0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0)  kelvin = 298.0
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
            atmsph = -1.0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0)  atmsph = 1.0
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
            kelvin = -1.0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0)  kelvin = 298.0
  320          continue
            end do
         end if
      end if
c
c     setup dynamics
c

!$acc enter data create(epotpi_loc,eksumpi_loc,etotpi_loc,ekinpi_loc)
!$acc&           create(eintrapi_loc,einterpi_loc,temppi)

      call mdinit(dt)      
      call allocpi()      

      do ibead = 1, nbeadsloc
        call mdinitbead(beadsloc(ibead)%ibead_glob,dt,restart) 
        call dedvcalc()
        call initbead(0,beadsloc(ibead),.FALSE.)
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
          aextvol = 0
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

      time0 = 0 !mpi_wtime()
      do istep = 1, nstep
        call timer_enter( timer_timestep )
        ! perform a step of PIMD integration
        call integrate_pimd(istep,dt)
        call timer_exit( timer_timestep )

        ! perform a step of PIMD integration
        if (mod(istep,iprint).eq.0) then
          time1 = timer_get_total( timer_timestep )
          timestep = time1-time0
          time0 = time1
          if (ranktot.eq.0) then
            call MPI_REDUCE(MPI_IN_PLACE,timestep,1,MPI_REAL8,MPI_SUM,
     $      0,MPI_COMM_WORLD,ierr)
          else
            call MPI_REDUCE(timestep,timestep,1,MPI_REAL8,MPI_SUM,0,
     $      MPI_COMM_WORLD,ierr)
          end if
          if (verbose.and.ranktot.eq.0) then
            write(6,1000) iprint, (timestep)/nproctot,(timestep)/
     $       dble(nproctot*iprint)
            write(6,1010) 86400*dt*real(iprint*nproctot,t_p)/
     $                                   (1000*timestep)
          end if
        end if

        ! Abort if any problem detected
        if (abort) call fatal
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
        use boxes, only: volbox
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
        real(r_p), intent(in) :: dt
        integer :: ibead,i,j,iglob
        real(8) time0,time1,time00,time01
        real(8) timebeads,timepush,timereass,timeinit
        real(8) timededvcalc,timebaoab2,timeaoa,time_com
        TYPE(POLYMER_COMM_TYPE) :: polymer,polymer_ctr
        logical :: skip_parameters_copy

        skip_parameters_copy=(nbeadsloc+nbeadsloc_ctr).eq.1

        timepush=0d0
        timereass=0d0
        timeinit=0d0
        time_com=0.d0

        time1 = mpi_wtime()
c       GATHER polymer info at ranktot 0 (arrays are allocated inside the subroutine)
        !write(0,*) "ranktot",ranktot,"entering gather_polymer"
        call gather_polymer(polymer,beadsloc,.TRUE.,.TRUE.,.TRUE.)

        ! CONTRACTIONS
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
!$acc data present(epotpi_loc,eksumpi_loc,etotpi_loc,ekinpi_loc)
!$acc&     present(eintrapi_loc,temppi,x,y,z,a)
            call prepare_loaded_bead(istep)            
            call compute_grad_slow()
!$acc end data
            call initbead(istep,beadsloc_ctr(ibead)
     &                           ,.false.)
            call savebeadnl(istep,beadsloc_ctr(ibead))
            
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
          !write(0,*) "start of bead",ibead
          time00=mpi_wtime()
          ! LOAD current bead
          call pushbead(istep,beadsloc(ibead)
     &          ,skip_parameters_copy)
          time01=mpi_wtime()
          timepush=timepush+time01-time00

          
          !write(0,*) "before prepare"
!$acc data present(epotpi_loc,eksumpi_loc,etotpi_loc,ekinpi_loc)
!$acc&     present(eintrapi_loc,x,y,z,v,a)
          ! REASSIGN all positions in the loaded bead
          ! and reconstruction of neighborlist
          call prepare_loaded_bead(istep)
                 

          ! APPLY B (compute gradient)
          !write(0,*) "before apply_b_pi"
          time00=mpi_wtime()
          call apply_b_pi(istep,dt)
          time01=mpi_wtime()
          timebaoab2=timebaoab2+time01-time00
          time00=mpi_wtime()

          !write(0,*) "before kinetic"
          ! COMPUTE kinetic energy of loaded bead
          call kineticgpu (eksumpi_loc,ekinpi_loc,temppi)
!$acc wait
          ! COMPUTE VIRIAL dE/dV

          !write(0,*) "ranktot",ranktot,"entering dedvcalc"
          !write(0,*) "before dedvcalc"
          call dedvcalc()
          time01=mpi_wtime()
          timededvcalc=timededvcalc+time01-time00  

          !write(0,*) "before mdsavebeads"
          ! SAVE bead trajectory
          call mdsavebeads (istep,dt)    

!$acc end data              

          !write(0,*) "before initbead"
          ! SAVE current bead
          call initbead(istep,beadsloc(ibead)
     &          ,skip_parameters_copy)
          time01=mpi_wtime()
          timeinit=timeinit+time01-time00

          time00=mpi_wtime()
          ! STORE new neighborlist if necessary
          !write(0,*) "before savebeadnl"
          if(.not.skip_parameters_copy) then
            call savebeadnl(istep,beadsloc(ibead))
          endif


        end do


        call deallocate_polymer(polymer)

        ! STEP DONE

        !write(0,*) "ranktot",ranktot,"step done"
c        if(ranktot.eq.0) then
c            write(*,*) 'Time pushbeads etc ', timebeads
c            write(*,*) 'Time pushbead ', timepush
c            write(*,*) 'Time reassigpi etc ', timereass
c            write(*,*) 'Time initbead etc ', timeinit
c            write(*,*) 'Time baoab 2 (gradient)  ', timebaoab2
c            write(*,*) 'Time dedv ', timededvcalc
c            write(*,*) 'Time aoa ', timeaoa
c      endif



      end subroutine

      subroutine initialize_pimd_contractions()
        use atomsMirror
        use atmtyp
        use usage
        use bath
        use beads
        use bound
        use boxes, only: volbox
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
        real(t_p), allocatable :: buffer(:,:,:)
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
            allocate(buffer(3,n,nbeadscomm))
            do k=1,nbeadscomm
              ibead = polymer_ctr%globbeadcomm(k,iproc+1)
              buffer(:,:,k) = polymer_ctr%pos(:,:,ibead)
            enddo
            call MPI_ISEND(buffer,3*n*nbeadscomm,MPI_TPREC
     $        ,iproc,tagmpi,MPI_COMM_WORLD,reqsend(iproc),ierr)
            call MPI_WAIT(reqsend(iproc),status,ierr)
            deallocate(buffer)
          enddo

          call deallocate_polymer(polymer_ctr)
        else
          !RECEIVE contracted beads
          tagmpi=ranktot+1
          allocate(buffer(3,n,nbeadsloc_ctr))
          call MPI_IRECV(buffer,3*n*nbeadsloc_ctr,MPI_TPREC
     $      ,0,tagmpi,MPI_COMM_WORLD,reqrec(1),ierr)
          call MPI_WAIT(reqrec(1),status,ierr)
          do k=1,nbeadsloc_ctr
            beadsloc_ctr(k)%x=buffer(1,:,k)
            beadsloc_ctr(k)%y=buffer(2,:,k)
            beadsloc_ctr(k)%z=buffer(3,:,k)
          enddo
        endif

        
        do k=1,nbeadsloc_ctr
          call resize_nl_arrays_bead(beadsloc_ctr(k))   
          call pushbead(0,beadsloc_ctr(k),.FALSE.)
!$acc data  present(x,y,z,xold,yold,zold,v,a,mass,glob,use)
          call ddpme3d
          call reinitnl(0)
          call reassignpme(.FALSE.)
          call mechanicstep(0)
          call nblist(0)
!$acc end data
          call initbead(0,beadsloc_ctr(k),.FALSE.)
          call savebeadnl(0,beadsloc_ctr(k))
        enddo        
        
        
        if(ranktot.eq.0) then
          write(*,*) "initialize_pimd_contractions done."
        !  write(0,*) "Error: PIMD contractions not implemented yet !"
        !  call fatal
        endif

      end subroutine initialize_pimd_contractions
        


      
