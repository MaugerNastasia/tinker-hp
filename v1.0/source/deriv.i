c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  deriv.i  --  Cartesian coordinate derivative components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     desum   total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     deba    stretch-bend Cartesian coordinate derivatives
c     deub    Urey-Bradley Cartesian coordinate derivatives
c     deaa    angle-angle Cartesian coordinate derivatives
c     deopb   out-of-plane bend Cartesian coordinate derivatives
c     deopd   out-of-plane distance Cartesian coordinate derivatives
c     deid    improper dihedral Cartesian coordinate derivatives
c     deit    improper torsion Cartesian coordinate derivatives
c     det     torsional Cartesian coordinate derivatives
c     dept    pi-orbital torsion Cartesian coordinate derivatives
c     debt    stretch-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     dec     charge-charge Cartesian coordinate derivatives
c     decd    charge-dipole Cartesian coordinate derivatives
c     ded     dipole-dipole Cartesian coordinate derivatives
c     dem     multipole Cartesian coordinate derivatives
c     dep     polarization Cartesian coordinate derivatives
c     der     reaction field Cartesian coordinate derivatives
c     des     solvation Cartesian coordinate derivatives
c     delf    metal ligand field Cartesian coordinate derivatives
c     deg     geometric restraint Cartesian coordinate derivatives
c     dex     extra energy term Cartesian coordinate derivatives
c
c     dotstgrad : flag when the main program is testgrad (communication
c      of the forces one by one)
c
c
      real*8,pointer :: desum(:,:),deb(:,:),dea(:,:),deba(:,:)
      real*8,pointer :: deub(:,:),deaa(:,:),deopb(:,:),deopd(:,:)
      real*8,pointer :: deid(:,:),det(:,:),dept(:,:),deit(:,:)
      real*8,pointer :: debt(:,:),dett(:,:),dev(:,:)
      real*8,pointer :: ded(:,:),dem(:,:),dep(:,:),der(:,:)
      real*8,pointer :: des(:,:),delf(:,:),deg(:,:),dex(:,:)
      real*8, pointer :: demrec(:,:),deprec(:,:)
      real*8, pointer :: torquerec(:,:),torquedir(:,:)
      real*8, pointer :: torquerecp(:,:),torquedirp(:,:)
      real*8, pointer :: debond(:,:)
      logical dotstgrad
       common /deriv/ dem,dep,desum,
     $               dev,deb,dea,
     $               deba,deub,deaa,
     $               deopb,deopd,deid,
     $               deit,det,dept,
     $               debt,dett,deg,
     $               demrec,deprec,
     $               torquerec,torquedir,
     $               torquedirp,torquerecp,
     $               debond,dotstgrad
