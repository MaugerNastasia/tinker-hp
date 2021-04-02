c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine field  --  get the potential energy functions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "field" sets the force field potential energy functions from
c     a parameter file and modifications specified in a keyfile
c
c
#include "tinker_precision.h"
      subroutine field
      use keys
      use potent
      use uprior ,only: use_pred
      implicit none
      integer i
      character*240 record
c
c
c     set the default values for the active potentials
c
      use_bond       = .true.
      use_angle      = .true.
      use_strbnd     = .true.
      use_urey       = .true.
      use_angang     = .true.
      use_opbend     = .true.
      use_opdist     = .true.
      use_improp     = .true.
      use_imptor     = .true.
      use_tors       = .true.
      use_pitors     = .true.
      use_strtor     = .true.
      use_angtor     = .true.
      use_tortor     = .true.
      use_vdw        = .true.
      use_vdwshort   = .false.
      use_vdwlong    = .false.
      use_charge     = .true.
      use_clong      = .false.
      use_creal      = .true.
      use_crec       = .true.
      use_cself      = .true.
      use_mpole      = .true.
      use_mpoleshortreal = .false.
      use_mpolelong      = .false.
      use_mreal      = .true.
      use_mrec       = .true.
      use_mself      = .true.
      use_polar      = .true.
      use_preal      = .true.
      use_prec       = .true.
      use_pself      = .true.
      use_polarshortreal = .false.
      use_solv       = .false.
      use_geom       = .true.
      use_extra      = .false.
      use_emtp       = .false.
      use_ctransfer  = .false.
      use_dispersion = .false.
      use_repulsion  = .false.
c
c     Flag for predictor-Corrector
c
      use_pred       = .true.
c
c     Set default values of Force field potential
c
      PotentialAll         = .true.
      PotentialAmoeba      = .false.
      PotentialAmoeba18    = .false.
      PotentialWaterAmoeba = .false.
      PotentialWaterCharmm = .false.
      PotentialCharmm      = .false.
c
c     read the potential energy force field parameter file
c
      call getprm
c
c     check keywords for potential function control parameters
c
      do i = 1, nkey
         record = keyline(i)
         call prmkey (record)
      end do
!$acc update device(use_mpole,use_polar)
      end
