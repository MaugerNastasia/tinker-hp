c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c   module stat : average values during a dynamic
c
c
#include "tinker_precision.h"
      module stat
      implicit none
      real(r_p) etot_sum,etot2_sum
      real(r_p) eint_sum,eint2_sum
      real(r_p) epot_sum,epot2_sum
      real(r_p) ekin_sum,ekin2_sum
      real(r_p) temp_sum,temp2_sum
      real(r_p) pres_sum,pres2_sum
      real(r_p) dens_sum,dens2_sum
      real(r_p) vol_sum,vol2_sum
      real(r_p) pistontemp_sum,pistontemp2_sum
      real(r_p) etotpi_sum,etot2pi_sum
      real(r_p) eintpi_sum,eint2pi_sum
      real(r_p) epotpi_sum,epot2pi_sum
      real(r_p) ekinpi_sum,ekin2pi_sum
      real(r_p) temppi_sum,temp2pi_sum
      real(r_p) prespi_sum,pres2pi_sum
      real(r_p) denspi_sum,dens2pi_sum
      real(r_p) volpi_sum,vol2pi_sum
      real(r_p) pistontemppi_sum,pistontemp2pi_sum
      save
      end
