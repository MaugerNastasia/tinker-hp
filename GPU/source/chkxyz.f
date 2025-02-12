c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine chkxyz  --  check for coincident coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "chkxyz" finds any pairs of atoms with identical Cartesian
c     coordinates, and prints a warning message
c
c
#include "tinker_precision.h"
      subroutine chkxyz (clash)
      use sizes
      use atomsMirror
      use iounit
      implicit none
      integer i,j
      real(r_p) xi,yi,zi
      real(r_p) eps,r2
      logical clash
      logical header
c
c
c     initialize atom collision flag and distance tolerance
c
      clash = .false.
      eps   = 1d-6
c
c     loop over atom pairs testing for identical coordinates
c
      header = .true.
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do j = i+1, n
            r2 = (x(j)-xi)**2 + (y(j)-yi)**2 + (z(j)-zi)**2
            if (r2 .lt. eps) then
               clash = .true.
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format ()
               end if
               write (iout,20)  i,j
   20          format (' CHKXYZ  --  Warning, Atoms',i6,' and',i6,
     &                    ' have Identical Coordinates')
            end if
         end do
      end do
      end
