c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  vdwpot.i  --  specifics of van der Waals functional form  ##
c     ##                                                            ##
c     ################################################################
c
c
c     abuck       value of "A" constant in Buckingham vdw potential
c     bbuck       value of "B" constant in Buckingham vdw potential
c     cbuck       value of "C" constant in Buckingham vdw potential
c     ghal        value of "gamma" in buffered 14-7 vdw potential
c     dhal        value of "delta" in buffered 14-7 vdw potential
c     v2scale     factor by which 1-2 vdw interactions are scaled
c     v3scale     factor by which 1-3 vdw interactions are scaled
c     v4scale     factor by which 1-4 vdw interactions are scaled
c     v5scale     factor by which 1-5 vdw interactions are scaled
c     igauss      coefficients of Gaussian fit to vdw potential
c     ngauss      number of Gaussians used in fit to vdw potential
c     use_vcorr   flag to use long range vdw der Waals correction
c     vdwindex    indexing mode (atom type or class) for vdw parameters
c     vdwtyp      type of van der Waals potential energy function
c     radtyp      type of parameter (sigma or R-min) for atomic size
c     radsiz      atomic size provided as radius or diameter
c     radrule     combining rule for atomic size parameters
c     epsrule     combining rule for vdw well depth parameters
c     gausstyp    type of Gaussian fit to van der Waals potential
c
c
      integer maxgauss
      parameter (maxgauss=10)
      integer ngauss
      real*8 abuck,bbuck,cbuck
      real*8 ghal,dhal
      real*8 v2scale,v3scale
      real*8 v4scale,v5scale
      real*8 igauss
      logical use_vcorr
      character*5 vdwindex
      character*5 radtyp
      character*8 radsiz,gausstyp
      character*10 radrule,epsrule
      character*13 vdwtyp
      common /vdwpot/ abuck,bbuck,cbuck,ghal,dhal,v2scale,v3scale,
     &                v4scale,v5scale,igauss(2,maxgauss),ngauss,
     &                use_vcorr,vdwindex,vdwtyp,radtyp,radsiz,radrule,
     &                epsrule,gausstyp
