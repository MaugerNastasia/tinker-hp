c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module charge  --  partial charges for the current structure  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     nion      total number of partial charges in system
c     nionloc   local number of partial charges in system
c     nionbloc   local+neighbors number of partial charges in system
c     nionlocnl  localnl number of partial charges in system
c     nionrecloc   local reciprocal number of partial charges in system
c     iion      number of the atom site for each partial charge
c     winiion    window object corresponding to iion
c     jion      neighbor generation site for each partial charge
c     winjion    window object corresponding to jion
c     kion      cutoff switching site for each partial charge
c     winkion    window object corresponding to kion
c     chglist   partial charge site for each atom (0=no charge)
c     winchglist    window object corresponding to chglist
c     nbchg     number of charges before each index
c     winnbchg    window object corresponding to nbchg
c     chgloc    global-local charge correspondance
c     winchgloc    window object corresponding to chgloc
c     chglocnl  global-localnl charge correspondance
c     winchglocnl    window object corresponding to chglocnl
c     chgrecloc  global-local reciprocal charge correspondance
c     winchgrecloc    window object corresponding to chgrecloc
c     pchg      magnitude of the partial charges (e-)
c     winpchg    window object corresponding to pchg
c
c
      module charge
      implicit none
      integer nion,nionloc,nionbloc,nionlocnl,nionrecloc
      integer, allocatable :: chgloc(:),chglocnl(:)
      integer, allocatable :: chgrecloc(:)
      integer, pointer :: chglist(:)
      integer, pointer :: iion(:)
      integer, pointer :: jion(:),kion(:)
      integer, pointer :: nbchg(:)
      integer :: winiion,winjion,winkion,winchglist,winnbchg
      integer :: winchgloc,winchglocnl,winchgrecloc
      real*8, pointer ::  pchg(:)
      integer :: winpchg
      save 
      end
