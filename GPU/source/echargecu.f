c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "echargecu" : driver for calculation of the point charge
c     energy and derivatives with respect to Cartesian coordinates on device
c
c
#ifndef TINKER_CUF
#define TINKER_CUF
#include "tinker_precision.h"
#include "tinker_types.h"
#include "tinker_cudart.h"
      module echargecu
        use utilcu  ,only: nproc,ndir,BLOCK_DIM,ALL_LANES,use_virial
        use utilgpu ,only: real3,real6,mdyn3_r,rpole_elt
     &              ,BLOCK_SIZE,RED_BUFF_SIZE,WARP_SIZE
        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "pair_charge.f.inc"

        attributes(global) subroutine ecreal1d_core_cu
     &        ( iion, cglob, loc, ieblst, eblst
     &        , nionlocnlb, nionlocnlb_pair, nionbloc, n
     &        , x, y, z, pchg
     &        , off2, f, aewald, ebuffer
     &        , dec, ec_buff, vir_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &        , inter
#endif
     &        )
        implicit none
        integer,value,intent(in)::nionlocnlb,nionbloc,n
     &         ,nionlocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,aewald,f,ebuffer
        integer,device,intent(in)::iion(*),cglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),pchg(*)
        real(t_p),device:: vir_buff(*)
        ener_rtyp,device:: ec_buff(*)
        mdyn_rtyp,device:: dec(3,*)
#ifdef TINKER_DEBUG
        integer  ,device:: inter(*)
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iichg,iglob,icloc,kchg,kcloc
#ifdef TINKER_DEBUG
        integer kdx_,istat
#endif
        integer location
        integer,shared,dimension(BLOCK_DIM)::kglob
        integer,parameter::no_scaling=0
        real(t_p) xk_,yk_,zk_,d2,fi
        ener_rtyp ec_
        real(t_p) rstat,zero
        type(real3) posi,pos
        type(real3) frc
        type(mdyn3_r) frc_i
        real(t_p)      ,shared::   fk(BLOCK_DIM)
        type(mdyn3_r)  ,shared::frc_k(BLOCK_DIM)
        type(real3)    ,shared:: posk(BLOCK_DIM)
        real(t_p) vir_(6)
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nionlocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kchg    = cglob(kdx)
           kglob(threadIdx%x)   = iion (kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           fk  (threadIdx%x)    = pchg(kchg)
           call syncwarp(ALL_LANES)

           !  Load atom block i parameters
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iichg   = cglob(idx)
           iglob   = iion (idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           fi      = f*pchg(iichg)

           ! zero data to compute
           frc_i%x = 0.0;
           frc_i%y = 0.0;
           frc_i%z = 0.0;
           frc_k(threadIdx%x)%x = 0.0;
           frc_k(threadIdx%x)%y = 0.0;
           frc_k(threadIdx%x)%z = 0.0;

           !* set compute Data to 0
           ec_ = 0
           vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
           vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_     = __shfl(kdx ,srclane)
#endif
              if (ndir.gt.1) then
                 xk_   = posk(klane)%x
                 yk_   = posk(klane)%y
                 zk_   = posk(klane)%z
                 pos%x = posi%x - xk_
                 pos%y = posi%y - yk_
                 pos%z = posi%z - zk_
                 call midpointimage_inl(xk_,yk_,zk_,pos%x,pos%y,pos%z)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 pos%x = posi%x - posk(klane)%x 
                 pos%y = posi%y - posk(klane)%y 
                 pos%z = posi%z - posk(klane)%z 
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                        ,same_block)

              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call charge_couple(d2,pos%x,pos%y,pos%z,ebuffer
     &                             ,fi*fk(klane),aewald,1.0_ti_p
     &                             ,ec_,frc,frc_k(klane),no_scaling)

                 vir_(1) = vir_(1) + pos%x * frc%x
                 vir_(2) = vir_(2) + pos%y * frc%x
                 vir_(3) = vir_(3) + pos%z * frc%x
                 vir_(4) = vir_(4) + pos%y * frc%y
                 vir_(5) = vir_(5) + pos%z * frc%y
                 vir_(6) = vir_(6) + pos%z * frc%z

                 frc_i%x = frc_i%x + tp2mdr( frc%x )
                 frc_i%y = frc_i%y + tp2mdr( frc%y )
                 frc_i%z = frc_i%z + tp2mdr( frc%z )

#ifdef TINKER_DEBUG
                 if (iglob<kglob(klane)) then
                    istat=AtomicAdd(inter(iglob),1)
                 else
                    istat=AtomicAdd(inter(kglob(klane)),1)
                 end if
#endif
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           ! Update energy buffer
           rstat = atomicAdd(ec_buff(location), ec_)

           ! Update virial buffer
           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))

           ! Update forces
           rstat = atomicAdd( dec(1,i   ), frc_i%x )
           rstat = atomicAdd( dec(2,i   ), frc_i%y )
           rstat = atomicAdd( dec(3,i   ), frc_i%z )
           call syncwarp(ALL_LANES)
           rstat = atomicAdd( dec(1,kbis), frc_k(threadIdx%x)%x )
           rstat = atomicAdd( dec(2,kbis), frc_k(threadIdx%x)%y )
           rstat = atomicAdd( dec(3,kbis), frc_k(threadIdx%x)%z )

        end do

        end subroutine

        attributes(global) subroutine ecreal3d_core_cu
     &        ( iion, cglob, loc, ieblst, eblst
     &        , nionlocnlb, nionlocnlb_pair, nionbloc, n
     &        , x, y, z, pchg
     &        , off2, f, aewald, ebuffer
     &        , ec_buff, nec_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &        )
        implicit none
        integer,value,intent(in)::nionlocnlb,nionbloc,n
     &         ,nionlocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,aewald,f,ebuffer
        integer,device,intent(in)::iion(*),cglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),pchg(*)
        integer  ,device:: nec_buff(*)
        ener_rtyp,device::  ec_buff(*)
#ifdef TINKER_DEBUG
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iichg,iglob,icloc,kchg,kglob,kcloc
#ifdef TINKER_DEBUG
        integer kglob_,kdx_,istat
#endif
        integer location
        integer,parameter::no_scaling=0
        integer nec_,istat
        real(t_p) xk_,yk_,zk_,d2,fi
        ener_rtyp ec_
        real(t_p) rstat,zero
        type(real3) posi,pos
        real(t_p)      ,shared::   fk(BLOCK_DIM)
        type(real3)    ,shared:: posk(BLOCK_DIM)
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nionlocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iichg   = cglob(idx)
           iglob   = iion (idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           fi      = f*pchg(iichg)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kchg    = cglob(kdx)
           kglob   = iion (kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x = x(kdx)
           posk(threadIdx%x)%y = y(kdx)
           posk(threadIdx%x)%z = z(kdx)
           fk  (threadIdx%x)   = pchg(kchg)

           !* set compute Data to 0
           ec_  = 0
           nec_ = 0

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_     = __shfl(kdx ,srclane)
              kglob_   = __shfl(kglob ,srclane)
#endif
              if (ndir.gt.1) then
                 xk_   = posk(klane)%x
                 yk_   = posk(klane)%y
                 zk_   = posk(klane)%z
                 pos%x = posi%x - xk_
                 pos%y = posi%y - yk_
                 pos%z = posi%z - zk_
                 call midpointimage_inl(xk_,yk_,zk_,pos%x,pos%y,pos%z)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 pos%x = posi%x - posk(klane)%x 
                 pos%y = posi%y - posk(klane)%y 
                 pos%z = posi%z - posk(klane)%z 
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)

              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call charge3_couple(d2,pos%x,pos%y,pos%z,ebuffer
     &                             ,fi*fk(klane),aewald,1.0_ti_p
     &                             ,ec_,no_scaling)
                 nec_ = nec_ + 1

#ifdef TINKER_DEBUG
#endif
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           ! Update energy buffer
           rstat = atomicAdd( ec_buff(location), ec_)
           istat = atomicAdd(nec_buff(location), nec_)
        end do

        end subroutine

        attributes(global) subroutine ecreal_scaling_cu
     &            ( ccorrect_ik,ccorrect_scale,loc,x,y,z
     &            , dec,ec_buff,vir_buff,n,nbloc,n_cscale
     &            , f,aewald,ebuffer,off2 )
        implicit none
        integer  ,value,intent(in)::n,nbloc,n_cscale
        real(t_p),value,intent(in):: f,aewald,ebuffer,off2
        integer  ,device,intent(in)::ccorrect_ik(n_cscale,2)
     &           , loc(n)
        real(t_p),device,intent(in):: ccorrect_scale(*)
     &           , x(n),y(n),z(n)
        real(t_p),device:: vir_buff(6*RED_BUFF_SIZE)
        ener_rtyp,device:: ec_buff(RED_BUFF_SIZE)
        mdyn_rtyp,device:: dec(3,nbloc)

        integer ii,iglob,kglob,i,k,ithread,lot
        real(t_p) rstat,r2
        real(t_p) xi,yi,zi
        real(t_p) xr,yr,zr
        real(t_p) scale_f,fik
        ener_rtyp e
        type(real3) ded
        type(mdyn3_r) dedc

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        do ii = ithread, n_cscale, blockDim%x*gridDim%x
           iglob   = ccorrect_ik(ii,1)
           kglob   = ccorrect_ik(ii,2)
           scale_f =   ccorrect_scale(2*ii+1)
           fik     = f*ccorrect_scale(2*ii+2)
           xi      = x(iglob)
           yi      = y(iglob)
           zi      = z(iglob)
c
c       compute the energy contribution for this interaction
c
           xr      = xi - x(kglob)
           yr      = yi - y(kglob)
           zr      = zi - z(kglob)
c
c       find energy for interactions within real space cutoff
c
           call image_inl (xr,yr,zr)
           r2 = xr*xr + yr*yr + zr*zr
           if (r2 .le. off2) then
              i  = loc(iglob)
              k  = loc(kglob)
              e  = 0
              dedc%x=0; dedc%y=0; dedc%z=0;
              ! correct pair
              call charge_couple(r2,xr,yr,zr,ebuffer
     &                          ,fik,aewald,scale_f
     &                          ,e,ded,dedc,1)

              lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
              ! increment the overall energy and derivative expressions
              rstat = atomicAdd( ec_buff(lot),e )
              rstat = atomicAdd( dec(1,i),-dedc%x )
              rstat = atomicAdd( dec(2,i),-dedc%y )
              rstat = atomicAdd( dec(3,i),-dedc%z )
              rstat = atomicAdd( dec(1,k),+dedc%x )
              rstat = atomicAdd( dec(2,k),+dedc%y )
              rstat = atomicAdd( dec(3,k),+dedc%z )
              ! increment the internal virial tensor components
              rstat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+lot),xr*ded%x)
              rstat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+lot),yr*ded%x)
              rstat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+lot),zr*ded%x)
              rstat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+lot),yr*ded%y)
              rstat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+lot),zr*ded%y)
              rstat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+lot),zr*ded%z)
           end if
        end do
        end subroutine

      end module
#endif
