!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_micro_interface

  use mod_realkinds
  use mod_service
  use mod_stdio
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_regcm_types
  use mod_runparams
  use mod_micro_nogtom
  use mod_micro_subex
  use mod_micro_wsm5
  use mod_cloud_subex
  use mod_cloud_xuran
  use mod_cloud_thomp
  use mod_cloud_guli2007
  use mod_cloud_texeira
  use mod_cloud_tompkins
  use mod_cloud_echam5

  implicit none

  private

  public :: allocate_micro , init_micro , microscheme , cldfrac , condtq

  type(nogtom_stats) , public :: ngs

  type(mod_2_micro) :: mo2mc
  type(micro_2_mod) :: mc2mo

  real(rkx) , pointer , dimension(:,:) :: rh0
  ! rh0adj - Adjusted relative humidity threshold
  real(rkx) , pointer , dimension(:,:,:) :: totc , rh0adj

  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_qcn
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_qin
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_qrn
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_qsn
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_phs
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_qvn
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_qs
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_rh
  real(rkx) , pointer , dimension(:,:,:) :: mc2mo_fcc
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_t
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_rho
  integer(ik4) , pointer , dimension(:,:) :: mo2mc_ldmsk
  real(rkx) , pointer , dimension(:,:,:) :: mo2mc_z
  real(rkx) , pointer , dimension(:,:) :: mo2mc_ps2
  integer(ik4) , pointer , dimension(:,:) :: mo2mc_iveg
  real(rkx) , pointer , dimension(:,:) :: atms_th700
  real(rkx) , pointer , dimension(:,:,:) :: atms_th3d

  real(rkx) , parameter :: alphaice = d_one

  integer(ik4) , parameter :: nchi = 256
  real(rkx) , dimension(0:nchi-1) :: chis

  public :: qck1 , cgul , rh0 , cevap , xcevap , caccr

  contains

  subroutine allocate_micro
    implicit none
    integer(ik4) :: i
    real(rkx) :: cf
    if ( ipptls == 1 ) then
      call allocate_subex
      call getmem3d(rh0adj,jci1,jci2,ici1,ici2,1,kz,'micro:rh0adj')
!$acc enter data create(rh0adj)
    else if ( ipptls == 2 ) then
      call allocate_mod_nogtom
#ifdef DEBUG
      if ( stats ) then
        call getmem3d(ngs%statssupw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statssupw')
        call getmem3d(ngs%statssupc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statssupc')
        call getmem3d(ngs%statsdetrw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsdetrw')
        call getmem3d(ngs%statsdetrc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsdetrc')
        call getmem3d(ngs%statserosw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statserosw')
        call getmem3d(ngs%statserosc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statserosc')
        call getmem3d(ngs%statsevapw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsevapw')
        call getmem3d(ngs%statsevapc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsevapc')
        call getmem3d(ngs%statscond1w,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statscond1w')
        call getmem3d(ngs%statscond1c,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statscond1c')
        call getmem3d(ngs%statsdepos,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsdepos')
        call getmem3d(ngs%statsmelt,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsmelt')
        call getmem3d(ngs%statsfrz,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsfrz')
        call getmem3d(ngs%statsrainev,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsrainev')
        call getmem3d(ngs%statssnowev,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statssnowev')
        call getmem3d(ngs%statsautocvw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsautocvw')
        call getmem3d(ngs%statsautocvc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsautocvc')
      end if
#endif
    else if ( ipptls == 3 ) then
      call allocate_mod_wsm5
    end if
    call getmem2d(rh0,jci1,jci2,ici1,ici2,'subex:rh0')
    call getmem3d(totc,jci1,jci2,ici1,ici2,1,kz,'subex:totc')
!$acc enter data create(totc)
!$acc enter data create(chis)
!$acc parallel present(chis)
!$acc loop
    do i = 1 , nchi
      cf = real(i-1,rkx)/real(nchi-1,rkx)
      chis(i-1) = 0.97_rkx*exp(-((cf-0.098_rkx)**2)/0.0365_rkx)+0.255_rkx
    end do
!$acc end parallel
  end subroutine allocate_micro

  subroutine init_micro
    use mod_atm_interface
    use mod_che_interface
    implicit none

    call assignpnt(atms%th700,atms_th700)
!$acc enter data create(atms_th700)
    call assignpnt(atms%th3d,atms_th3d)
!$acc enter data create(atms_th3d)
    call assignpnt(mddom%ldmsk,mo2mc%ldmsk)
    call assignpnt(mo2mc%ldmsk,mo2mc_ldmsk)
!$acc enter data create(mo2mc_ldmsk)
    call assignpnt(mddom%iveg,mo2mc%iveg)
    call assignpnt(mo2mc%iveg,mo2mc_iveg)
!$acc enter data create(mo2mc_iveg)
    call assignpnt(mddom%xlat,mo2mc%xlat)
    call assignpnt(sfs%psb,mo2mc%psb)
    call assignpnt(atms%pb3d,mo2mc%phs)
    call assignpnt(mo2mc%phs,mo2mc_phs)
!$acc enter data create(mo2mc_phs)
    call assignpnt(atms%pf3d,mo2mc%pfs)
    call assignpnt(atms%tb3d,mo2mc%t)
    call assignpnt(mo2mc%t,mo2mc_t)
!$acc enter data create(mo2mc_t)
    call assignpnt(atms%za,mo2mc%z)
    call assignpnt(mo2mc%z,mo2mc_z)
!$acc enter data create(mo2mc_z)
    call assignpnt(atms%dzq,mo2mc%delz)
    call assignpnt(atms%wpx3d,mo2mc%pverv)
    call assignpnt(atms%wb3d,mo2mc%verv)
    call assignpnt(atms%qxb3d,mo2mc%qxx)
    call assignpnt(atms%rhob3d,mo2mc%rho)
    call assignpnt(mo2mc%rho,mo2mc_rho)
!$acc enter data create(mo2mc_rho)
    call assignpnt(atms%rhb3d,mo2mc%rh)
    call assignpnt(mo2mc%rh,mo2mc_rh)
!$acc enter data create(mo2mc_rh)
    call assignpnt(atms%qsb3d,mo2mc%qs)
    call assignpnt(mo2mc%qs,mo2mc_qs)
!$acc enter data create(mo2mc_qs)
    call assignpnt(atms%ps2d,mo2mc%ps2)
    call assignpnt(mo2mc%ps2,mo2mc_ps2)
!$acc enter data create(mo2mc_ps2)
    call assignpnt(heatrt,mo2mc%heatrt)
    call assignpnt(q_detr,mo2mc%qdetr)
    call assignpnt(cldfra,mo2mc%cldf)

    call assignpnt(atms%qxb3d,mo2mc%qvn,iqv)
    call assignpnt(mo2mc%qvn,mo2mc_qvn)
!$acc enter data create(mo2mc_qvn)
    call assignpnt(atms%qxb3d,mo2mc%qcn,iqc)
    call assignpnt(mo2mc%qcn,mo2mc_qcn)
!$acc enter data create(mo2mc_qcn)
    if ( ipptls > 1 ) then
      call assignpnt(atms%qxb3d,mo2mc%qin,iqi)
      call assignpnt(mo2mc%qin,mo2mc_qin)
!$acc enter data create(mo2mc_qin)
      call assignpnt(atms%qxb3d,mo2mc%qsn,iqs)
      call assignpnt(mo2mc%qsn,mo2mc_qsn)
!$acc enter data create(mo2mc_qsn)
      call assignpnt(atms%qxb3d,mo2mc%qrn,iqr)
      call assignpnt(mo2mc%qrn,mo2mc_qrn)
!$acc enter data create(mo2mc_qrn)
    end if

    if ( ichem == 1 ) then
      if ( iaerosol == 1 .and. iindirect == 2 ) then
        call assignpnt(ccn,mo2mc%ccn)
      end if
      if ( idiag > 0 ) then
        call assignpnt(qdiag%qcr,mc2mo%dia_qcr)
        call assignpnt(qdiag%qcl,mc2mo%dia_qcl)
        call assignpnt(qdiag%acr,mc2mo%dia_acr)
      end if
    end if

    call assignpnt(fcc,mc2mo%fcc)
    call assignpnt(mc2mo%fcc,mc2mo_fcc)
!$acc enter data create(mc2mo_fcc)
    if ( idynamic == 3 ) then
      call assignpnt(mo_atm%tten,mc2mo%tten)
      call assignpnt(mo_atm%qxten,mc2mo%qxten)
    else
      call assignpnt(aten%t,mc2mo%tten,pc_physic)
      call assignpnt(aten%qx,mc2mo%qxten,pc_physic)
    end if
    call assignpnt(sfs%rainnc,mc2mo%rainnc)
    call assignpnt(sfs%snownc,mc2mo%snownc)
    call assignpnt(pptnc,mc2mo%lsmrnc)
    call assignpnt(rain_ls,mc2mo%rainls)
    call assignpnt(remrat,mc2mo%remrat)
    call assignpnt(rembc,mc2mo%rembc)
    call assignpnt(ncrrate,mc2mo%trrate)

    select case ( ipptls )
      case (1)
        call init_subex(mddom%xlat)
      case (2)
        call init_nogtom(mddom%ldmsk)
      case(3)
        call init_wsm5
      case default
        return
    end select
!$acc update device(mo2mc_qcn, mo2mc_qin, mo2mc_qrn, mo2mc_qsn, mo2mc_phs, mo2mc_qvn, mo2mc_qs, mo2mc_rh, mc2mo_fcc, mo2mc_t, mo2mc_rho, mo2mc_ldmsk, mo2mc_z, mo2mc_ps2, mo2mc_iveg, atms_th700, atms_th3d)
  end subroutine init_micro

  subroutine microscheme
    implicit none
    select case ( ipptls )
      case (1)
        call subex(mo2mc,mc2mo)
      case (2)
#ifdef DEBUG
        call nogtom(mo2mc,mc2mo,ngs)
#else
        call nogtom(mo2mc,mc2mo)
#endif
      case (3)
        call wsm5(mo2mc,mc2mo)
      case default
        return
    end select
  end subroutine microscheme
  !
  ! This subroutine computes the fractional cloud coverage and
  ! liquid water content (in cloud value).  Both are use in
  ! radiation.
  !
  subroutine cldfrac(cldlwc,cldfra)
    use mod_atm_interface , only : atms
    implicit none
    real(rkx) :: exlwc
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: cldlwc , cldfra
    integer(ik4) :: i , j , k , ichi

    if ( ipptls > 1 ) then
      if ( icldfrac == 3 ) then
!$acc parallel present(totc, mo2mc_qcn, mo2mc_qin, mo2mc_qrn, mo2mc_qsn)
!$acc loop collapse(3)        
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              totc(j,i,k) = mo2mc_qcn(j,i,k) + mo2mc_qin(j,i,k) + &
                            mo2mc_qrn(j,i,k) + mo2mc_qsn(j,i,k)
            end do
          end do
        end do
!$acc end parallel
      else
!$acc parallel present(totc, mo2mc_qcn, mo2mc_qin)
!$acc loop collapse(3)
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              totc(j,i,k) = (mo2mc_qcn(j,i,k) + alphaice*mo2mc_qin(j,i,k))
            end do
          end do
        end do
!$acc end parallel
      end if
    else
!$acc parallel present(totc, mo2mc_qcn)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            totc(j,i,k) = mo2mc_qcn(j,i,k)
          end do
        end do
      end do
!$acc end parallel
    end if

    select case ( icldfrac )
      case (1)
        call xuran_cldfrac(mo2mc_phs,totc,mo2mc_qvn, &
                           mo2mc_qs,mo2mc_rh,mc2mo_fcc)
      case (2)
        call thomp_cldfrac(mo2mc_phs,mo2mc_t,mo2mc_rho,mo2mc_qvn, &
                           totc,mo2mc_qsn,mo2mc_qin,mo2mc_ldmsk,  &
                           ds,mc2mo_fcc)
      case (3)
        call gulisa_cldfrac(totc,mo2mc_z,mc2mo_fcc)
      case (4)
        call texeira_cldfrac(totc,mo2mc_qs,mo2mc_rh,rh0,mc2mo_fcc)
      case (5)
        call tompkins_cldfrac(totc,mo2mc_rh,mo2mc_phs,mo2mc_ps2,mc2mo_fcc)
      case (6)
        call echam5_cldfrac(totc,mo2mc_rh,mo2mc_phs,mo2mc_ps2,mc2mo_fcc)
      case default
        call subex_cldfrac(mo2mc_t,mo2mc_phs,mo2mc_qvn, &
                           totc,mo2mc_rh,tc0,rh0,mc2mo_fcc)
    end select

    !------------------------------------------
    ! 1a. Determine Marine stratocumulus clouds
    !------------------------------------------

    if ( icldmstrat == 1 ) then
!$acc parallel present(mo2mc_iveg, mo2mc_phs, mc2mo_fcc, atms_th700, atms_th3d)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( mo2mc_iveg(j,i) == 15 ) then
              if ( mo2mc_phs(j,i,k) >= 70000.0_rkx ) then
                ! Klein, S. A., and D. L. Hartmann,
                ! The seasonal cycle of low stratiform clouds,
                ! J. Climate, 6, 1587-1606, 1993
                mc2mo_fcc(j,i,k) = max(mc2mo_fcc(j,i,k), &
                      min(((atms_th700(j,i)-atms_th3d(j,i,k)) * &
                           0.057_rkx) - 0.5573_rkx,1.0_rkx))
              end if
            end if
          end do
        end do
      end do
!$acc end parallel
    end if

!$acc parallel present(mc2mo_fcc)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          mc2mo_fcc(j,i,k) = max(min(mc2mo_fcc(j,i,k),hicld),d_zero)
        end do
      end do
    end do
!$acc end parallel

    !-----------------------------------------------------------------
    ! 2.  Combine large-scale and convective fraction and liquid water
    !     to be passed into radiation.
    !-----------------------------------------------------------------

    if ( iconvlwp == 1 ) then
!$acc parallel present(mc2mo_fcc, mo2mc_t, cldfra, cldlwc) private(exlwc)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            exlwc = d_zero
            ! Cloud Water Volume
            if ( mc2mo_fcc(j,i,k) > lowcld ) then
              ! Apply the parameterisation based on temperature to the
              ! the large scale clouds.
              exlwc = clwfromt(mo2mc_t(j,i,k))
            end if
            if ( cldfra(j,i,k) > lowcld ) then
              ! get maximum cloud fraction between cumulus and large scale
              cldlwc(j,i,k) = (exlwc * mc2mo_fcc(j,i,k) + &
                              cldlwc(j,i,k) * cldfra(j,i,k)) / &
                              (cldfra(j,i,k) + mc2mo_fcc(j,i,k))
              cldfra(j,i,k) = max(cldfra(j,i,k),mc2mo_fcc(j,i,k))
            else
              cldfra(j,i,k) = mc2mo_fcc(j,i,k)
              cldlwc(j,i,k) = exlwc
            end if
            if ( cldlwc(j,i,k) > d_zero ) then
              cldfra(j,i,k) = min(max(cldfra(j,i,k),d_zero),hicld)
            else
              cldfra(j,i,k) = d_zero
            end if
          end do
        end do
      end do
!$acc end parallel
    else
      if ( any(icup > 1) ) then
!$acc parallel present(mc2mo_fcc, totc, mo2mc_rho, chis, cldfra, cldlwc) private(exlwc, ichi)
!$acc loop collapse(3)
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              exlwc = d_zero
              ! Cloud Water Volume
              ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
              if ( mc2mo_fcc(j,i,k) > lowcld ) then
                exlwc = (totc(j,i,k)*d_1000)*mo2mc_rho(j,i,k)
                ! NOTE : IN CLOUD HERE IS NEEDED !!!
                exlwc = exlwc/mc2mo_fcc(j,i,k)
              end if
              ! Scaling for CF
              ! Implements CF scaling as in Liang GRL 32, 2005
              ! doi: 10.1029/2004GL022301
              ichi = int(mc2mo_fcc(j,i,k)*real(nchi-1,rkx))
              exlwc = exlwc * chis(ichi)
              if ( cldfra(j,i,k) > lowcld ) then
                ! get maximum cloud fraction between cumulus and large scale
                cldlwc(j,i,k) = (exlwc * mc2mo_fcc(j,i,k) + &
                                cldlwc(j,i,k) * cldfra(j,i,k)) / &
                                (cldfra(j,i,k) + mc2mo_fcc(j,i,k))
                cldfra(j,i,k) = max(cldfra(j,i,k),mc2mo_fcc(j,i,k))
              else
                cldfra(j,i,k) = mc2mo_fcc(j,i,k)
                cldlwc(j,i,k) = exlwc
              end if
              if ( cldlwc(j,i,k) > d_zero ) then
                cldfra(j,i,k) = min(max(cldfra(j,i,k),d_zero),hicld)
              else
                cldfra(j,i,k) = d_zero
              end if
            end do
          end do
        end do
!$acc end parallel
      else
!$acc parallel present(mc2mo_fcc, totc, mo2mc_rho, chis, cldlwc, cldfra) private(exlwc, ichi)
!$acc loop collapse(3)
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              exlwc = d_zero
              ! Cloud Water Volume
              ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
              if ( mc2mo_fcc(j,i,k) > lowcld ) then
                exlwc = (totc(j,i,k)*d_1000)*mo2mc_rho(j,i,k)
                ! NOTE : IN CLOUD HERE IS NEEDED !!!
                exlwc = exlwc/mc2mo_fcc(j,i,k)
              end if
              ! Scaling for CF
              ! Implements CF scaling as in Liang GRL 32, 2005
              ! doi: 10.1029/2004GL022301
              ichi = int(mc2mo_fcc(j,i,k)*real(nchi-1,rkx))
              exlwc = exlwc * chis(ichi)
              cldlwc(j,i,k) = exlwc
              if ( cldlwc(j,i,k) > d_zero ) then
                cldfra(j,i,k) = min(max(mc2mo_fcc(j,i,k),d_zero),hicld)
              else
                cldfra(j,i,k) = d_zero
              end if
            end do
          end do
        end do
!$acc end parallel
      end if
    end if

    if ( ipptls == 1 ) then
!$acc kernels present(rh0adj, mo2mc_rh, cldfra)
      rh0adj = d_one - (d_one-mo2mc_rh)/(d_one-cldfra)**2
      rh0adj = max(0.0_rkx,min(rh0adj,0.99999_rkx))
!$acc end kernels
    end if

    contains

#include <clwfromt.inc>

  end subroutine cldfrac
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                 c
  ! This subroutine computes the condensational or evaporational    c
  ! heating term from the explicit moisture scheme.                 c
  !                                                                 c
  ! ---the condensational or evaporational term are one step        c
  !    adjustment based on asai (1965, j. meteo. soc. japan).       c
  !                                                                 c
  ! ---modified to include the effects of partial cloud cover       c
  !    (see Pal et al 2000).  When partial clouds exist, the qxten  c
  !    in/out of the clear and cloudy portions of the grid cell is  c
  !    assumed to be at the same rate (i.e., if there is 80% cloud  c
  !    cover, .2 of qxten goes to raising qv in the clear region    c
  !    and .8 goes to condensation or evaporation of qc in the      c
  !    cloudy portion).                                             c
  !                                                                 c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine condtq
    use mod_atm_interface , only : mo_atm , atm0 , atm2 , sfs , aten
    implicit none
    !
    ! rhc    - Relative humidity at ktau+1
    !
    real(rkx) :: qccs , qvcs , tmp1 , tmp2 , tmp3
    real(rkx) :: dqv , exces , pres , qvc_cld , qvs , fccc , &
               r1 , rhc , rlv , cpm
    integer(ik4) :: i , j , k

    !---------------------------------------------------------------------
    !     1.  Compute t, qv, and qc at tau+1 without condensational term
    !---------------------------------------------------------------------
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( idynamic == 3 ) then
            tmp3 = mo_atm%t(j,i,k) + dt*mo_atm%tten(j,i,k)
            qvcs = max(mo_atm%qx(j,i,k,iqv) + dt*mo_atm%qxten(j,i,k,iqv),minqq)
            qccs = max(mo_atm%qx(j,i,k,iqc) + dt*mo_atm%qxten(j,i,k,iqc),d_zero)
            pres = mo_atm%p(j,i,k)
            qvc_cld = max((mo2mc_qs(j,i,k) + &
                       dt * mc2mo%qxten(j,i,k,iqv)),minqq)
          else
            tmp3 = (atm2%t(j,i,k)+dt*aten%t(j,i,k,pc_total))/sfs%psc(j,i)
            qvcs = atm2%qx(j,i,k,iqv) + dt*aten%qx(j,i,k,iqv,pc_total)
            qccs = atm2%qx(j,i,k,iqc) + dt*aten%qx(j,i,k,iqc,pc_total)
            qvc_cld = max((mo2mc_qs(j,i,k) + &
                       dt * mc2mo%qxten(j,i,k,iqv)/sfs%psc(j,i)),minqq)
            if ( idynamic == 1 ) then
              pres = (hsigma(k)*sfs%psc(j,i)+ptop)*d_1000
            else
              pres = atm0%pr(j,i,k) + &
                 (atm2%pp(j,i,k)+dt*aten%pp(j,i,k,pc_total))/sfs%psc(j,i)
            end if
            if ( qvcs < minqq * sfs%psc(j,i) ) then
              qvcs = minqq * sfs%psc(j,i)
            end if
            if ( qccs < dlowval * sfs%psc(j,i) ) then
              qccs = d_zero
            end if
            qvcs = qvcs /sfs%psc(j,i)
            qccs = qccs /sfs%psc(j,i)
          end if
#ifdef DEBUG
          if ( tmp3 < d_zero ) then
            write(stderr,*) 'Time step = ', rcmtimer%lcount
            write(stderr,*) 'Consistency TEMPERATURE ERROR in condtq (T < 0K)'
            write(stderr,*) 'At global J : ',j
            write(stderr,*) 'At global I : ',i
            write(stderr,*) 'At global K : ',k
          end if
#endif
          !
          ! 2.  Compute the cloud condensation/evaporation term.
          !
          ! 2a. Calculate the saturation mixing ratio and relative humidity
          qvs = pfwsat(tmp3,pres)
          rlv = wlh(tmp3)
          cpm = cpd*(d_one-qvcs) + cpv*qvcs
          r1 = d_one/(d_one+rlv*rlv*qvs/(rwat*cpm*tmp3*tmp3))
          rhc = min(max(qvcs/qvs,d_zero),d_one)
          ! 2b. Compute the relative humidity threshold at ktau+1
          if ( rhc < rh0adj(j,i,k) ) then  ! Low cloud cover
            dqv = conf * (qvcs - qvs)
          else if ( rhc > 0.99999_rkx ) then
            dqv = conf * (qvcs - qvs)      ! High cloud cover
          else
            fccc = d_one-sqrt((d_one-rhc)/(d_one-rh0adj(j,i,k)))
            fccc = min(max(fccc,d_zero),d_one)
            ! qv diff between predicted qv_c
            dqv = conf * fccc * (qvc_cld - qvs)
          end if

          ! 2c. Compute the water vapor in excess of saturation
          tmp1 = r1*dqv               ! grid cell average

          ! 2d. Compute the new cloud water + old cloud water
          exces = qccs + tmp1
          if ( exces >= d_zero ) then ! Some cloud is left
            tmp2 = tmp1/dt
          else                        ! The cloud evaporates
            tmp2 = -qccs/dt
          end if
          !
          ! 3. Compute the tendencies.
          !
          if ( abs(tmp2) > dlowval ) then
            if ( idynamic == 3 ) then
              mo_atm%qxten(j,i,k,iqv) = mo_atm%qxten(j,i,k,iqv) - tmp2
              mo_atm%qxten(j,i,k,iqc) = mo_atm%qxten(j,i,k,iqc) + tmp2
              mo_atm%tten(j,i,k) = mo_atm%tten(j,i,k) + tmp2*rlv/cpm
            else
              aten%qx(j,i,k,iqv,pc_physic) = &
                  aten%qx(j,i,k,iqv,pc_physic) - sfs%psc(j,i)*tmp2
              aten%qx(j,i,k,iqc,pc_physic) = &
                  aten%qx(j,i,k,iqc,pc_physic) + sfs%psc(j,i)*tmp2
              aten%t(j,i,k,pc_physic) = &
                  aten%t(j,i,k,pc_physic) + sfs%psc(j,i)*tmp2*rlv/cpm
            end if
          end if
        end do
      end do
    end do

    contains

#include <pfesat.inc>
#include <pfwsat.inc>
#include <wlh.inc>

  end subroutine condtq

end module mod_micro_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
