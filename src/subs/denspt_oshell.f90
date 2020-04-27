










!---------------------------------------------------------------------!
! Created by Madu Manathunga on 04/19/2020                            !
!                                                                     !
! Previous contributors: Yipu Miao, Xio He, Alessandro Genoni,        !
!                         Ken Ayers & Ed Brothers                     !
!                                                                     ! 
! Copyright (C) 2020-2021 Merz lab                                    !
! Copyright (C) 2020-2021 Götz lab                                    !
!                                                                     !
! This Source Code Form is subject to the terms of the Mozilla Public !
! License, v. 2.0. If a copy of the MPL was not distributed with this !
! file, You can obtain one at http://mozilla.org/MPL/2.0/.            !
!_____________________________________________________________________!

   subroutine denspt_oshell &
   (gridx,gridy,gridz,densitya,densityb,gax,gay,gaz,gbx,gby,gbz,Ibin)
   use allmod
   implicit none
   ! Given a point in space, this function calculates the densities and
   ! gradient  at that point.  The gradients are stored in the common block
   ! three element arrays ga and gb for alpha and beta electron gradients. Thus
   ! the x component of the alpha density is stored in ga(1).
   
   
   ! INPUT PARAMETERS
   double precision :: gridx,gridy,gridz
   double precision :: densitya,densityb
   double precision :: gax,gay,gaz
   double precision :: gbx,gby,gbz
   integer :: Ibin

   ! INNER VARIBLES
   double precision :: dphidx,dphidy,dphidz
   double precision :: phi,phi2
   double precision :: densebij,denseij
   integer :: Ibas,Jbas, icount, jcount

   densitya=0.0d0
   densityb=0.0d0

   gax=0.0d0
   gay=0.0d0
   gaz=0.0d0
   gbx=0.0d0
   gby=0.0d0
   gbz=0.0d0

   icount=quick_dft_grid%basf_counter(Ibin)+1
   do while (icount < quick_dft_grid%basf_counter(Ibin+1)+1)
      Ibas=quick_dft_grid%basf(icount)+1

      DENSEIJ=quick_qm_struct%dense(Ibas,Ibas)
      DENSEBIJ=quick_qm_struct%denseb(Ibas,Ibas)

      if((DABS(quick_qm_struct%dense(Ibas,Ibas)) < quick_method%DMCutoff) .and. &
        (DABS(quick_qm_struct%denseb(Ibas,Ibas)) < quick_method%DMCutoff)) then
!      if((DABS(quick_qm_struct%dense(Ibas,Ibas))+DABS(quick_qm_struct%denseb(Ibas,Ibas))) < quick_method%DMCutoff) then
         continue
      else

         phi=phixiao(Ibas)
         dphidx=dphidxxiao(Ibas)
         dphidy=dphidyxiao(Ibas)
         dphidz=dphidzxiao(Ibas)

         if (DABS(dphidx+dphidy+dphidz+phi) < quick_method%DMCutoff ) then
            continue
         else
            densitya=densitya+DENSEIJ*phi*phi
            densityb=densityb+DENSEBIJ*phi*phi
            
              ! write(*,*) "a",densitya,DENSEIJ,phi,phi,quick_qm_struct%dense(Ibas,Ibas)

            gax=gax+DENSEIJ*phi*dphidx
            gay=gay+DENSEIJ*phi*dphidy
            gaz=gaz+DENSEIJ*phi*dphidz

            gbx=gbx+DENSEBIJ*phi*dphidx
            gby=gby+DENSEBIJ*phi*dphidy
            gbz=gbz+DENSEBIJ*phi*dphidz

            jcount = icount+1
            do while( jcount<quick_dft_grid%basf_counter(Ibin+1)+1)
               Jbas = quick_dft_grid%basf(jcount)+1

               DENSEIJ=quick_qm_struct%dense(Jbas,Ibas)
               DENSEBIJ=quick_qm_struct%denseb(Jbas,Ibas)
               phi2=phixiao(Jbas)
               
               densitya=densitya+2.0d0*DENSEIJ*phi*phi2
               densityb=densityb+2.0d0*DENSEBIJ*phi*phi2

               gax=gax+DENSEIJ*(phi*dphidxxiao(Jbas)+phi2*dphidx)
               gay=gay+DENSEIJ*(phi*dphidyxiao(Jbas)+phi2*dphidy)
               gaz=gaz+DENSEIJ*(phi*dphidzxiao(Jbas)+phi2*dphidz)

               gbx=gbx+DENSEBIJ*(phi*dphidxxiao(Jbas)+phi2*dphidx)
               gby=gby+DENSEBIJ*(phi*dphidyxiao(Jbas)+phi2*dphidy)
               gbz=gbz+DENSEBIJ*(phi*dphidzxiao(Jbas)+phi2*dphidz)

               jcount=jcount+1
            enddo
         endif
      endif

      icount=icount+1 
   enddo

   gax=2.0d0*gax
   gay=2.0d0*gay
   gaz=2.0d0*gaz

   gbx=2.0d0*gbx
   gby=2.0d0*gby
   gbz=2.0d0*gbz

   end subroutine denspt_oshell
