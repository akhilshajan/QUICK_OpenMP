#include "config.h"

!  Created by Madu Manathunga on 04/07/2020 
!  Copyright 2020 Michigan State University. All rights reserved.
!
!-------------------------------------------------------
!  scfoperator
!-------------------------------------------------------
!  04/07/2020 Madu Manathunga: Reorganized and improved content 
!                             written by previous authors, implemented
!                             new dft version
!  03/21/2007 Alessandro Genoni: Implemented ECP integral contribution
!                       for operator matrix
!  11/27/2001 Ed Brothers: wrote the original code
!-------------------------------------------------------

subroutine uscf_operator(oneElecO,deltaO)
   use allmod
   use quick_oshell_module
   use quick_cutoff_module
   implicit double precision(a-h,o-z)

   double precision :: oneElecO(nbasis,nbasis)
   logical :: deltaO

! The purpose of this subroutine is to form the operator matrices
! for a uscf calculation, i.e. alpha & beta Fock matrices.  The
! Fock matrix is as follows:

! O(I,J) =  F(I,J) = KE(I,J) + IJ attraction to each atom + repulsion_prim
! with each possible basis  - 1/2 exchange with each
! possible basis.
! Note that the Fock matrix is symmetric.

!-----------------------------------------------------------------
!  Step 1. evaluate 1e integrals
!-----------------------------------------------------------------
   call copyDMat(oneElecO,quick_qm_struct%o,nbasis)
   call copyDMat(oneElecO,quick_qm_struct%ob,nbasis)
 
   if(quick_method%printEnergy) call get1eEnergy_uscf(oneElecO)

! Alessandro GENONI 03/21/2007
! Sum the ECP integrals to the partial Fock matrix
!
!    if (quick_method%ecp) then
!      call ecpoperator
!    end if

!  if only calculate operation difference
   if (deltaO) then
!     save density matrix
      call CopyDMat(quick_qm_struct%dense,quick_qm_struct%denseSave,nbasis)
      call CopyDMat(quick_qm_struct%oSave,quick_qm_struct%o,nbasis)

      call CopyDMat(quick_qm_struct%denseb,quick_qm_struct%densebSave,nbasis)
      call CopyDMat(quick_qm_struct%obSave,quick_qm_struct%ob,nbasis)

      do I=1,nbasis; do J=1,nbasis
         quick_qm_struct%dense(J,I)=quick_qm_struct%dense(J,I)-quick_qm_struct%denseOld(J,I)
         quick_qm_struct%denseb(J,I)=quick_qm_struct%denseb(J,I)-quick_qm_struct%densebOld(J,I)
      enddo; enddo

   endif

!-----------------------------------------------------------------
! Step 2. evaluate 2e integrals
!-----------------------------------------------------------------
!
! The previous two terms are the one electron part of the Fock matrix.
! The next two terms define the two electron part.
!-----------------------------------------------------------------

   call oshell_density_cutoff

   call cpu_time(timer_begin%T2e)

#ifdef CUDA
   if (quick_method%bCUDA) then

      if(quick_method%HF)then
         call gpu_upload_method(0, 1.0d0)
      elseif(quick_method%uselibxc)then
        call gpu_upload_method(3, quick_method%x_hybrid_coeff)
      elseif(quick_method%BLYP)then
         call gpu_upload_method(2, 0.0d0)
      elseif(quick_method%B3LYP)then
         call gpu_upload_method(1, 0.2d0)
      endif

      call gpu_upload_calculated(quick_qm_struct%o,quick_qm_struct%co, &
      quick_qm_struct%vec,quick_qm_struct%dense)
      call gpu_upload_cutoff(cutmatrix,quick_method%integralCutoff,quick_method%primLimit)

      call gpu_upload_calculated_beta(quick_qm_struct%ob, quick_qm_struct%denseb)

      call gpu_get_oshell_eri(quick_qm_struct%o, quick_qm_struct%ob)

   endif   
   
#else

    do II=1,jshell
        call get_oshell_eri(II)
    enddo

#endif

!  Remember the operator is symmetric
   call copySym(quick_qm_struct%o,nbasis)
   call copySym(quick_qm_struct%ob,nbasis)

   if(quick_method%printEnergy) call get_oshell_eri_energy

   call cpu_time(timer_end%T2e)

   timer_cumer%T2e=timer_cumer%T2e+timer_end%T2e-timer_begin%T2e

!-----------------------------------------------------------------
!  Step 3. If DFT, evaluate the exchange/correlation contribution 
!          to the operator
!-----------------------------------------------------------------

   if (quick_method%DFT) then

#ifdef MPIV
   if(master) then
#endif
!  Start the timer for exchange correlation calculation
      call cpu_time(timer_begin%TEx)
#ifdef MPIV
   endif
#endif

!  Calculate exchange correlation contribution & add to operator    
      call getxc_oshell

!  Remember the operator is symmetric
      call copySym(quick_qm_struct%o,nbasis)
      call copySym(quick_qm_struct%ob,nbasis)
#ifdef MPIV
   if(master) then
#endif

!  Stop the exchange correlation timer
      call cpu_time(timer_end%TEx)

!  Add time total time
      timer_cumer%TEx=timer_cumer%TEx+timer_end%TEx-timer_begin%TEx
      write(*,*) "XC time for this step (s):", timer_end%TEx-timer_begin%TEx
   endif

#ifdef MPIV
   endif
#endif

end subroutine uscf_operator



subroutine getxc_oshell
!----------------------------------------------------------------
!  The purpose of this subroutine is to calculate the exchange
!  correlation contribution to the Fock operator. 
!  The angular grid code came from CCL.net.  The radial grid
!  formulas (position and wieghts) is from Gill, Johnson and Pople,
!  Chem. Phys. Lett. v209, n 5+6, 1993, pg 506-512.  The weighting scheme
!  is from Stratmann, Scuseria, and Frisch, Chem. Phys. Lett., v 257,
!  1996, pg 213-223.
!
!  The actual element is:
!  F alpha mu nu = Integral((df/drhoa Phimu Phinu)+
!  (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))
!  where F alpha mu nu is the the alpha spin portion of the operator matrix
!  element mu, nu,
!  df/drhoa is the derivative of the functional by the alpha density,
!  df/dgaa is the derivative of the functional by the alpha gradient
!  invariant, i.e. the dot product of the gradient of the alpha
!  density with itself.
!  df/dgab is the derivative of the functional by the dot product of
!  the gradient of the alpha density with the beta density.
!  Grad(Phimu Phinu) is the gradient of Phimu times Phinu. 
!----------------------------------------------------------------
   use allmod
   use xc_f90_types_m
   use xc_f90_lib_m
   implicit double precision(a-h,o-z)

#ifdef MPIV
   include "mpif.h"
#endif

   double precision, dimension(2) :: libxc_rho
   double precision, dimension(3) :: libxc_sigma
   double precision, dimension(1) :: libxc_exc
   double precision, dimension(2) :: libxc_vrho
   double precision, dimension(3) :: libxc_vsigma

   type(xc_f90_pointer_t), dimension(quick_method%nof_functionals) :: xc_func
   type(xc_f90_pointer_t), dimension(quick_method%nof_functionals) :: xc_info
#ifdef MPIV
   double precision, allocatable:: temp2d(:,:)
   integer :: ierror

   double precision :: Eelxc, Eelxcslave
   allocate(temp2d(nbasis,nbasis))

!  Braodcast libxc information to slaves
   call MPI_BCAST(quick_method%nof_functionals,1,mpi_integer,0,MPI_COMM_WORLD,mpierror)
   call MPI_BCAST(quick_method%functional_id,size(quick_method%functional_id),mpi_integer,0,MPI_COMM_WORLD,mpierror)
   call MPI_BCAST(quick_method%xc_polarization,1,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
#endif

   quick_qm_struct%aelec=0.d0
   quick_qm_struct%belec=0.d0

#ifdef MPIV
!  Set the values of slave operators to zero
   if (.not.master) then
      call zeroMatrix(quick_qm_struct%o, nbasis)
      Eelxc=0
   endif
   call zeroMatrix(temp2d, nbasis)
#endif

#ifdef CUDA

   if(quick_method%bCUDA) then
      call gpu_upload_calculated(quick_qm_struct%o,quick_qm_struct%co, &
            quick_qm_struct%vec,quick_qm_struct%dense)

      call gpu_getxc_new_imp(Eelxc, quick_qm_struct%aelec, quick_qm_struct%belec, quick_qm_struct%o, &
      quick_method%nof_functionals, quick_method%functional_id, quick_method%xc_polarization)

   endif
#else

   if(quick_method%uselibxc) then
!  Initiate the libxc functionals
      do ifunc=1, quick_method%nof_functionals
         call xc_f90_func_init(xc_func(ifunc), xc_info(ifunc), &
              quick_method%functional_id(ifunc),XC_POLARIZED)
      enddo
   endif


#ifdef MPIV
      if(bMPI) then
         irad_init = quick_dft_grid%igridptll(mpirank+1)
         irad_end = quick_dft_grid%igridptul(mpirank+1)
      else
         irad_init = 1
         irad_end = quick_dft_grid%nbins
      endif

   do Ibin=irad_init, irad_end

#else
    do Ibin=1, quick_dft_grid%nbins
#endif

        Igp=quick_dft_grid%bin_counter(Ibin)+1

        do while(Igp < quick_dft_grid%bin_counter(Ibin+1)+1)

           gridx=quick_dft_grid%gridxb(Igp)
           gridy=quick_dft_grid%gridyb(Igp)
           gridz=quick_dft_grid%gridzb(Igp)

           sswt=quick_dft_grid%gridb_sswt(Igp)
           weight=quick_dft_grid%gridb_weight(Igp)
           Iatm=quick_dft_grid%gridb_atm(Igp)

            if (weight < quick_method%DMCutoff ) then
               continue
            else

               icount=quick_dft_grid%basf_counter(Ibin)+1
               do while (icount < quick_dft_grid%basf_counter(Ibin+1)+1)
               Ibas=quick_dft_grid%basf(icount)+1
                  call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                  dphidz,Ibas)
                  phixiao(Ibas)=phi
                  dphidxxiao(Ibas)=dphidx
                  dphidyxiao(Ibas)=dphidy
                  dphidzxiao(Ibas)=dphidz

                  icount=icount+1
               enddo

!  Next, evaluate the densities at the grid point and the gradient
!  at that grid point.

               call denspt_oshell(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
               gbx,gby,gbz,Ibin)

!                write(*,*) gridx,gridy,gridz,density,densityb,gax,gbx,gay,gby,gaz,gbz

!               if ((density+densityb) < quick_method%DMCutoff) then
               if ((density < quick_method%DMCutoff) .and. (densityb < quick_method%DMCutoff)) then
                  continue
               else

!  This allows the calculation of the derivative of the functional with regard to the 
!  density (dfdr), with regard to the alpha-alpha density invariant (df/dgaa), and the
!  alpha-beta density invariant.

                  gaa = (gax*gax+gay*gay+gaz*gaz)
                  gbb = (gbx*gbx+gby*gby+gbz*gbz)
                  gab = (gax*gbx+gay*gby+gaz*gbz)

                  libxc_rho(1)=density
                  libxc_rho(2)=densityb

                  libxc_sigma(1)=gaa
                  libxc_sigma(2)=gab
                  libxc_sigma(3)=gbb

                  excpp=0.0d0
                  dfdr=0.0d0
                  dfdrb=0.0d0

                  dfdgaa=0.0d0
                  dfdgab=0.0d0
                  dfdgbb=0.0d0

                  if(quick_method%uselibxc) then
                     do ifunc=1, quick_method%nof_functionals
                        select case(xc_f90_info_family(xc_info(ifunc)))
                           case(XC_FAMILY_LDA)
!                              call xc_f90_lda_exc_vxc(xc_func(ifunc),1,libxc_rho(1), &
!                              libxc_exc(1), libxc_vrhoa(1))

                              libxc_vsigma(1) = 0.0d0
                              libxc_vsigma(2) = 0.0d0
                              libxc_vsigma(3) = 0.0d0

                           case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                              call xc_f90_gga_exc_vxc(xc_func(ifunc),1,libxc_rho(1),libxc_sigma(1), &
                              libxc_exc(1),libxc_vrho(1),libxc_vsigma(1))
                        end select

                        excpp=excpp+libxc_exc(1)
                        dfdr=dfdr+libxc_vrho(1)
                        dfdrb=dfdrb+libxc_vrho(2)

                        dfdgaa=dfdgaa+libxc_vsigma(1)
                        dfdgab=dfdgab+libxc_vsigma(2)
                        dfdgbb=dfdgbb+libxc_vsigma(3)

                     enddo

                     zkec=(density+densityb)*excpp

!  Calculate the first term in the dot product shown above,
!  i.e.: (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                     xdot = 2.d0*dfdgaa*gax + dfdgab*gbx
                     ydot = 2.d0*dfdgaa*gay + dfdgab*gby
                     zdot = 2.d0*dfdgaa*gaz + dfdgab*gbz

                     xdotb = 2.d0*dfdgbb*gbx + dfdgab*gax
                     ydotb = 2.d0*dfdgbb*gby + dfdgab*gay
                     zdotb = 2.d0*dfdgbb*gbz + dfdgab*gaz

                  elseif(quick_method%BLYP) then

!                     call becke_E(density, densityb, gax, gay, gaz, gbx, gby,gbz, Ex)
!                     call lyp_e(density, densityb, gax, gay, gaz, gbx, gby, gbz,Ec)

!                     zkec=Ex+Ec

!write(*,*) density, densityb, gax, gay, gaz, gbx, gby,gbz, Ex, Ec

!                     call becke(density, gax, gay, gaz, gbx, gby, gbz, dfdr, dfdgaa, dfdgab)
!                     call lyp(density, densityb, gax, gay, gaz, gbx, gby, gbz, dfdr2, dfdgaa2, dfdgab2)

!                     call becke(densityb, gax, gay, gaz, gbx, gby, gbz, dfdrb,dfdgbb, dfdgabb)
!                     call lyp(density, densityb, gax, gay, gaz, gbx, gby, gbz,dfdrb2, dfdgbb2, dfdgabb2)

!                     dfdr = dfdr + dfdr2
!                     dfdgaa = dfdgaa + dfdgaa2
!                     dfdgab = dfdgab + dfdgab2

!                     dfdrb = dfdrb + dfdrb2
!                     dfdgbb = dfdgbb + dfdgbb2
!                     dfdgabb = dfdgabb + dfdgabb2

!                     xdot = 2.d0*dfdgaa*gax + dfdgab*gbx
!                     ydot = 2.d0*dfdgaa*gay + dfdgab*gby
!                     zdot = 2.d0*dfdgaa*gaz + dfdgab*gbz

!                     xdotb = 2.d0*dfdgbb*gbx + dfdgabb*gax
!                     ydotb = 2.d0*dfdgbb*gby + dfdgabb*gay
!                     zdotb = 2.d0*dfdgbb*gbz + dfdgabb*gaz

                  elseif(quick_method%B3LYP) then

!                     call b3lyp_e(density, sigma, zkec)
!                     call b3lypf(density, sigma, dfdr, xiaodot)
!
!                     call b3lyp_e(densityb, sigmab, zkecb)
!                     call b3lypf(densityb, sigmab, dfdrb, xiaodotb)
!
!write(*,*) density,densityb,sigma,sigmab,zkec,zkecb,dfdr,dfdrb,xiaodot,xiaodotb
!
!                     zkec=zkec+zkecb
!
!                     xdot=xiaodot*gax
!                     ydot=xiaodot*gay
!                     zdot=xiaodot*gaz
!
!                     xdotb=xiaodotb*gbx
!                     ydotb=xiaodotb*gby
!                     zdotb=xiaodotb*gbz
!
                  endif

                  Eelxc = Eelxc + zkec*weight

                  quick_qm_struct%aelec = weight*density+quick_qm_struct%aelec
                  quick_qm_struct%belec = weight*densityb+quick_qm_struct%belec

!  Now loop over basis functions and compute the addition to the matrix element.
!                  do Ibas=1,nbasis
                  icount=quick_dft_grid%basf_counter(Ibin)+1
                  do while (icount < quick_dft_grid%basf_counter(Ibin+1)+1)
                  Ibas=quick_dft_grid%basf(icount)+1

                     phi=phixiao(Ibas)
                     dphidx=dphidxxiao(Ibas)
                     dphidy=dphidyxiao(Ibas)
                     dphidz=dphidzxiao(Ibas)
                     quicktest = DABS(dphidx+dphidy+dphidz+phi)

                     if (quicktest < quick_method%DMCutoff ) then
                        continue
                     else
                        jcount=icount
                        do while(jcount<quick_dft_grid%basf_counter(Ibin+1)+1)
                        Jbas = quick_dft_grid%basf(jcount)+1
                           phi2=phixiao(Jbas)
                           dphi2dx=dphidxxiao(Jbas)
                           dphi2dy=dphidyxiao(Jbas)
                           dphi2dz=dphidzxiao(Jbas)
                           temp = phi*phi2
                           tempgx = phi*dphi2dx + phi2*dphidx
                           tempgy = phi*dphi2dy + phi2*dphidy
                           tempgz = phi*dphi2dz + phi2*dphidz

                           quick_qm_struct%o(Jbas,Ibas)=quick_qm_struct%o(Jbas,Ibas)+(temp*dfdr+&
                           xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight

                           quick_qm_struct%ob(Jbas,Ibas)=quick_qm_struct%ob(Jbas,Ibas)+(temp*dfdrb+&
                           xdotb*tempgx+ydotb*tempgy+zdotb*tempgz)*weight                           

                           jcount=jcount+1
                        enddo
                     endif
                     icount=icount+1
                  enddo
               endif
            endif

         Igp=Igp+1
      enddo
   enddo

   if(quick_method%uselibxc) then
!  Uninitilize libxc functionals
      do ifunc=1, quick_method%nof_functionals
         call xc_f90_func_end(xc_func(ifunc))
      enddo
   endif
#endif

#ifdef MPIV
   if(.not. master) then
!  Send the Exc energy value
      Eelxcslave=Eelxc
      call MPI_SEND(Eelxcslave,1,mpi_double_precision,0,mpirank,MPI_COMM_WORLD,IERROR)
      call copyDMat(quick_qm_struct%o,temp2d,nbasis)
      call MPI_SEND(temp2d,nbasis*nbasis,mpi_double_precision,0,mpirank,MPI_COMM_WORLD,IERROR)
   else

!  Master node will receive infos from every nodes
      do i=1,mpisize-1
!  Receive exchange correlation energy from slaves
         call MPI_RECV(Eelxcslave,1,mpi_double_precision,i,i,MPI_COMM_WORLD,MPI_STATUS,IERROR)
         Eelxc=Eelxc+Eelxcslave
!  Receive opertors from slave nodes
         call  MPI_RECV(temp2d,nbasis*nbasis,mpi_double_precision,i,i,MPI_COMM_WORLD,MPI_STATUS,IERROR)
!  Sum them into operator
         do ii=1,nbasis
            do jj=1,nbasis
               quick_qm_struct%o(ii,jj)=quick_qm_struct%o(ii,jj)+temp2d(ii,jj)
            enddo
         enddo
      enddo
   endif
   call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
#endif

#ifdef MPIV
   if(master) then
#endif
!  Add the exchange correlation energy to total electronic energy
   quick_qm_struct%Eel=quick_qm_struct%Eel+Eelxc

!   if(quick_method%debug) then
      write(*,*) "Eelex=",Eelxc
      write(*,*) "E1+E2+Eelxc=",quick_qm_struct%Eel
!   endif
#ifdef MPIV
   endif
#endif

! call exit

   return

end subroutine getxc_oshell

