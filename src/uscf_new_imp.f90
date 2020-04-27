#include "config.h"

! Madu Manathunga 03/31/2020

subroutine uscf_new_imp2(failed)
   !-------------------------------------------------------
   ! this subroutine is to do scf job for restricted system
   !-------------------------------------------------------
   use allmod
   implicit double precision(a-h,o-z)

   logical :: done,failed
   done=.false.

   !-----------------------------------------------------------------
   ! The purpose of this subroutine is to perform scf cycles.  At this
   ! point, X has been formed. The remaining steps are:
   ! 1)  Form operator matrix.
   ! 2)  Calculate O' = Transpose[X] O X
   ! 3)  Diagonalize O' to obtain C' and eigenvalues.
   ! 4)  Calculate C = XC'
   ! 5)  Form new density matrix.
   ! 6)  Check for convergence.
   !-----------------------------------------------------------------

   ! Each location in the code that the step is occurring will be marked.
   ! The cycles stop when prms  is less than pmaxrms or when the maximum
   ! number of scfcycles has been reached.
   jscf=0

   if(master) then
      write(ioutfile,'(40x," USCF ENERGY")')
      if (quick_method%printEnergy) then
         write(ioutfile,'(120("-"))')
      else
         write(ioutfile,'(90("-"))')
      endif
      write(ioutfile,'("NCYC",6x)',advance="no")
      if (quick_method%printEnergy) write(ioutfile,'(" ENERGY",8x,"DELTA_E",5x)',advance="no")
      write(ioutfile,'(" SCF_TIME",2x,"DII_CYC",2x," DII_TIME ",2x,"O_TIME",2x,&
            "DIAG_TIME",4x,"MAX_ERR",4x,"RMS_CHG",4x,"MAX_CHG")')
      if (quick_method%printEnergy) then
         write(ioutfile,'(120("-"))')
      else
         write(ioutfile,'(90("-"))')
      endif
   endif

   ! Alessandro GENONI 03/21/2007
   ! ECP integrals computation exploiting Alexander V. Mitin Subroutine
   ! Note: the integrals are stored in the array ecp_int that corresponds
   !       to the lower triangular matrix of the ECP integrals
   !if (quick_method%ecp) call ecpint

   ! if not direct SCF, generate 2e int file
   !if (quick_method%nodirect) call aoint

   if (quick_method%diisscf .and. .not. quick_method%divcon) call uelectdiis_new_imp2(jscf)       ! normal scf

   jscf=jscf+1

   failed = failed.and.(jscf.gt.quick_method%iscf)
   if (quick_method%debug)  call debug_SCF(jscf)

   return

end subroutine uscf_new_imp2


subroutine uelectdiis_new_imp2(jscf)
   use allmod
   implicit none

   integer :: jscf
   double precision :: oneElecO(nbasis,nbasis)

   logical :: diisdone = .false. ! flag to indicate if diis is done
   logical :: deltaO   = .false. ! delta Operator
   integer :: idiis = 0          ! diis iteration
   integer :: lsolerr = 0
   integer :: IDIIS_Error_Start, IDIIS_Error_End,IDIISfinal,iidiis,current_diis
   
   double precision :: Sum2Mat,BIJ,CIJ,DENSEIJ,rms,errormax,HOLDIJ,OIJ,OJK,OLDPRMS,PCHANGE,PRMS,PRMS2,tmp,temp
   double precision :: t1,t2, tempij,DENSEJI
   integer :: i,j,IERROR,k

   double precision :: alloperatorB(quick_method%maxdiisscf,nbasis,nbasis)
   double precision :: B(quick_method%maxdiisscf+1,quick_method%maxdiisscf+1)
   double precision :: BSAVE(quick_method%maxdiisscf+1,quick_method%maxdiisscf+1)
   double precision :: BCOPY(quick_method%maxdiisscf+1,quick_method%maxdiisscf+1)
   double precision :: W(quick_method%maxdiisscf+1), V2(3,nbasis)
   double precision :: COEFF(quick_method%maxdiisscf+1),RHS(quick_method%maxdiisscf+1)

   double precision :: allerror(quick_method%maxdiisscf,nbasis,nbasis)
   double precision :: alloperator(quick_method%maxdiisscf,nbasis,nbasis)

   double precision :: TEne_begin,TEne_end, step_TDiag ! Operator and diagonalization times for a single step
   double precision :: oldEnergy            ! To keep track of step energy change

   ! The purpose of this subroutine is to utilize Pulay's accelerated
   ! scf convergence as detailed in J. Comp. Chem, Vol 3, #4, pg 566-60, 1982.
   ! At the beginning of this process, their is an approximate density
   ! matrix.
   ! This is the unrestricted (Pople-Neesbitt) version of the code.
   ! The step in the procedure are:
   ! 1)  Form the alpha operator matrix for step i, O(i).  (Store in
   ! alloperator array.)
   ! 2)  Form alpha error matrix for step i.
   ! e(i) = Oa Da S - S Da Oa
   ! 3)  Form the beta operator matrix for step i, O(i).  (Store in
   ! alloperatorb array.)
   ! 4)  Form beta error matrix for step i.
   ! e(i) = e(i,alpha part)+Ob Db S - S Db Ob
   ! 5)  Move e to an orthogonal basis.  e'(i) = Transpose[X] .e(i). X
   ! 6)  Store the e'(I) in allerror.
   ! 7)  Form matrix B, which is:

   ! _                                                 _
   ! |                                                   |
   ! |  B(1,1)      B(1,2)     . . .     B(1,J)      -1  |
   ! |  B(2,1)      B(2,2)     . . .     B(2,J)      -1  |
   ! |  .            .                     .          .  |
   ! B = |  .            .                     .          .  |
   ! |  .            .                     .          .  |
   ! |  B(I,1)      B(I,2)     . . .     B(I,J)      -1  |
   ! | -1            -1        . . .      -1          0  |
   ! |_                                                 _|


   ! Where B(i,j) = Trace(e(i) Transpose(e(j)) )

   ! 8)  Solve B*COEFF = RHS which is:

   ! _                                             _  _  _     _  _
   ! |                                               ||    |   |    |
   ! |  B(1,1)      B(1,2)     . . .     B(1,J)  -1  ||  C1|   |  0 |
   ! |  B(2,1)      B(2,2)     . . .     B(2,J)  -1  ||  C2|   |  0 |
   ! |  .            .                     .      .  ||  . |   |  0 |
   ! |  .            .                     .      .  ||  . | = |  0 |
   ! |  .            .                     .      .  ||  . |   |  0 |
   ! |  B(I,1)      B(I,2)     . . .     B(I,J)  -1  ||  Ci|   |  0 |
   ! | -1            -1        . . .      -1      0  || -L |   | -1 |
   ! |_                                             _||_  _|   |_  _|


   ! 9) Form a new alpha operator matrix based on
   ! O(new) = [Sum over i] c(i)O(i)

   ! 10) Diagonalize the operator matrix to form a new density matrix.

   ! 11) Form a new beta operator matrix based on
   ! O(new) = [Sum over i] c(i)O(i)

   ! 12) Diagonalize the operator matrix to form a new density matrix.

   ! As in scf.F, each step wil be reviewed as we pass through the code.

   ! 1)  Form the alpha operator matrix for step i, O(i).  (Store in
   ! alloperator array.)

   call get1e(oneElecO)

#ifdef CUDA
   if(quick_method%bCUDA) then
      if (quick_method%DFT) then
      call gpu_upload_dft_grid(quick_dft_grid%gridxb, quick_dft_grid%gridyb,quick_dft_grid%gridzb, quick_dft_grid%gridb_sswt, &
      quick_dft_grid%gridb_weight, quick_dft_grid%gridb_atm,quick_dft_grid%dweight, quick_dft_grid%basf, quick_dft_grid%primf, &
      quick_dft_grid%basf_counter, quick_dft_grid%primf_counter,quick_dft_grid%gridb_count, quick_dft_grid%nbins,&
      quick_dft_grid%nbtotbf, quick_dft_grid%nbtotpf, quick_method%isg, sigrad2)
      endif
   endif
#endif

   do i=1,nbasis
       do j=1,nbasis
           temp=quick_qm_struct%dense(j,i)/2.0d0
           quick_qm_struct%dense(j,i)=temp
           quick_qm_struct%denseb(j,i)=temp
       enddo
   enddo

   do while (.not.diisdone)

      call cpu_time(timer_begin%TSCF)

      ! Save current total energy 
      oldEnergy=quick_qm_struct%Eel+quick_qm_struct%Ecore


      temp=Sum2Mat(quick_qm_struct%dense,quick_qm_struct%s,nbasis)


      !--------------------------------------------
      ! 1)  Form the operator matrix for step i, O(i).
      !--------------------------------------------

      idiis=idiis+1
      jscf=jscf+1
      
      if(idiis.le.quick_method%maxdiisscf)then
         IDIISfinal=idiis; iidiis=idiis
      else
         IDIISfinal=quick_method%maxdiisscf; iidiis=1
      endif

      !-----------------------------------------------
      ! Before Delta Densitry Matrix, normal operator is implemented here
      !-----------------------------------------------

      ! Triger Operator timer
      call cpu_time(timer_begin%TOp)
     
      ! if want to calculate operator difference?
      if(jscf.ge.quick_method%ncyc) deltaO = .true. 

      call uscf_operator(oneElecO, deltaO)

      ! Terminate Operator timer
      call cpu_time(timer_end%TOp)

      call cpu_time(timer_begin%TDII)

      call CopyDMat(quick_qm_struct%o,quick_qm_struct%oSave,nbasis)
      call CopyDMat(quick_qm_struct%dense,quick_qm_struct%denseOld,nbasis)

      !-----------------------------------------------
      ! 2)  Form error matrix for step i.
      ! e(i) = ODS - SDO
      !-----------------------------------------------
      ! The matrix multiplier comes from Steve Dixon. It calculates
      ! C = Transpose(A) B.  Thus to utilize this we have to make sure that the
      ! A matrix is symetric. First, calculate DENSE*S and store in the scratch
      ! matrix hold. Then calculate O*(DENSE*S).  As the operator matrix is symmetric, the
      ! above code can be used. Store this (the ODS term) in the all error matrix.

#ifdef CUDA
      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%dense, &
            nbasis, quick_qm_struct%s, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%o, &
            nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#else

      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%dense, &
               nbasis, quick_qm_struct%s, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%o, &
               nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#endif

      do i = 1, nbasis
         do j = 1, nbasis
            allerror(iidiis, i, j) = quick_scratch%hold2( i, j)
         enddo
      enddo

      ! Calculate D O. then calculate S (do) and subtract that from the
      ! allerror matrix. This means we now have the e(i) matrix.
      ! allerror=ODS-SDO

#ifdef CUDA
      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%dense, &
            nbasis, quick_qm_struct%o, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%s, &
            nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#else
      quick_scratch%hold=MATMUL(quick_qm_struct%dense,quick_qm_struct%o)
      quick_scratch%hold2=MATMUL(quick_qm_struct%s,quick_scratch%hold)
#endif

         !errormax = 0.d0
      do I=1,nbasis
         do J=1,nbasis
            allerror(iidiis,I,J) = allerror(iidiis,I,J) - quick_scratch%hold2(i,j) !e=ODS=SDO
            !errormax = max(allerror(iidiis,I,J),errormax)
         enddo
      enddo

      ! 3)  Form the beta operator matrix for step i, O(i).  (Store in
      ! alloperatorb array.)

      call CopyDMat(quick_qm_struct%ob,quick_qm_struct%obSave,nbasis)
      call CopyDMat(quick_qm_struct%denseb,quick_qm_struct%densebOld,nbasis)

      ! 4)  Form beta error matrix for step i.
      ! e(i) = e(i,alpha part)+Ob Db S - S Db Ob

      ! First, calculate quick_qm_struct%denseb*S and store in the scratch
      ! matrix hold. Then calculate O*(quick_qm_struct%denseb*S).  As the operator matrix is
      ! symmetric, the above code can be used. Add this (the ODS term) into the allerror
      ! matrix.

#ifdef CUDA
      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%denseb, &
            nbasis, quick_qm_struct%s, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%ob, &
            nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#else
      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%denseb, &
               nbasis, quick_qm_struct%s, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%ob, &
               nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#endif

      do i = 1, nbasis
         do j = 1, nbasis
             allerror(iidiis, i, j) = allerror(iidiis, i, j)+quick_scratch%hold2( i, j)   
         enddo
      enddo

      ! Calculate Db O.Then calculate S (DbO) and subtract that from the allerror matrix.
      ! This means we now have the complete e(i) matrix.
#ifdef CUDA

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%denseb, &
            nbasis, quick_qm_struct%ob, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%s, &
            nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#else
      quick_scratch%hold=MATMUL(quick_qm_struct%denseb,quick_qm_struct%ob)
      quick_scratch%hold2=MATMUL(quick_qm_struct%s,quick_scratch%hold)
#endif

      errormax = 0.d0
      do I=1,nbasis
         do J=1,nbasis
            allerror(iidiis,I,J) = allerror(iidiis,I,J) - quick_scratch%hold2(i,j) !e=ODS=SDO
            errormax = max(allerror(iidiis,I,J),errormax)
         enddo
      enddo

      ! 5)  Move e to an orthogonal basis.  e'(i) = Transpose[X] .e(i). X
      ! X is symmetric, but we do not know anything about the symmetry of e.
      ! The easiest way to do this is to calculate e(i) . X , store
      ! this in HOLD, and then calculate Transpose[X] (.e(i) . X)

      do i = 1, nbasis 
         do j = 1, nbasis
            quick_scratch%hold2( i, j) = allerror(iidiis, i, j)
         enddo
      enddo

#ifdef CUDA

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_scratch%hold2, &
            nbasis, quick_qm_struct%x, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x, &
            nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#else           
      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_scratch%hold2, &
           nbasis, quick_qm_struct%x, nbasis, 0.0d0, quick_scratch%hold,nbasis)
            
      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x,&
           nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_scratch%hold2,nbasis)
#endif

      do i = 1, nbasis
         do j = 1, nbasis
            allerror(iidiis, i, j) = quick_scratch%hold2( i, j)
         enddo
      enddo

      !-----------------------------------------------
      ! 4)  Store the e'(I) and O(i).
      ! e'(i) is already stored.  Simply store the operator matrix in
      ! all operator.
      !-----------------------------------------------

      if(idiis.le.quick_method%maxdiisscf)then
         call CopyDMat(quick_qm_struct%o,alloperator(iidiis,1:nbasis,1:nbasis),nbasis)
         call CopyDMat(quick_qm_struct%ob,alloperatorB(iidiis,1:nbasis,1:nbasis),nbasis)
      else
         do K=1,quick_method%maxdiisscf-1
            call CopyDMat(alloperator(K+1,1:nbasis,1:nbasis),alloperator(K,1:nbasis,1:nbasis),nbasis)
            call CopyDMat(alloperatorB(K+1,1:nbasis,1:nbasis),alloperatorB(K,1:nbasis,1:nbasis),nbasis)
         enddo
         call CopyDMat(quick_qm_struct%o,alloperator(quick_method%maxdiisscf,1:nbasis,1:nbasis),nbasis)
         call CopyDMat(quick_qm_struct%ob,alloperatorB(quick_method%maxdiisscf,1:nbasis,1:nbasis),nbasis)
      endif

      ! 7)  Form matrix B, which is:

      ! _                                                 _
      ! |                                                   |
      ! |  B(1,1)      B(1,2)     . . .     B(1,J)      -1  |
      ! |  B(2,1)      B(2,2)     . . .     B(2,J)      -1  |
      ! |  .            .                     .          .  |
      ! B = |  .            .                     .          .  |
      ! |  .            .                     .          .  |
      ! |  B(I,1)      B(I,2)     . . .     B(I,J)      -1  |
      ! | -1            -1        . . .      -1          0  |
      ! |_                                                 _|


      ! Where B(i,j) = Trace(e(i) Transpose(e(j)))
      ! According to an example done in mathematica, B12 = B21.  Note that
      ! the rigorous proof of this phenomenon is left as an exercise for the
      ! reader.  Thus the first step is copying BCOPY to B.  In this way we
      ! only have to recalculate the new elements.

      do I=1,IDIISfinal
         do J=1,IDIISfinal
            B(J,I) = BCOPY(J,I)
         enddo
      enddo

      if(IDIIS.gt.quick_method%maxdiisscf)then
         do I=1,IDIISfinal-1
            do J=1,IDIISfinal-1
               B(J,I) = BCOPY(J+1,I+1)
            enddo
         enddo
      endif


      ! Now copy the current matrix into HOLD2 transposed.  This will be the
      ! Transpose[ej] used in B(i,j) = Trace(e(i) Transpose(e(j)))

      call CopyDMat(allerror(iidiis,1:nbasis,1:nbasis),quick_scratch%hold2,nbasis)

         do I=1,IDIISfinal
            ! Copy the transpose of error matrix I into HOLD.
            call CopyDMat(allerror(I,1:nbasis,1:nbasis),quick_scratch%hold,nbasis)

            ! Calculate and sum together the diagonal elements of e(i) Transpose(e(j))).
            BIJ=Sum2Mat(quick_scratch%hold2,quick_scratch%hold,nbasis)

            ! Now place this in the B matrix.
            if(idiis.le.quick_method%maxdiisscf)then
               B(iidiis,I) = BIJ
               B(I,iidiis) = BIJ
            else
               if(I.gt.1)then
                  B(quick_method%maxdiisscf,I-1)=BIJ
                  B(I-1,quick_method%maxdiisscf)=BIJ
               else
                  B(quick_method%maxdiisscf,quick_method%maxdiisscf)=BIJ
               endif
            endif
         enddo


         if(idiis.gt.quick_method%maxdiisscf)then
            call CopyDMat(allerror(1,1:nbasis,1:nbasis),quick_scratch%hold,nbasis)
            do J=1,quick_method%maxdiisscf-1
               call CopyDMat(allerror(J+1,1:nbasis,1:nbasis),allerror(J,1:nbasis,1:nbasis),nbasis)
            enddo
            call CopyDMat(quick_scratch%hold,allerror(quick_method%maxdiisscf,1:nbasis,1:nbasis),nbasis)
         endif

      ! Now that all the BIJ elements are in place, fill in all the column
      ! and row ending -1, and fill up the rhs matrix.

         do I=1,IDIISfinal
            B(I,IDIISfinal+1) = -1.d0
            B(IDIISfinal+1,I) = -1.d0
         enddo
         do I=1,IDIISfinal
            RHS(I) = 0.d0
         enddo
         RHS(IDIISfinal+1) = -1.d0
         B(IDIISfinal+1,IDIISfinal+1) = 0.d0

         ! Now save the B matrix in Bcopy so it is available for subsequent
         ! iterations.
         do I=1,IDIISfinal
            do J=1,IDIISfinal
               BCOPY(J,I)=B(J,I)
            enddo
         enddo

      ! 8)  Solve B*COEFF = RHS which is:

      ! _                                             _  _  _     _  _
      ! |                                               ||    |   |    |
      ! |  B(1,1)      B(1,2)     . . .     B(1,J)  -1  ||  C1|   |  0 |
      ! |  B(2,1)      B(2,2)     . . .     B(2,J)  -1  ||  C2|   |  0 |
      ! |  .            .                     .      .  ||  . |   |  0 |
      ! |  .            .                     .      .  ||  . | = |  0 |
      ! |  .            .                     .      .  ||  . |   |  0 |
      ! |  B(I,1)      B(I,2)     . . .     B(I,J)  -1  ||  Ci|   |  0 |
      ! | -1            -1        . . .      -1      0  || -L |   | -1 |
      ! |_                                             _||_  _|   |_  _|

         call CopyDMat(B,BSAVE,IDIISfinal+1)
         call LSOLVE(IDIISfinal+1,quick_method%maxdiisscf+1,B,RHS,W,quick_method%DMCutoff,COEFF,LSOLERR)

         IDIIS_Error_Start = 1
         IDIIS_Error_End   = IDIISfinal
         111     IF (LSOLERR.ne.0)then
            IDIISfinal=Idiisfinal-1
            do I=1,IDIISfinal+1
               do J=1,IDIISfinal+1
                  B(I,J)=BSAVE(I+IDIIS_Error_Start,J+IDIIS_Error_Start)
               enddo
            enddo
            IDIIS_Error_Start = IDIIS_Error_Start + 1

            do i=1,IDIISfinal
               RHS(i)=0.0d0
            enddo

            RHS(IDIISfinal+1)=-1.0d0


            call LSOLVE(IDIISfinal+1,quick_method%maxdiisscf+1,B,RHS,W,quick_method%DMCutoff,COEFF,LSOLERR)

            goto 111
         endif


      ! 9) Form a new alpha operator matrix based on
      ! O(new) = [Sum over i] c(i)O(i)
      ! If the solution to step eight failed, skip this step and revert
      ! to a standard scf cycle.

      if (LSOLERR == 0) then
         do J=1,nbasis
            do K=1,nbasis
               OJK=0.d0
               do I=IDIIS_Error_Start, IDIIS_Error_End
                  OJK = OJK + COEFF(I-IDIIS_Error_Start+1) * alloperator(I,K,J)
               enddo
               quick_qm_struct%o(J,K) = OJK
            enddo
         enddo
      endif


      ! 10) Diagonalize the alpha operator matrix to form a new alpha
      ! density matrix.

      ! First you have to transpose this into an orthogonal basis, which
      ! is accomplished by calculating Transpose[X] . O . X.

#ifdef CUDA
      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%o, &
            nbasis, quick_qm_struct%x, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x, &
            nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_qm_struct%o,nbasis)
#else

      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%o,&
               nbasis, quick_qm_struct%x, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x,&
               nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_qm_struct%o,nbasis)
#endif

      ! Now diagonalize the operator matrix.
      call cpu_time(timer_begin%TDiag)
      CALL DIAG(nbasis,quick_qm_struct%o,nbasis,quick_method%DMCutoff,V2,quick_qm_struct%E, &
            quick_qm_struct%idegen,quick_qm_struct%vec, &
            IERROR)
      call cpu_time(timer_end%TDiag)

      step_TDiag=timer_end%TDiag-timer_begin%TDiag

      ! Calculate C = XC' and form a new density matrix.
      ! The C' is from the above diagonalization.  Also, save the previous
      ! Density matrix to check for convergence.

#ifdef CUDA
      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x, &
               nbasis, quick_qm_struct%vec, nbasis, 0.0d0, quick_qm_struct%co,nbasis)
#else
      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x,&
           nbasis, quick_qm_struct%vec, nbasis, 0.0d0, quick_qm_struct%co,nbasis)
#endif

      call CopyDMat(quick_qm_struct%dense,quick_scratch%hold,nbasis) 

      do I=1,nbasis
         do J=1,nbasis
            DENSEJI = 0.d0
            do K=1,quick_molspec%nelec
               DENSEJI = DENSEJI + (quick_qm_struct%co(J,K)*quick_qm_struct%co(I,K))
            enddo
            quick_qm_struct%dense(J,I) = DENSEJI
         enddo
      enddo

      ! Now check for convergence of the alpha matrix.

      PCHANGE=0.d0
      do I=1,nbasis
         do J=1,nbasis
            PCHANGE=MAX(PCHANGE,ABS(quick_qm_struct%dense(J,I)-quick_scratch%hold(J,I)))
         enddo
      enddo

      PRMS = rms(quick_qm_struct%dense,quick_scratch%hold,nbasis)

      ! 11) Form a new BETA operator matrix based on
      ! O(new) = [Sum over i] c(i)O(i)
      ! If the solution to step eight failed, skip this step and revert
      ! to a standard scf cycle.

         if (LSOLERR == 0) then
            do J=1,nbasis
               do K=1,nbasis
                  OJK=0.d0
                  do I=IDIIS_Error_Start, IDIIS_Error_End
                     OJK = OJK + COEFF(I-IDIIS_Error_Start+1) * alloperatorb(I,K,J)
                  enddo
                  quick_qm_struct%ob(J,K) = OJK
               enddo
            enddo
         endif

      ! 8) Diagonalize the beta operator matrix to form a new beta
      ! density matrix.

      ! First you have to transpose this into an orthogonal basis, which
      ! is accomplished by calculating Transpose[X] . O . X.

#ifdef CUDA
      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%ob, &
            nbasis, quick_qm_struct%x, nbasis, 0.0d0, quick_scratch%hold,nbasis)

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x, &
            nbasis, quick_scratch%hold, nbasis, 0.0d0, quick_qm_struct%ob,nbasis)
#else
      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%ob,&
            nbasis, quick_qm_struct%x, nbasis, 0.0d0,quick_scratch%hold,nbasis)

      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x,&
            nbasis, quick_scratch%hold, nbasis, 0.0d0,quick_qm_struct%ob,nbasis)
#endif

      ! Now diagonalize the operator matrix.
      call cpu_time(timer_begin%TDiag)
      CALL DIAG(nbasis,quick_qm_struct%ob,nbasis,quick_method%DMCutoff,V2, &
        quick_qm_struct%EB,quick_qm_struct%idegen,quick_qm_struct%vec, IERROR)
      call cpu_time(timer_end%TDiag)

      step_TDiag=step_TDiag+timer_end%TDiag-timer_begin%TDiag

      ! Calculate C = XC' and form a new density matrix.
      ! The C' is from the above diagonalization.  Also, save the previous
      ! Density matrix to check for convergence.

#ifdef CUDA

      call cublas_DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x, &
            nbasis, quick_qm_struct%vec, nbasis, 0.0d0, quick_qm_struct%cob,nbasis)
#else
      call DGEMM ('n', 'n', nbasis, nbasis, nbasis, 1.0d0, quick_qm_struct%x, &
            nbasis, quick_qm_struct%vec, nbasis, 0.0d0, quick_qm_struct%cob,nbasis)
#endif

      call CopyDMat(quick_qm_struct%denseb,quick_scratch%hold,nbasis)

      do I=1,nbasis
         do J=1,nbasis
            DENSEJI = 0.d0
            do K=1,quick_molspec%nelecb
               DENSEJI = DENSEJI + (quick_qm_struct%cob(J,K)*quick_qm_struct%cob(I,K))
            enddo
            quick_qm_struct%denseb(J,I) = DENSEJI
         enddo
      enddo
      call cpu_time(timer_end%TDII)

      ! Now check for convergence.

      do I=1,nbasis
         do J=1,nbasis
            PCHANGE=MAX(PCHANGE,ABS(quick_qm_struct%denseb(J,I)-quick_scratch%hold(J,I)))
         enddo
      enddo

      PRMS2 = rms(quick_qm_struct%denseb,quick_scratch%hold,nbasis)

      PRMS = MAX(PRMS,PRMS2)

      tmp = quick_method%integralCutoff

      call adjust_cutoff(PRMS,PCHANGE,quick_method)  !from quick_method_module

      !call cpu_time(TEne_begin)

      !call UHFEnergy

      !call cpu_time(TEne_end)

      call cpu_time(timer_end%TSCF)

      ! Add to cumulative times
      timer_cumer%TOp=timer_cumer%TOp+timer_end%TOp-timer_begin%TOp
      timer_cumer%TSCF=timer_cumer%TSCF+timer_end%TSCF-timer_begin%TSCF
      timer_cumer%TDII=timer_cumer%TDII+timer_end%TDII-timer_begin%TDII
      timer_cumer%TDiag=step_TDiag+timer_cumer%TDiag

      current_diis=mod(idiis-1,quick_method%maxdiisscf)
      current_diis=current_diis+1

      write (ioutfile,'(I3,1x)',advance="no") jscf
      if(quick_method%printEnergy)then
         write (ioutfile,'(F16.9,2x)',advance="no") quick_qm_struct%Eel+quick_qm_struct%Ecore
         if (jscf.ne.1) then
            write(ioutFile,'(E12.6,2x)',advance="no") oldEnergy-quick_qm_struct%Eel-quick_qm_struct%Ecore
         else
            write(ioutFile,'(4x,"------",4x)',advance="no")
         endif
         oldEnergy=quick_qm_struct%Eel+quick_qm_struct%Ecore
      endif

      write(*,'(I3,2X,f50.40)') jscf,quick_qm_struct%Eel+quick_qm_struct%Ecore
      write (ioutfile,'(F10.3,4x)',advance="no") timer_end%TSCF-timer_begin%TSCF
      write (ioutfile,'(I2,4x,F8.2,2x,F8.2,2x)',advance="no") current_diis,timer_end%TDII-timer_begin%TDII, &
            timer_end%TOp-timer_begin%TOp
      !write (ioutfile,'(I2,4x,F8.2,2x,F8.2,2x)',advance="no") idiis,timer_end%TDII-timer_begin%TDII, &
      !      step_TOp      
      write (ioutfile,'(F8.2,4x)',advance="no") step_TDiag
      write (ioutfile,'(E10.4,2x)',advance="no") errormax
      write (ioutfile,'(E10.4,2x,E10.4)')  PRMS,PCHANGE

      if (lsolerr /= 0) write (ioutfile,'("DIIS FAILED !!", &
               & " PERFORM NORMAL SCF. (NOT FATAL.)")')

      if (PRMS < quick_method%pmaxrms .and. pchange < quick_method%pmaxrms*100.d0)then

         if (quick_method%printEnergy) then
            write(ioutfile,'(120("-"))')
         else
            write(ioutfile,'(90("-"))')
         endif
         write (ioutfile,'(" REACH CONVERGENCE AFTER ",i3," CYLCES")') jscf
         write (ioutfile,'(" MAX ERROR = ",E12.6,2x," RMS CHANGE = ",E12.6,2x," MAX CHANGE = ",E12.6)') &
               errormax,prms,pchange
         write (ioutfile,*) '-----------------------------------------------'
         if (quick_method%DFT .or. quick_method%SEDFT) then
            write (ioutfile,'("ALPHA ELECTRON DENSITY    =",F16.10)') quick_qm_struct%aelec
            write (ioutfile,'("BETA ELECTRON DENSITY     =",F16.10)') quick_qm_struct%belec
         endif

         if (quick_method%prtgap) then 
                write (ioutfile,'("HOMO-LUMO GAP (EV) =",11x,F12.6)') &
               (quick_qm_struct%E((quick_molspec%nelec/2)+1) - quick_qm_struct%E(quick_molspec%nelec/2))*AU_TO_EV
         endif      

         diisdone=.true.

      endif

      if(jscf >= quick_method%iscf-1) then
         write (ioutfile,'("RAN OUT OF CYCLES.  NO CONVERGENCE.")')
         write (ioutfile,'("PERFORM FINAL NO INTERPOLATION ITERATION")')
         diisdone=.true.
      endif

         diisdone = idiis.gt.MAX_DII_CYCLE_TIME*quick_method%maxdiisscf .or.diisdone

       if((tmp .ne. quick_method%integralCutoff).and. .not.diisdone) then
            write(ioutfile, '(4x, "--------------- 2E-INT CUTOFF CHANGE TO ",E10.4, " -------------")') quick_method%integralCutoff
       endif

      !write (ioutfile,'(/,"SCF CYCLE      = ",I8, &
      !      & "      TIME      = ",F6.2)') &
      !      jscf,T2-T1
      !write (ioutfile,'("DIIS CYCLE     = ",I8, &
      !      & "      MAX ERROR = ",E12.6)') &
      !      idiis,errormax
      !write (ioutfile,'("RMS CHANGE     = ",E12.6, &
      !      & "  MAX CHANGE= ",E12.6)') &
      !      PRMS,PCHANGE

      !if (quick_method%DFT .OR. quick_method%SEDFT) then
      !   write (ioutfile,'("ALPHA ELECTRON DENSITY    =",F16.10)') &
      !         quick_qm_struct%aelec
      !   write (ioutfile,'("BETA ELECTRON DENSITY     =",F16.10)') &
      !         quick_qm_struct%belec
      !endif

      !if (quick_method%prtgap) then
      !   write (ioutfile,'("ALPHA HOMO-LUMO GAP (EV) =", &
      !         & 11x,F12.6)') (quick_qm_struct%E(quick_molspec%nelec+1) - quick_qm_struct%E(quick_molspec%nelec))*27.2116d0
      !   write (ioutfile,'("BETA HOMO-LUMO GAP (EV)  =", &
      !         & 11x,F12.6)') (quick_qm_struct%EB(quick_molspec%nelecb+1) - quick_qm_struct%EB(quick_molspec%nelecb))*27.2116d0
      !endif

      !if (PRMS < quick_method%pmaxrms .and. pchange < quick_method%pmaxrms*100.d0)then
      !   write (ioutfile,' &
      !         & ("PREPARING FOR FINAL NO INTERPOLATION ITERATION")')
      !   diisdone=.true.
      !   elseif(OLDPRMS <= PRMS) then
      !   write (ioutfile,' &
      !         & ("DIIS NOT IMPROVING. RETURN TO TRADITIONAL SCF.")')
      !   diisdone=.true.
      !endif
      !if(jscf >= quick_method%iscf-1) then
      !   write (ioutfile,'("RAN OUT OF CYCLES.  NO CONVERGENCE.")')
      !   write (ioutfile,' &
      !         & ("PERFORM FINAL NO INTERPOLATION ITERATION")')
      !   diisdone=.true.
      !endif
!      diisdone = idiis.eq.quick_method%maxdiisscf.or.diisdone

   enddo

#ifdef CUDA
   if(quick_method%bCUDA) then
      if (quick_method%DFT) then
         call gpu_delete_dft_grid()
      endif
   endif
#endif

   return

end subroutine uelectdiis_new_imp2

