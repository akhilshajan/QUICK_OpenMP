! Ed Brothers. November 27, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine uhfenergy
    use allmod
    implicit double precision(a-h,o-z)

! The purpose of this subroutine is to calculate the total energy
! of a molecule at the unrestricted hartree fock level given an optimized
! density matrix.


    quick_qm_struct%Eel=0.d0
    DO Ibas=1,nbasis
        DO Icon=1,ncontract(Ibas)
            DO Jcon=1,ncontract(Ibas)

            ! Kinetic energy.

                quick_qm_struct%Eel=quick_qm_struct%Eel+(quick_qm_struct%DENSE(Ibas,Ibas)+quick_qm_struct%DENSEB(Ibas,Ibas))* &
                dcoeff(Jcon,Ibas)*dcoeff(Icon,Ibas)* &
                ekinetic(aexp(Jcon,Ibas),aexp(Icon,Ibas), &
                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                xyz(1,quick_basis%ncenter(Ibas)),xyz(2,quick_basis%ncenter(Ibas)), &
                xyz(3,quick_basis%ncenter(Ibas)),xyz(1,quick_basis%ncenter(Ibas)), &
                xyz(2,quick_basis%ncenter(Ibas)),xyz(3,quick_basis%ncenter(Ibas)))

!                print*,eel
!                print*,Eel,dcoeff(Jcon,Ibas),ekinetic(aexp(Jcon,Ibas),aexp(Icon,Ibas), &
!                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                xyz(1,ncenter(Ibas)),xyz(2,ncenter(Ibas)), &
!                xyz(3,ncenter(Ibas)),xyz(1,ncenter(Ibas)), &
!                xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

            ! Nuclear attraction.

                DO iatom = 1,natom
                    quick_qm_struct%Eel=quick_qm_struct%Eel+(quick_qm_struct%DENSE(Ibas,Ibas)+ &
                    quick_qm_struct%DENSEB(Ibas,Ibas))*dcoeff(Jcon,Ibas)*dcoeff(Icon,Ibas)* &
                    attraction(aexp(Jcon,Ibas),aexp(Icon,Ibas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,quick_basis%ncenter(Ibas)),xyz(2,quick_basis%ncenter(Ibas)), &
                    xyz(3,quick_basis%ncenter(Ibas)),xyz(1,quick_basis%ncenter(Ibas)), &
                    xyz(2,quick_basis%ncenter(Ibas)),xyz(3,quick_basis%ncenter(Ibas)), &
                    xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
                    quick_molspec%chg(iatom))
                ENDDO
            ENDDO
        ENDDO
    ENDDO

    DO Ibas=1,nbasis
        DO Jbas=Ibas+1,nbasis
            DO Icon=1,ncontract(ibas)
                DO Jcon=1,ncontract(jbas)

                ! Kinetic energy.

                    quick_qm_struct%Eel=quick_qm_struct%Eel+(quick_qm_struct%DENSE(Jbas,Ibas)+ &
                    quick_qm_struct%DENSEB(Jbas,Ibas))* &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    2.d0*ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,quick_basis%ncenter(Jbas)),xyz(2,quick_basis%ncenter(Jbas)), &
                    xyz(3,quick_basis%ncenter(Jbas)),xyz(1,quick_basis%ncenter(Ibas)), &
                    xyz(2,quick_basis%ncenter(Ibas)),xyz(3,quick_basis%ncenter(Ibas)))

                ! Nuclear attraction.

                    DO iatom = 1,natom
                        quick_qm_struct%Eel=quick_qm_struct%Eel+(quick_qm_struct%DENSE(Jbas,Ibas)+ &
                        quick_qm_struct%DENSEB(Jbas,Ibas))*dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                        2.d0*attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                        xyz(1,quick_basis%ncenter(Jbas)),xyz(2,quick_basis%ncenter(Jbas)), &
                        xyz(3,quick_basis%ncenter(Jbas)),xyz(1,quick_basis%ncenter(Ibas)), &
                        xyz(2,quick_basis%ncenter(Ibas)),xyz(3,quick_basis%ncenter(Ibas)), &
                        xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
                        quick_molspec%chg(iatom))
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO

!
! Alessandro GENONI 03/21/2007
! Computation of the ECP contribution to the Electronic Energy
!
    if (quick_method%ecp) then
      !call ecp_energy
    end if

!
! The previous two terms are the one electron part of the Fock matrix.
! The next two terms define the two electron part.
    DO I=1,nbasis
    ! Set some variables to reduce access time for some of the more
    ! used quantities.

        xI = xyz(1,quick_basis%ncenter(I))
        yI = xyz(2,quick_basis%ncenter(I))
        zI = xyz(3,quick_basis%ncenter(I))
        itype1I=itype(1,I)
        itype2I=itype(2,I)
        itype3I=itype(3,I)
        DENSEII=quick_qm_struct%DENSE(I,I) + quick_qm_struct%DENSEB(I,I)
        DENSEIIA=quick_qm_struct%DENSE(I,I)
        DENSEIIB=quick_qm_struct%DENSEB(I,I)

    ! DO all the (ii|ii) integrals.
        Ibas=I
        Jbas=I
        IIbas=I
        JJbas=I
        repint=0.d0
        DO Icon=1,ncontract(ibas)
            DO Jcon=1,ncontract(jbas)
                DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                        repint = repint+ &
                        dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                        *dcoeff(IIcon,IIbas)*dcoeff(JJcon,JJbas)* &
                        (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                        aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                        itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                        itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                        xI,yI,zI,xI,yI,zI,xI,yI,zI,xI,yI,zI))
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        quick_qm_struct%Eel = quick_qm_struct%Eel+.5d0*(DENSEII*DENSEII &
        -DENSEIIA*DENSEIIA-DENSEIIB*DENSEIIB)*repint

        DO J=I+1,nbasis
        ! Set some variables to reduce access time for some of the more
        ! used quantities. (AGAIN)

            xJ = xyz(1,quick_basis%ncenter(J))
            yJ = xyz(2,quick_basis%ncenter(J))
            zJ = xyz(3,quick_basis%ncenter(J))
            itype1J=itype(1,J)
            itype2J=itype(2,J)
            itype3J=itype(3,J)
            DENSEJI=quick_qm_struct%DENSE(J,I)+quick_qm_struct%DENSEB(J,I)
            DENSEJJ=quick_qm_struct%DENSE(J,J)+quick_qm_struct%DENSEB(J,J)
            DENSEJIA=quick_qm_struct%DENSE(J,I)
            DENSEJJA=quick_qm_struct%DENSE(J,J)
            DENSEJIB=quick_qm_struct%DENSEB(J,I)
            DENSEJJB=quick_qm_struct%DENSEB(J,J)

        ! Find  all the (ii|jj) integrals.
            Ibas=I
            Jbas=I
            IIbas=J
            JJbas=J
            repint=0.d0
            DO Icon=1,ncontract(ibas)
                DO Jcon=1,ncontract(jbas)
                    DO IIcon=1,ncontract(iibas)
                        DO JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ))
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            quick_qm_struct%Eel = quick_qm_struct%Eel+(DENSEII*DENSEJJ &
            -DENSEJIA*DENSEJIA-DENSEJIB*DENSEJIB)*repint

        ! Find  all the (ij|jj) integrals.
            Ibas=I
            Jbas=J
            IIbas=J
            JJbas=J
            repint=0.d0
            DO Icon=1,ncontract(ibas)
                DO Jcon=1,ncontract(jbas)
                    DO IIcon=1,ncontract(iibas)
                        DO JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ,xJ,yJ,zJ))

                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            quick_qm_struct%Eel = quick_qm_struct%Eel + 2.d0*(DENSEJI*DENSEJJ &
            -DENSEJIA*DENSEJJA-DENSEJIB*DENSEJJB)*repint

        ! Find  all the (ii|ij) integrals.
            Ibas=I
            Jbas=I
            iiBAS=i
            JJbas=J
            repint=0.d0
            DO Icon=1,ncontract(ibas)
                DO Jcon=1,ncontract(jbas)
                    DO IIcon=1,ncontract(iibas)
                        DO JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xI,yI,zI,xI,yI,zI,xJ,yJ,zJ))

                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            quick_qm_struct%Eel = quick_qm_struct%Eel+2.d0*(DENSEJI*DENSEII &
            -DENSEJIA*DENSEIIA-DENSEJIB*DENSEIIB)*repint



        ! Find all the (ij|ij) integrals
            Ibas=I
            Jbas=J
            IIbas=I
            JJbas=J
            repint=0.d0
            DO Icon=1,ncontract(ibas)
                DO Jcon=1,ncontract(jbas)
                    DO IIcon=1,ncontract(iibas)
                        DO JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xJ,yJ,zJ))

                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            quick_qm_struct%Eel = quick_qm_struct%Eel + (2.d0*DENSEJI*DENSEJI &
            -DENSEJJA*DENSEIIA-DENSEJJB*DENSEIIB &
            -DENSEJIA*DENSEJIA-DENSEJIB*DENSEJIB)*repint

            DO K=J+1,nbasis
            ! Set some variables to reduce access time for some of the more
            ! used quantities. (AGAIN)

                xK = xyz(1,quick_basis%ncenter(K))
                yK = xyz(2,quick_basis%ncenter(K))
                zK = xyz(3,quick_basis%ncenter(K))
                itype1K=itype(1,K)
                itype2K=itype(2,K)
                itype3K=itype(3,K)
                DENSEKI=quick_qm_struct%DENSE(K,I)+quick_qm_struct%DENSEB(K,I)
                DENSEKJ=quick_qm_struct%DENSE(K,J)+quick_qm_struct%DENSEB(K,J)
                DENSEKK=quick_qm_struct%DENSE(K,K)+quick_qm_struct%DENSEB(K,K)
                DENSEKIA=quick_qm_struct%DENSE(K,I)
                DENSEKJA=quick_qm_struct%DENSE(K,J)
                DENSEKKA=quick_qm_struct%DENSE(K,K)
                DENSEKIB=quick_qm_struct%DENSEB(K,I)
                DENSEKJB=quick_qm_struct%DENSEB(K,J)
                DENSEKKB=quick_qm_struct%DENSEB(K,K)

            ! Find all the (ij|ik) integrals where j>i,k>j
                Ibas=I
                Jbas=J
                IIbas=I
                JJbas=K
                repint=0.d0
                DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                        DO IIcon=1,ncontract(iibas)
                            DO JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xK,yK,zK))

                            ENDDO
                        ENDDO
                    ENDDO
                ENDDO
                quick_qm_struct%Eel = quick_qm_struct%Eel + 2.d0*(2.d0*DENSEJI*DENSEKI &
                - DENSEIIA*DENSEKJA-DENSEIIB*DENSEKJB &
                - DENSEJIA*DENSEKIA-DENSEJIB*DENSEKIB)*repint

            ! Find all the (ij|kk) integrals where j>i, k>j.
                Ibas=I
                Jbas=J
                IIbas=K
                JJbas=K
                repint=0.d0
                DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                        DO IIcon=1,ncontract(iibas)
                            DO JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                itype1K,itype2K,itype3K,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xK,yK,zK))

                            ENDDO
                        ENDDO
                    ENDDO
                ENDDO
                quick_qm_struct%Eel = quick_qm_struct%Eel +2.d0*(DENSEJI*DENSEKK &
                - DENSEKIA*DENSEKJA-DENSEKIB*DENSEKJB)*repint

            ! Find all the (ik|jj) integrals where j>i, k>j.
                Ibas=I
                Jbas=K
                IIbas=J
                JJbas=J
                repint=0.d0
                DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                        DO IIcon=1,ncontract(iibas)
                            DO JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                                itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                                xI,yI,zI,xK,yK,zK,xJ,yJ,zJ,xJ,yJ,zJ))

                            ENDDO
                        ENDDO
                    ENDDO
                ENDDO
                quick_qm_struct%Eel = quick_qm_struct%Eel + 2.d0*(DENSEKI*DENSEJJ &
                -DENSEKJA*DENSEJIA-DENSEKJB*DENSEJIB)*repint

            ! Find all the (ii|jk) integrals where j>i, k>j.
                Ibas=I
                Jbas=I
                IIbas=J
                JJbas=K
                repint=0.d0
                DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                        DO IIcon=1,ncontract(iibas)
                            DO JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                                itype1J,itype2J,itype3J,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xK,yK,zK))

                            ENDDO
                        ENDDO
                    ENDDO
                ENDDO
                quick_qm_struct%Eel = quick_qm_struct%Eel + 2.d0*(DENSEKJ*DENSEII &
                -DENSEJIA*DENSEKIA-DENSEJIB*DENSEKIB)*repint

            ENDDO

            DO K=I+1,nbasis-1
                xK = xyz(1,quick_basis%ncenter(K))
                yK = xyz(2,quick_basis%ncenter(K))
                zK = xyz(3,quick_basis%ncenter(K))
                itype1K=itype(1,K)
                itype2K=itype(2,K)
                itype3K=itype(3,K)
                DENSEKI=quick_qm_struct%DENSE(K,I)+quick_qm_struct%DENSEB(K,I)
                DENSEKJ=quick_qm_struct%DENSE(K,J)+quick_qm_struct%DENSEB(K,J)
                DENSEKK=quick_qm_struct%DENSE(K,K)+quick_qm_struct%DENSEB(K,K)
                DENSEKIA=quick_qm_struct%DENSE(K,I)
                DENSEKJA=quick_qm_struct%DENSE(K,J)
                DENSEKKA=quick_qm_struct%DENSE(K,K)
                DENSEKIB=quick_qm_struct%DENSEB(K,I)
                DENSEKJB=quick_qm_struct%DENSEB(K,J)
                DENSEKKB=quick_qm_struct%DENSEB(K,K)

                DO L=K+1,nbasis
                    xL = xyz(1,quick_basis%ncenter(L))
                    yL = xyz(2,quick_basis%ncenter(L))
                    zL = xyz(3,quick_basis%ncenter(L))
                    itype1L=itype(1,L)
                    itype2L=itype(2,L)
                    itype3L=itype(3,L)
                    DENSELJ=quick_qm_struct%DENSE(L,J)+quick_qm_struct%DENSEB(L,J)
                    DENSELI=quick_qm_struct%DENSE(L,I)+quick_qm_struct%DENSEB(L,I)
                    DENSELK=quick_qm_struct%DENSE(L,K)+quick_qm_struct%DENSEB(L,K)
                    DENSELJA=quick_qm_struct%DENSE(L,J)
                    DENSELIA=quick_qm_struct%DENSE(L,I)
                    DENSELKA=quick_qm_struct%DENSE(L,K)
                    DENSELJB=quick_qm_struct%DENSEB(L,J)
                    DENSELIB=quick_qm_struct%DENSEB(L,I)
                    DENSELKB=quick_qm_struct%DENSEB(L,K)

                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    Ibas=I
                    Jbas=J
                    IIbas=K
                    JJbas=L
                    repint=0.d0
                    DO Icon=1,ncontract(ibas)
                        DO Jcon=1,ncontract(jbas)
                            DO IIcon=1,ncontract(iibas)
                                DO JJcon=1,ncontract(jjbas)
                                    repint = repint+ &
                                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                    *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                    (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                    aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                    itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                    itype1K,itype2K,itype3K,itype1L,itype2L,itype3L, &
                                    xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xL,yL,zL))


                                ENDDO
                            ENDDO
                        ENDDO
                    ENDDO
                    quick_qm_struct%Eel = quick_qm_struct%Eel + 2.d0*(2.d0*DENSEJI*DENSELK &
                    -DENSEKIA*DENSELJA-DENSEKIB*DENSELJB &
                    -DENSELIA*DENSEKJA-DENSELIB*DENSEKJB)*repint
                ENDDO
            ENDDO
        ENDDO
    ENDDO

! Xiao tests
!   Do i=1,nbasis
!     Do j=1,nbasis
!       print*,i,j,DENSE(i,j)
!     enddo
!   enddo


    END subroutine uhfenergy

