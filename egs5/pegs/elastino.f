!-----------------------------------------------------------------------
!                             SUBROUTINE ELASTINO  
!  Version: 060313-1235
!  Reference:  Based on code developed and provided by F. Salvat
!              for computing GS Mult Scat with PW cross sections
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

      subroutine elastino
C
C  This program computes multiple elastic scattering distributions of
C  electrons and positrons in matter. The DCS is obtained by log-log
C  cubic spline interpolation of atomic DCSs from an elastic scatter-
C  ing database, which was generated by running the code ELSEPA with
C  Desclaux's multiconfiguration Dirac-Fock atomic electron density and
C  a Fermi nuclear charge distribution. For electrons, the exchange
C  potential of Furness and McCarthy was used.
C
C                    Francesc Salvat. Barcelona/Ann Arbor. June, 2001.
C
      implicit none

      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'include/egs5_coefgs.f'
      include 'include/egs5_mscon.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_media.f'

      include 'pegscommons/dcsstr.f'

      double precision Y(NREDA),RA(NREDA),RB(NREDA),RC(NREDA),RD(NREDA)
      double precision Y1(NREDA),Y2(NREDA)

      integer ntb
      PARAMETER (NTB=NFIT1)
      double precision  XGR(NTB),YC(NTB),Y0(NTB),
     1          YY(NREDA),TA(NREDA),TB(NREDA),TC(NREDA),TD(NREDA)

      logical printDiag, longDiag, otherM, uniInt
      integer i, j, n, ik1, iener, ipart, ielec, nleg, ntermm, nterm
      integer ia, didGS
      integer nleg0
      double precision e, s, thdeg,raddeg
      double precision sumt, sum, sump, dsum, xl, xu, xx, st, dmu
      double precision k1start(2), dk1(2), k1, ecut, emax

      double precision DCSEL

      raddeg = 180/3.14159
      longDiag = .false.
      printDiag = .true.
      otherM = .false.
      uniInt = .false.

!  rewrite the pre-amble to the msfit file to note the new values
!  of the characteristic dimension, etc.

      didGs = 1
      write(17,*) nsdcs
      do n = 1, nsdcs
        write(17,5001) (mednam(n,i),i=1,24)
        write(17,*) didGS,charD(n),efrch(n),efrcl(n),
     *              egrdhi(n),egrdlo(n)
      end do
5001  format(' MEDIUM=',24A1)

!  open a new pegs5 listing file for diagnostics
!  ***  Get the K1 range constants and write the ladder

      open(UNIT=26,FILE='pgs5job.GSlst',STATUS='unknown')

      dk1(1) = dexp( (1.d0/(NK1-1.d0)) * dlog(k1maxe/k1mine) )
      dk1(2) = dexp( (1.d0/(NK1-1.d0)) * dlog(k1maxp/k1minp) )
      k1start(1) = k1mine
      k1start(2) = k1minp

      write(17,*) 'Scattering strength ladder for this run'
      write(17,*) 'IK1       K1(e-)         K1(e+)'
      do ik1 = 1, NK1
        write(17,*) ik1, 
     &              k1start(1) * dk1(1) ** (ik1-1),
     &              k1start(2) * dk1(2) ** (ik1-1)
      end do

      if(printDiag) write(26,1001)
1001  format('In Elastino to get Multiple Scattering distribution')

!  ***  Loop over all the materials in the problem
!
      do n = 1, nsdcs

!  ***  load the cross section data for this material

      call dcsload(n,ecut,emax,nleg0)
      call inigrd(ecut,emax,0)

!  ***  Loop over the two particle types
!
      do ipart = 1,2
!
        NLEG = MIN(nleg0,NGT)
!
        IF(ipart.EQ.1) THEN
          ielec = -1
          if(printDiag) WRITE(26,*) 'Projectile: electron.'
        ELSE
          ielec = +1
          if(printDiag) WRITE(26,*) 'Projectile: positron.'
        ENDIF
!
!  Loop over the energy steps 
!
        do iener = 1,nmscate
!
          e = mscate(iener)
          if(e .gt. 1.d8) go to 1005
C
C  ****  Now, we generate a table of the DCS for the selected projectile
C        and energy, by log-log cubic spline interpolation in E.
C
          CALL DCSTAB(E,IELEC)
C
C  ****  Here we compute the Goudsmit-Saunderson transport coefficients,
C        from which we can obtain the angular distribution after a given
C        path length (in units of the mean free path).
C
          CALL GSCOEF(NLEG)
C
C  store the scattering power (units on tcs1 are cm^2 from GSCOEF)
C
          if(printDiag) then
          WRITE(26,*) 'Kinetic energy (MeV) =',E*1.d-6
          WRITE(26,*) ' Total cross section (cm^-1) =',CS*atomd(n)
          WRITE(26,*) ' 1st transport cross section =',TCS1*atomd(n)
          if(longDiag) WRITE(26,*) ' 2nd transport cross section =',TCS2
          endif

!  loop over all the values of mubar that we're likely to
!  see.  This is actually a loop over a range of distances...

          do ik1 = 1, NK1

!  ****  get s from k1 

            k1 = k1start(ipart) * dk1(ipart) ** (ik1-1)

!            CALL GSTEST  ! Numerical check not available this version
!  ****  We are now ready to evaluate GS distributions for any path
!        length.  Get the pathlength (in mfp's) at this energy

            s = k1 * cs / tcs1
!
            if(printDiag) then
            WRITE(26,*) ' For this step, K1 =',k1
            WRITE(26,*) '  At this K1, path length (mfps)=',S
            WRITE(26,*) '  At this K1, path length (cm) =',
     &                                                     S/CS/atomd(n)
            endif
            if(longDiag) then
              WRITE(26,*) '# '
              WRITE(26,*) '# Differential cross section'
              WRITE(26,*) '# theta(deg)       mu        dcs(cm^2/sr)'
              DO IA=1,NREDA
                WRITE(26,'(1P,4E14.6)') TH(IA),XMU(IA),DCSI(IA)
              ENDDO
              WRITE(26,*) '# Goudsmit-Saunderson distribution, PDF(mu)'
              IF(IELEC.EQ.1) THEN
                WRITE(26,*) '# Projectile: positron'
              ELSE
                WRITE(26,*) '# Projectile: electron'
              ENDIF
              WRITE(26,*) '# Number of terms in Legendre series =',NLEG
              WRITE(26,*) '# '
              WRITE(26,*) '# theta(deg)       mu           PDF(mu)',
     1       '       S*DCS       NTERM '
            endif

            nterm = nleg
            CALL GSDIST(S,0.0D0,Y(1),NTERM)
            NTERMM=NTERM
            DO I=1,NREDA
              CALL GSDIST(S,XMU(I),Y(I),NTERM)
              if(longDiag) then
                THDEG=ACOS(1.0D0-2.0D0*XMU(I))*RADDEG
                WRITE(26,'(1P,4E14.6,I7)') THDEG,XMU(I),Y(I),
     1         2.0D0*S*DCSEL(XMU(I))/CS0,NTERM
              endif
              NTERMM=MAX(NTERMM,NTERM)
            ENDDO
    2       CONTINUE
            if(printDiag) WRITE(26,*) ' Number of terms used: ',NTERMM
            !  this is a temp fix until single scattering
            if(nterm.lt.0) nterm = ntermm
C
C  ****  Finally, we check that the normalization is correct.
C
C  ... we first compute the integral of the continuous distribution,
            CALL SPLINE(XMU,Y,RA,RB,RC,RD,0.0D0,0.0D0,NREDA)
            CALL INTEG(XMU,RA,RB,RC,RD,0.0D0,XMU(NREDA),SUMP,NREDA)
C  ... and we add the probability of no scattering, EXP(-S).
            if(s.lt.100)  then
              SUM=SUMP+EXP(-S)
            else
              sum = sump
            endif
C
            if(printDiag) then
            WRITE(26,*) ' Normalization =',SUM
            WRITE(26,*) ' No-scatter probability = ',sum-sump
            endif
            probns(ipart,iener,ik1) = sum-sump
C
C  ****  Cumulative probability distribution and moments.
C        Equiprobable intervals.
C
            DSUM=SUMP/DBLE(NBFIT)
            XGR(1)=0.0D0
            DO I=2,NBFIT
              SUMT=(I-1)*DSUM
              XL=XGR(I-1)
              XU=1.0D0
1234          XX=0.5D0*(XL+XU)
              CALL INTEG(XMU,RA,RB,RC,RD,0.0D0,XX,ST,NREDA)
              IF(ST.GT.SUMT) THEN
                XU=XX
              ELSE
                XL=XX
              ENDIF
              IF(ABS(XU-XL).GT.1.0D-6*XX) GO TO 1234
              XGR(I)=XX
            ENDDO
C
C  Now get extra bins in the last angle, set uniformly
C
            dmu = (1.0d0 - xgr(NBFIT)) / NEXFIT
            do i=1,NEXFIT
              j = NBFIT + i
              xgr(j) = xgr(NBFIT) + i * dmu
            end do
C
C now get the cum dist for all angles
C
            DO I=1,NREDA
              YY(I)=Y(I)/SUM
            ENDDO
            CALL SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
            nterm = nleg
            DO I=1,NTB
              CALL GSDIST(S,XGR(I),YC(I),NTERM)
              if(i.le.NBFIT) then
                y0(i) = (i-1)*dsum/sum
              else 
                CALL INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y0(I),NREDA)
              endif
            ENDDO
            
            call fitms(xgr,yc,y0,ntb,ipart,iener,ik1) 

           ! other moments...
            if(otherM) then
              DO I=1,NREDA
                YY(I)=YY(I)*XMU(I)
              ENDDO
              CALL SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
              DO I=1,NTB
                CALL INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y1(I),NREDA)
              ENDDO

              DO I=1,NREDA
                YY(I)=YY(I)*XMU(I)
              ENDDO
              CALL SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
              DO I=1,NTB
                CALL INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y2(I),NREDA)
              ENDDO

              if(printDiag) then
              WRITE(26,*) 
     1 '# Cumulative Goudsmit-Saunderson distributions'
              IF(IELEC.EQ.1) THEN
                WRITE(26,*) '# Projectile: positron'
              ELSE
                WRITE(26,*) '# Projectile: electron'
              ENDIF
              WRITE(26,*) '# Kinetic energy (eV) =',E
              WRITE(26,*) '# Path length/mfp =',S
              WRITE(26,*) '# Number of terms in Legendre series =',NLEG
              if(longDiag) then
              WRITE(26,*) '# '
              WRITE(26,3001)
 3001 FORMAT(' # theta(deg)       mu           PDF(mu)',
     1  '        CDF         <MU*PDF>     <MU*MU*PDF>')
              DO I=1,NTB
                THDEG=ACOS(1.0D0-2.0D0*XGR(I))*RADDEG
                WRITE(26,'(1P,6E24.16)') 
     1                 THDEG,XGR(I),YC(I),Y0(I),Y1(I),Y2(I)
              ENDDO
              endif
              endif
            endif
C
C  ****  Cumulative probability distribution and moments
C        for Uniform intervals.
            if(uniInt) then

              DMU=1.0D0/DBLE(NTB-1)
              DO I=1,NTB
                XGR(I)=(I-1)*DMU
              ENDDO

              DO I=1,NREDA
                YY(I)=Y(I)/SUM
              ENDDO
              CALL SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
              DO I=1,NTB
                CALL GSDIST(S,XGR(I),YC(I),NTERM)
                CALL INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y0(I),NREDA)
              ENDDO

              call fitms(xgr,yc,y0,ntb,ipart,iener,ik1) 

C      other moments...
              if(otherM) then
                DO I=1,NREDA
                  YY(I)=YY(I)*XMU(I)
                ENDDO
                CALL SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
                DO I=1,NTB
                  CALL INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y1(I),NREDA)
                ENDDO

                DO I=1,NREDA
                  YY(I)=YY(I)*XMU(I)
                ENDDO
                CALL SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
                DO I=1,NTB
                  CALL INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y2(I),NREDA)
                ENDDO
C
                if(printDiag) then
                WRITE(26,*)
     1 '# Cumulative Goudsmit-Saunderson distributions'
                IF(IELEC.EQ.1) THEN
                  WRITE(26,*) '# Projectile: positron'
                ELSE
                  WRITE(26,*) '# Projectile: electron'
                ENDIF
                WRITE(26,*) '# Kinetic energy (eV) =',E
                WRITE(26,*) '# Path length/mfp =',S
                WRITE(26,*) '# Number of terms in Legendre series=',NLEG
                if(longDiag) then
                WRITE(26,*) '# '
                WRITE(26,3002)
 3002 FORMAT(' # theta(deg)       mu           PDF(mu)',
     1   '        CDF         <MU*PDF>     <MU*MU*PDF>')
                DO I=1,NTB
                  THDEG=ACOS(1.0D0-2.0D0*XGR(I))*RADDEG
                  WRITE(26,'(1P,6E24.16)') 
     1                 THDEG,XGR(I),YC(I),Y0(I),Y1(I),Y2(I)
                ENDDO
                endif
                endif
              endif        !  get other moments
            endif          !  use uniform ints

          end do	   !  end path length loop
        end do	           !  end energy loop
      end do               !  end particle types loop

1005  call wmsfit(k1start,dk1)

1010  end do               !  end material loop

      close(26)

      return
      END
