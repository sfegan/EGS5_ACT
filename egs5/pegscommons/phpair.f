!----------------------------phpair.f-----------------------------------
!  Version: 051219-1435
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

! ----------------------------------------
! PEGS5 common file for PHPAIR
! ----------------------------------------

      integer nphe
      double precision phe,phd,ekedge,pre,prd
      COMMON/PHPAIR/PHE(61,100),PHD(61,100),EKEDGE(100),PRE(17),
     &              PRD(17,100),NPHE(100)

!--------------------------last line------------------------------------
