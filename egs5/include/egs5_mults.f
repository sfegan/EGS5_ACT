!-----------------------------egs5_mults.f------------------------------
! Version: 051219-1435
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

      common/MULTS/                                ! Multiple scattering
     * B0G21,B1G21,G210(7),G211(7),G212(7),
     * B0G22,B1G22,G220(8),G221(8),G222(8),
     * B0G31,B1G31,G310(11),G311(11),G312(11),
     * B0G32,B1G32,G320(25),G321(25),G322(25),
     * B0BGB,B1BGB,BGB0(8),BGB1(8),BGB2(8),
     * NG21,NG22,NG31,NG32,NBGB

      real*8
     * B0G21,B1G21,G210,G211,G212,
     * B0G22,B1G22,G220,G221,G222,
     * B0G31,B1G31,G310,G311,G312,
     * B0G32,B1G32,G320,G321,G322,
     * B0BGB,B1BGB,BGB0,BGB1,BGB2

      integer
     * NG21,NG22,NG31,NG32,NBGB

!------------------------last line of egs5_mults.f----------------------
