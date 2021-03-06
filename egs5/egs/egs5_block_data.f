!---------------------------egs5_block_data.f---------------------------
! Version: 060802-1335
! Reference: SLAC-R-730/KEK-2005-8
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

      block data

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_elecin.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
      include 'include/egs5_media.f'
      include 'include/egs5_mults.f'
      include 'include/egs5_thresh.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/randomm.f'

      character*4 media1(24)                            ! Local variable
      equivalence(media1(1),media(1,1))

      data                                              ! common/ELECIN/
     * ekelim/0./,
     * icomp/1/

      data                                              ! common/EPCONT/
     * iausfl/5*1,MXAUSM5*0/,
     * rhof/1.0/

      data                                               ! common/MEDIA/
     * nmed/1/,
     * media1/'N','A','I',' ',' ',' ',' ',' ',' ',' ',' ',' ',
     *        ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/,
     * charD/MXMED*0.d0/,
     * iraylm/MXMED*0/,
     * incohm/MXMED*0/,
     * iprofm/MXMED*0/,
     * impacm/MXMED*0/,
     * useGSD/MXMED*0/

      data                                               ! common/MULTS/
     * NG21/  7/,B0G21/ 2.0000E+00/,B1G21/ 5.0000E+00/,
     * G210( 1),G211( 1),G212( 1)/-9.9140E-04, 2.7672E+00,-1.1544E+00/,
     * G210( 2),G211( 2),G212( 2)/-9.9140E-04, 2.7672E+00,-1.1544E+00/,
     * G210( 3),G211( 3),G212( 3)/-7.1017E-02, 3.4941E+00,-3.0773E+00/,
     * G210( 4),G211( 4),G212( 4)/-7.3556E-02, 3.5487E+00,-3.1989E+00/,
     * G210( 5),G211( 5),G212( 5)/ 3.6658E-01, 2.1162E+00,-2.0311E+00/,
     * G210( 6),G211( 6),G212( 6)/ 1.4498E+00,-5.9717E-01,-3.2951E-01/,
     * G210( 7),G211( 7),G212( 7)/ 1.4498E+00,-5.9717E-01,-3.2951E-01/
      data
     * NG22/  8/,B0G22/ 2.0000E+00/,B1G22/ 6.0000E+00/,
     * G220( 1),G221( 1),G222( 1)/-5.2593E-04, 1.4285E+00,-1.2670E+00/,
     * G220( 2),G221( 2),G222( 2)/-5.2593E-04, 1.4285E+00,-1.2670E+00/,
     * G220( 3),G221( 3),G222( 3)/-6.4819E-02, 2.2033E+00,-3.6399E+00/,
     * G220( 4),G221( 4),G222( 4)/ 3.7427E-02, 1.6630E+00,-2.9362E+00/,
     * G220( 5),G221( 5),G222( 5)/ 6.1955E-01,-6.2713E-01,-6.7859E-01/,
     * G220( 6),G221( 6),G222( 6)/ 1.7584E+00,-4.0390E+00, 1.8810E+00/,
     * G220( 7),G221( 7),G222( 7)/ 2.5694E+00,-6.0484E+00, 3.1256E+00/,
     * G220( 8),G221( 8),G222( 8)/ 2.5694E+00,-6.0484E+00, 3.1256E+00/
      data
     * NG31/ 11/,B0G31/ 2.0000E+00/,B1G31/ 9.0000E+00/,
     * G310( 1),G311( 1),G312( 1)/ 4.9437E-01, 1.9124E-02, 1.8375E+00/,
     * G310( 2),G311( 2),G312( 2)/ 4.9437E-01, 1.9124E-02, 1.8375E+00/,
     * G310( 3),G311( 3),G312( 3)/ 5.3251E-01,-6.1555E-01, 4.5595E+00/,
     * G310( 4),G311( 4),G312( 4)/ 6.6810E-01,-2.2056E+00, 8.9293E+00/,
     * G310( 5),G311( 5),G312( 5)/-3.8262E+00, 2.5528E+01,-3.3862E+01/,
     * G310( 6),G311( 6),G312( 6)/ 4.2335E+00,-1.0604E+01, 6.6702E+00/,
     * G310( 7),G311( 7),G312( 7)/ 5.0694E+00,-1.4208E+01, 1.0456E+01/,
     * G310( 8),G311( 8),G312( 8)/ 1.4563E+00,-3.3275E+00, 2.2601E+00/,
     * G310( 9),G311( 9),G312( 9)/-3.2852E-01, 1.2938E+00,-7.3254E-01/,
     * G310(10),G311(10),G312(10)/-2.2489E-01, 1.0713E+00,-6.1358E-01/,
     * G310(11),G311(11),G312(11)/-2.2489E-01, 1.0713E+00,-6.1358E-01/
      data
     * NG32/ 25/,B0G32/ 2.0000E+00/,B1G32/ 2.3000E+01/,
     * G320( 1),G321( 1),G322( 1)/ 2.9907E-05, 4.7318E-01, 6.5921E-01/,
     * G320( 2),G321( 2),G322( 2)/ 2.9907E-05, 4.7318E-01, 6.5921E-01/,
     * G320( 3),G321( 3),G322( 3)/ 2.5820E-03, 3.5853E-01, 1.9776E+00/,
     * G320( 4),G321( 4),G322( 4)/-5.3270E-03, 4.9418E-01, 1.4528E+00/,
     * G320( 5),G321( 5),G322( 5)/-6.6341E-02, 1.4422E+00,-2.2407E+00/,
     * G320( 6),G321( 6),G322( 6)/-3.6027E-01, 4.7190E+00,-1.1380E+01/,
     * G320( 7),G321( 7),G322( 7)/-2.7953E+00, 2.6694E+01,-6.0986E+01/,
     * G320( 8),G321( 8),G322( 8)/-3.6091E+00, 3.4125E+01,-7.7512E+01/,
     * G320( 9),G321( 9),G322( 9)/ 1.2491E+01,-7.1103E+01, 9.4496E+01/,
     * G320(10),G321(10),G322(10)/ 1.9637E+01,-1.1371E+02, 1.5794E+02/,
     * G320(11),G321(11),G322(11)/ 2.1692E+00,-2.5019E+01, 4.5340E+01/,
     * G320(12),G321(12),G322(12)/-1.6682E+01, 6.2067E+01,-5.5257E+01/,
     * G320(13),G321(13),G322(13)/-2.1539E+01, 8.2651E+01,-7.7065E+01/,
     * G320(14),G321(14),G322(14)/-1.4344E+01, 5.5193E+01,-5.0867E+01/,
     * G320(15),G321(15),G322(15)/-5.4990E+00, 2.3874E+01,-2.3140E+01/,
     * G320(16),G321(16),G322(16)/ 3.1029E+00,-4.4708E+00, 2.1318E-01/,
     * G320(17),G321(17),G322(17)/ 6.0961E+00,-1.3670E+01, 7.2823E+00/,
     * G320(18),G321(18),G322(18)/ 8.6179E+00,-2.0950E+01, 1.2536E+01/,
     * G320(19),G321(19),G322(19)/ 7.5064E+00,-1.7956E+01, 1.0520E+01/,
     * G320(20),G321(20),G322(20)/ 5.9838E+00,-1.4065E+01, 8.0342E+00/,
     * G320(21),G321(21),G322(21)/ 4.4959E+00,-1.0456E+01, 5.8462E+00/,
     * G320(22),G321(22),G322(22)/ 3.2847E+00,-7.6709E+00, 4.2445E+00/,
     * G320(23),G321(23),G322(23)/ 1.9514E+00,-4.7505E+00, 2.6452E+00/,
     * G320(24),G321(24),G322(24)/ 4.8808E-01,-1.6910E+00, 1.0459E+00/,
     * G320(25),G321(25),G322(25)/ 4.8808E-01,-1.6910E+00, 1.0459E+00/
      data
     * NBGB/  8/,B0BGB/ 1.5714E+00/,B1BGB/ 2.1429E-01/,
     * BGB0( 1),BGB1( 1),BGB2( 1)/-1.0724E+00, 2.8203E+00,-3.5669E-01/,
     * BGB0( 2),BGB1( 2),BGB2( 2)/ 3.7136E-01, 1.4560E+00,-2.8072E-02/,
     * BGB0( 3),BGB1( 3),BGB2( 3)/ 1.1396E+00, 1.1910E+00,-5.2070E-03/,
     * BGB0( 4),BGB1( 4),BGB2( 4)/ 1.4908E+00, 1.1267E+00,-2.2565E-03/,
     * BGB0( 5),BGB1( 5),BGB2( 5)/ 1.7342E+00, 1.0958E+00,-1.2705E-03/,
     * BGB0( 6),BGB1( 6),BGB2( 6)/ 1.9233E+00, 1.0773E+00,-8.1806E-04/,
     * BGB0( 7),BGB1( 7),BGB2( 7)/ 2.0791E+00, 1.0649E+00,-5.7197E-04/,
     * BGB0( 8),BGB1( 8),BGB2( 8)/ 2.0791E+00, 1.0649E+00,-5.7197E-04/
 
      data                                              ! common/THRESH/
     * RMT2/1.021997804/            ! Two times electron rest mass (MeV)
     * RMSQ/0.261119878/           ! Electron rest mass squared (MeV**2)

      data                                              ! common/UPHIOT/
     * PI/3.1415926535897932d0/,                                    ! Pi
     * TWOPI/6.2831853071795864d0/,                         ! 2 times Pi
     * PI5D2/7.853981634/                   ! Five times Pi divided by 2

      data                                              ! common/USEFUL/
     * RM/0.510998902/                        ! Electron rest mass (MeV)

      data                                             ! common/RLUXDAT/
     * twom24/1./,
     * ndskip/0,24,73,199,365/,
     * luxlev/1/,
     * inseed/0/,
     * kount/0/,
     * mkount/0/,
     * isdext/25*0/,
     * rluxset/.false./

      end
!---------------------last line of egs5_block_data.f--------------------
