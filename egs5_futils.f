      subroutine fort_open_old(theunit,thefile)
        integer theunit
        character*(*) thefile
        open(UNIT=theunit,FILE=thefile,STATUS='old')
      end

      subroutine fort_open_new(theunit,thefile)
        integer theunit
        character*(*) thefile
        open(UNIT=theunit,FILE=thefile,STATUS='new')
      end

      subroutine fort_open_unknown(theunit,thefile)
        integer theunit
        character*(*) thefile
        open(UNIT=theunit,FILE=thefile,STATUS='unknown')
      end

      subroutine fort_close(theunit)
        integer theunit
        close(theunit)
      end

      subroutine sjf_test
        include 'egs5/include/egs5_h.f'
        include 'egs5/include/egs5_media.f'
        character*24 medarr(2)
        medarr(1)='AIR AT NTP              '
        do i=1,24
           media(i,1)=medarr(1)(i:i) 
        end do
      end
