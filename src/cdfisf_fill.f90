PROGRAM cdfisf_fill
  !!======================================================================
  !!                     ***  PROGRAM  cdfisf_fill  ***
  !!=====================================================================
  !!  ** Purpose : Build a file containing one value for each closed pools
  !!               seeded by a list of points.
  !!
  !!  ** Method  : flood filling algorithm
  !!               
  !! History : 3.0  : 04/2014  : Pierre Mathiot 
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  USE modutils
  USE cdftools
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class ice_shelf_processing
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jisf, ji, jj       ! dummy loop integer 
  INTEGER(KIND=4)                               :: ierr, ipos         ! working integer
  INTEGER(KIND=4)                               :: idep, idep_max     ! possible depth index, maximum
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt, nisf     ! size of the domain
  INTEGER(KIND=4)                               :: iunit=10           ! file unit for txt input file
  INTEGER(KIND=4)                               :: iunitu=11          ! file unit for txt output file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4)                               :: iiseed, ijseed
  INTEGER(KIND=4)                               :: ifill
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: itab, itabcnt

  REAL(KIND=4)                                  :: rlon, rlat         ! longitude and latitude of one point in ISF
  REAL(KIND=4)                                  :: rdraftmin, rdraftmax
  REAL(KIND=4)                                  :: rdum, rdum0, rdum1
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab, ssmask, e1, e2, area

  CHARACTER(LEN=256)                            :: cf_in              ! input file name
  CHARACTER(LEN=256)                            :: cf_isflist         ! input file name (txt)
  CHARACTER(LEN=256)                            :: cf_isflistup       ! output file name (update of input, with draftmin/max
  CHARACTER(LEN=256)                            :: cf_out='mskisf.nc'   ! output file for average
  CHARACTER(LEN=256)                            :: cf_boundary        ! boundary file
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cv_in              ! depth dimension name
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument
  CHARACTER(LEN=256)                            :: cisf               ! dummy string argument
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep            ! array of possible depth name (or 3rd dimension)

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values
  LOGICAL                                       :: lnc4 = .FALSE.     ! flag for netcdf4 chunk and deflation
  LOGICAL                                       :: lperio = .FALSE.   ! flag for input file periodicity
  LOGICAL                                       :: ldiag = .FALSE.
  LOGICAL                                       :: lboundf = .FALSE.

  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfisf_fill  -f ISF-file -v ISF-var -l ISF-list [-nc4 ] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE : Builds nc file with a single value for each pool around a list'
     PRINT *,'               of given point. A warning is given when neighbouring ice-shelves'
     PRINT *,'               cannot be discriminated (no gap in between). In this case, hand'
     PRINT *,'               edit on the ISF-file is required.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS : '
     PRINT *,'       [-l ISF-list] : text file containing at least the following information: '
     PRINT *,'                 1  NAME    LON  LAT I  J DUM DUM DUM DUM ldiag'
     PRINT *,'                 ...             '
     PRINT *,'                 i  NAMEi   LON  LAT I  J DUM DUM DUM DUM ldiag'
     PRINT *,'                 ...             '
     PRINT *,'                 EOF             '
     PRINT *,'                 No NAME  X    Y   I  J DUM DUM DUM DUM ldiag'
     PRINT *,'       [-bf txtfile] : txt file describing the section used in -fill'
     PRINT *,'                        Extra boundary could be set up in boundary.txt.'
     PRINT *,'                        Format of the file is on each line : '
     PRINT *,'                        NAME /n iimin iimax jjmin jjmax linc.'
     PRINT *,'                        Section is exclude from the selection if linc=F .'
     PRINT *,'       [ -ew ]       : input file are periodic along e/w direction'
     PRINT *,'      '
     PRINT *,'     OPTIONS : '
     PRINT *,'          -nc4 : use NetCDF4 chunking and deflation for the output'
     PRINT *,'          -o OUT-file : specify the name of the output file instead of ',TRIM(cf_out)
     PRINT *,'                 This file will be one of the input file for cdfmkforcingisf '
     PRINT *,'                 as the ISF-fill_file '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'              netcdf file : fill.nc '
     PRINT *,'              variable : sofillvar contains for all points in ice shelf NAME '
     PRINT *,'                         the value -i (negative value)'
     PRINT *,'              text file : <ISF-list>_zmin_zmax.txt '
     PRINT *,'                        this output file is similar to <ISF-list> but updated'
     PRINT *,'                        with the minimum and maximul value of ice-draft for '
     PRINT *,'                        each shelf.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'           cdfisf_forcing,  cdfisf_rnf , cdfisf_poolchk'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 
     SELECT CASE ( cldum)
     CASE ( '-ew') ; lperio = .true.
     CASE ( '-bf') ; CALL getarg(ijarg, cf_boundary) ; lboundf=.TRUE. ; ijarg = ijarg + 1
     CASE ( '-l' ) ; CALL getarg(ijarg, cf_isflist ) ;                ijarg = ijarg + 1
     CASE ( '-o' ) ; CALL getarg(ijarg, cf_out     ) ;                ijarg = ijarg + 1
     CASE ('-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT  ; PRINT *,' ERROR : ', TRIM(cldum) ,' : unknown option.'; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile (cn_fmsk) .OR. chkfile (cn_fzgr) .OR. chkfile (cf_isflist)  ) STOP 99 ! missing file

  ipos = INDEX(cf_isflist,'.')
  cldum=cf_isflist(ipos+1:)
  cf_isflistup=cf_isflist(1:ipos-1)//'_zmin_zmax.'//TRIM(cldum)

  npiglo = getdim (cn_fzgr, cn_x)
  npjglo = getdim (cn_fzgr, cn_y)
  npk    = getdim (cn_fzgr, cn_z) 

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'NPK    = ', npk

  ALLOCATE(dtab(npiglo, npjglo), ssmask(npiglo, npjglo), e1(npiglo, npjglo), e2(npiglo, npjglo), area(npiglo, npjglo))
  ALLOCATE(itab(npiglo, npjglo), itabcnt(npiglo, npjglo))

  ALLOCATE (stypvar(2))
  ALLOCATE (ipk(2),id_varout(2))

  CALL CreateOutput 

  ! initialize variable
  dtab(:,:) = 0.d0 
  
  ! get area
  e1 = getvar(cn_fhgr,cn_ve1t, 1 ,npiglo, npjglo )
  e2 = getvar(cn_fhgr,cn_ve1t, 1 ,npiglo, npjglo )
  area(:,:) = e1(:,:) * e2(:,:)

  ! read ice shelf draft data
  itab=1
  itabcnt = 0

  ! read ice shelf draft or top level data
  dtab = getvar(cn_fzgr, 'isfdraft', 1 ,npiglo, npjglo )
  ssmask = getvar(cn_fmsk, cn_tmask, 1 ,npiglo, npjglo )
  
  WHERE ( dtab <=0 ) itab=0
  PRINT *, 'Maximum of ISF-draft : ', MAXVAL(dtab),' m'

  ! open isf-list file
  OPEN(unit=iunit,  file=cf_isflist  , form='formatted', status='old')
  OPEN(unit=iunitu, file=cf_isflistup, form='formatted'              )
  ! get total number of isf
  nisf = 0
  cisf='XXX'
  DO WHILE ( TRIM(cisf) /= 'EOF')
     READ(iunit,*) cisf
     nisf=nisf+1
  END DO
  REWIND(iunit)

  nisf = nisf - 1
  PRINT *, '   Number of ISF found in file list : ', nisf

  ! loop over each ice shelf
  DO jisf=1,nisf
     ! get iiseed, ijseed, ice shelf number ifill
     READ(iunit,*) ifill, cisf, rlon, rlat, iiseed, ijseed, rdum, rdum, rdum0, rdum1, ldiag
     IF (dtab(iiseed, ijseed) < 0 ) THEN
        PRINT *,'  ==> WARNING: Likely a problem with ',TRIM(cldum)
        PRINT *,'               check separation with neighbours'
     ENDIF

     IF (ifill < 99) THEN
        ! add section boundary
        IF (lboundf) CALL update_mask(itab,-ifill)
        IF (itab(iiseed, ijseed) < 0 ) THEN
           PRINT *,'  ==> WARNING: Likely a problem with ',TRIM(cisf)
           PRINT *,'               check separation with neighbours'
        ENDIF
        IF (dtab(iiseed, ijseed) < 2) THEN
           PRINT *, '  ==> ERROR: Trouble with seed for isf : ',TRIM(cisf),' id : ',ifill
           STOP
        END IF

        ! fill ice shelf cavities
        CALL FillPool2D(iiseed, ijseed, itab, -ifill, lperio, ldiag)

        ! range of depth for each ice shelf
        rdraftmax=MAXVAL(dtab, (itab == -ifill) )
        rdraftmin=MINVAL(dtab, (itab == -ifill) )

        ! find ice shelf front cell
        DO ji=2,npiglo-1
           DO jj = 2, npjglo-1
              IF ( ssmask(ji,jj) == 1 .AND.                        &
                   & MINVAL(itab(ji-1:ji+1 , jj-1:jj+1)) == -ifill ) THEN
                 itabcnt(ji,jj)  = -ifill
              END IF
           END DO
        END DO
        IF ( lperio ) THEN
           itabcnt(1     ,:) = itabcnt(npiglo-1,:)
           itabcnt(npiglo,:) = itabcnt(2       ,:)
        END IF
 
        PRINT *,'Iceshelf   : ', TRIM(cisf)
        PRINT *,'  index    : ', ifill
        PRINT *,'  seed val.: ', INT(dtab(iiseed, ijseed ) )
        PRINT *,'  depmin [m]   : ', rdraftmin
        PRINT *,'  depmax [m]   : ', rdraftmax
        PRINT *,'  area   [km2] : ', SUM(area, (itab == -ifill) ) * 1e-6 ! kmsq
        PRINT *,'   '
     END IF
     WRITE(iunitu,'(i4,1x,a20,2f9.4,2i5,2f8.1,i8)') ifill,ADJUSTL(cisf),rlon, rlat, iiseed, ijseed,rdraftmin,rdraftmax,rdum0,rdum1,ldiag
  END DO
  WRITE(iunitu,'(a)') 'EOF  '
  WRITE(iunitu,'(a5,a20,2a9,2a5,2a8,a)' ) 'No ','NAME                           ',' X',' Y',' I ',' J ',' Zmin',' Zmax',' FWF'

  CLOSE(iunitu)
  CLOSE(iunit)

  ierr = putvar(ncout, id_varout(1), itab, 1, npiglo, npjglo)

 ierr = putvar(ncout, id_varout(2), itabcnt, 1, npiglo, npjglo)

  ! close file
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create the output file. This is done outside the main
    !!               in order to increase readability of the code. 
    !!
    !! ** Method  :  Use global variables, defined in main 
    !!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(1) :: dl_tim
    !!----------------------------------------------------------------------
    
    ! define new variables for output
    ipk(1) = 1  !  2D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'mask_isf'
    stypvar(1)%cunits            = 'N/A'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = -1000.
    stypvar(1)%valid_max         =  1000.
    stypvar(1)%clong_name        = 'Mask of each ice shelf (id different for each)'
    stypvar(1)%cshort_name       = 'mask_isf'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'
    stypvar(1)%cprecision        = 'i2'

    ! define new variables for output
    ipk(2) = 1  !  2D
    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = 'mask_isf_front'
    stypvar(2)%cunits            = 'N/A'
    stypvar(2)%rmissing_value    = 0.
    stypvar(2)%valid_min         = -1000.
    stypvar(2)%valid_max         =  1000.
    stypvar(2)%clong_name        = 'Mask of each ice shelf front (id different for each)'
    stypvar(2)%cshort_name       = 'mask_isf_front'
    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = 'TYX'
    stypvar(2)%cprecision        = 'i2'

    ! create output file taking the sizes in cf_in
    ncout  = create      (cf_out,  cn_fzgr,    npiglo, npjglo, 1, ld_nc4=lnc4)
    ierr   = createvar   (ncout ,  stypvar,  2,  ipk,    id_varout, ld_nc4=lnc4)
    ierr   = putheadervar(ncout,   cn_fzgr,    npiglo, npjglo, 1 )

    dl_tim(1)=0.d0
    ierr  = putvar1d(ncout, dl_tim, 1, 'T')

  END SUBROUTINE CreateOutput

  SUBROUTINE update_mask(imask, kfill)

    INTEGER(KIND=4), PARAMETER                   :: jseg=10000   ! dummy loop index
    INTEGER(KIND=4)                              :: ipos           ! working integer (position of ' ' in strings)
    INTEGER(KIND=4)                              :: ii, jk, js     ! working integer
    INTEGER(KIND=4)                              :: iunitb=12
    INTEGER(KIND=4)                              :: iimin, iimax, iipts      ! limit in i
    INTEGER(KIND=4)                              :: ijmin, ijmax, ijpts      ! limit in i
    INTEGER(KIND=4)                              :: nsec

    INTEGER(KIND=4), INTENT(in) :: kfill                                          ! filling value
    INTEGER(KIND=2), DIMENSION(:,:), INTENT(inout)  :: imask          ! mask
    REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE    :: rxx, ryy       ! working variables

    CHARACTER(LEN=256)                           :: cline          ! dummy char variable
    CHARACTER(LEN=256)                           :: csection       ! section names

    LOGICAL :: lsection

    ALLOCATE(rxx(npiglo+npjglo), ryy(npiglo+npjglo)) ! dimension specified in broken_line
    ! optimal dimension could be ABS(imax-imin +1) + ABS(jmax-jmin+1) - 1

    IF (lboundf) THEN
       PRINT *,''
       PRINT *,'Boundary file: ',TRIM(cf_boundary),' is used to close the basin'
       PRINT *,''

       IF ( chkfile(cf_boundary) ) STOP 99 ! missing file

       OPEN(iunitb, FILE=cf_boundary)
       lsection = .TRUE.
       DO WHILE (lsection)
          rxx(:)=1; ryy(:)=1

          ! read section name
          READ(iunitb,'(a)') csection
          IF (TRIM(csection) == 'EOF' ) THEN
             lsection = .FALSE.
          ELSEIF (TRIM(csection) == TRIM(cisf) ) THEN
             ! read section coordinates
             READ(iunitb,*) nsec
             DO js = 1,nsec
                READ(iunitb,*) iimin, iimax, ijmin, ijmax

                ! get index of cell included into the section
                CALL broken_line(iimin, iimax, ijmin, ijmax, rxx, ryy, npt, npiglo, npjglo)
 
                ! mask boundary and keep location in rmskline
                DO jk=1,npt
                   imask(rxx(jk),ryy(jk))=kfill*imask(rxx(jk),ryy(jk))
                END DO
             END DO
          ENDIF
       END DO
       CLOSE(iunitb)
    ELSE
       PRINT *,''
       PRINT *, 'NO BOUNDARIES ARE ADDED TO THE INPUT FILE TO CLOSE THE BASIN'
       PRINT *,''
    END IF

  END SUBROUTINE update_mask

END PROGRAM cdfisf_fill
