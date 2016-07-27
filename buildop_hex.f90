MODULE grid

IMPLICIT NONE

! COUNTING

! Number of grids in multigrid hierarchy
INTEGER :: ngrids

! nface        Number of faces on each grid
! nedge        Number of edges on each grid
! nvert        Number of vertices on each edge 
INTEGER, ALLOCATABLE :: nface(:), nedge(:), nvert(:)
INTEGER :: nfacex, nedgex, nvertx

! neoff        Number of edges and vertices of each face on each grid
! neofv        Number of edges of each vertex on each grid
INTEGER, ALLOCATABLE :: neoff(:,:), neofv(:,:)
INTEGER :: nefmx, nevmx


! CONNECTIVITY

! fnxtf        Faces next to each face on each grid
! eoff         Edges of each face on each grid
! voff         Vertices of each face on each grid
! fnxte        Faces either side of each edge on each grid
! vofe         Vertices at the ends of each edge on each grid
! fofv         Faces around each vertex on each grid
! eofv         Edges incident on each vertex on each grid
INTEGER, ALLOCATABLE :: fnxtf(:,:,:), eoff(:,:,:), voff(:,:,:), &
                        fnxte(:,:,:), vofe(:,:,:), &
                        fofv(:,:,:), eofv(:,:,:)

! Conventions
!
! 1. fnxtf and eoff are ordered anticlockwise and such that
! the i'th edge lies between the face in question and its i'th
! neighbour.
!
! 2. voff are ordered anticlockwise such that the k'th vertex
! is the common vertex of the k'th and (k+1)'th edge.
!
! 3. The positive normal direction n points from
! fnxte(:,1,:) -> fnxte(:,2,:)
! and the positive tangential direction t points from
! vofe(:,1,:) -> vofe(:,2,:)
! such that t = k x n (where k is vertical).


! eoffin(f,j,:)   Indicates whether the normal at the j'th edge is
!                 inward or outward relative to face f.
! eofvin(v,j,:)   Indicates whether the tangent at the j'th edge is
!                 inward or outward relative to vertex v.
INTEGER, ALLOCATABLE :: eoffin(:,:,:), eofvin(:,:,:)


! COORDINATES AND GEOMETRICAL INFORMATION

! flong        Longitude of faces on each grid
! flat         Latitude of faces on each grid
! vlong        Longitude of vertices on each grid
! vlat         Latitude of vertices on each grid
! farea        Area of faces on each grid
! ldist        Primal edge length, i.e. distance between neighbouring face centres
! ddist        Dual edge length, i.e. distance between neighbouring vertices
REAL*8, ALLOCATABLE :: flong(:,:), flat(:,:), &
                       vlong(:,:), vlat(:,:), &
                       farea(:,:), ldist(:,:), ddist(:,:)

! Conventions
!
! Latitude and longitude in radians.
! Area and lengths are for the unit sphere.


! HODGE STAR AND RELATED OPERATORS

! nisten        Number of faces in stencil for I operator
! isten         Stencil for I operator
! istar         Coefficients for I operator
! njsten        Number of vertices in stencil for J operator
! jsten         Stencil for J operator
! jstar         Coefficients for J operator
! nhsten        Number of edges in stencil for H operator
! hsten         Stencil for H operator
! hstar         Coefficients for H operator
! nrsten        Number of vertices in stencil for R operator (= neoff)
! rsten         Stencil for R operator (= voff)
! rcoeff        Coefficients for R operator
INTEGER, ALLOCATABLE :: nisten(:,:), njsten(:,:), nhsten(:,:), nrsten(:,:)
INTEGER, ALLOCATABLE :: isten(:,:,:), jsten(:,:,:), &
                        hsten(:,:,:), rsten(:,:,:)
REAL*8, ALLOCATABLE :: istar(:,:,:), jstar(:,:,:), &
                       hstar(:,:,:), rcoeff(:,:,:)
INTEGER :: nismx, njsmx, nhsmx, nrsmx


! INJECTION AND PROLONGATION OPERATORS FOR MULTIGRID

! ninj           Number of faces in stencil for injection operator
! injsten        Stencil for injection operator
! injwgt         Weights for injection operator
INTEGER, ALLOCATABLE :: ninj(:,:), injsten(:,:,:)
REAL*8, ALLOCATABLE :: injwgt(:,:,:)
INTEGER :: ninjmx


END MODULE grid

! ================================================================

MODULE channels

! Tidy list of all I/O channels in one place to avoid accidental
! overuse of any channel number


INTEGER, PARAMETER :: changin = 25          ! Input grid information
INTEGER, PARAMETER :: changout = 26         ! Output grid and operator information


END MODULE channels

! ===============================================================

PROGRAM buildop

! To take basic grid information from grid generator program
! and construct mimetic operators and other operators needed
! for a shallow water model.
!
! This version: Hexagonal-icoshedral C grid
!
! John Thuburn December 2011
!
! ---------------------------------------------------------------

! Read input grid data
CALL readgrid
PRINT *,'Input data read'

! Build I operator
CALL buildI
PRINT *,'I operator built'

! Build J operator
CALL buildJ
PRINT *,'J operator built'

! Build H operator
CALL buildH
PRINT *,'H operator built'

! Build R operator
CALL buildR
PRINT *,'R operator built'

! Build Injection operator
CALL buildinj
PRINT *,'Injection operator built'

! Write everything out
CALL writegrid
PRINT *,'Data written'

! ---------------------------------------------------------------

END PROGRAM buildop

! ===============================================================

SUBROUTINE buildI

! To build the mimetic I operator. This operator is diagonal
! with coefficient given by the inverse of the face area

USE grid

IMPLICIT NONE
INTEGER :: if0, igrid

! ---------------------------------------------------------------

! Allocate array for size of operator stencil
ALLOCATE(nisten(nfacex,ngrids))

! Stencil is a single cell
nisten = 1

! Find maximum values in order to allocate subsequent arrays
nismx = MAXVAL(nisten)

! Allocate arrays for operator stencils and coefficients
ALLOCATE(isten(nfacex,nismx,ngrids))
ALLOCATE(istar(nfacex,nismx,ngrids))

! Initialize to zero
isten = 0
istar = 0.0d0

! And define stencil and coefficient
DO igrid = 1, ngrids
  DO if0 = 1, nface(igrid)
    isten(if0,1,igrid) = if0
    istar(if0,1,igrid) = 1.0d0/farea(if0,igrid)
  ENDDO
ENDDO

! ---------------------------------------------------------------

END SUBROUTINE buildI

! ===============================================================

SUBROUTINE buildJ

! To build the mimetic J operator. This operator is diagonal
! with coefficient given by the inverse of the dual cell area

USE grid

IMPLICIT NONE
INTEGER :: iv0, igrid, if0
REAL*8 :: long, lat, x1, y1, z1, x2, y2, z2, x3, y3, z3, varea

! ---------------------------------------------------------------

! Allocate array for size of operator stencil
ALLOCATE(njsten(nvertx,ngrids))

! Stencil is a single vertex
njsten = 1

! Find maximum values in order to allocate subsequent arrays
njsmx = MAXVAL(njsten)

! Allocate arrays for operator stencils and coefficients
ALLOCATE(jsten(nvertx,njsmx,ngrids))
ALLOCATE(jstar(nvertx,njsmx,ngrids))

! Initialize to zero
jsten = 0
jstar = 0.0d0

! And define stencil and coefficient
DO igrid = 1, ngrids
  DO iv0 = 1, nvert(igrid)
    jsten(iv0,1,igrid) = iv0
    ! Dual cell is triangular so use formula for area
    ! of spherical triangle
    if0 = fofv(iv0,1,igrid)
    long = flong(if0,igrid)
    lat = flat(if0,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    if0 = fofv(iv0,2,igrid)
    long = flong(if0,igrid)
    lat = flat(if0,igrid)
    CALL ll2xyz(long,lat,x2,y2,z2)
    if0 = fofv(iv0,3,igrid)
    long = flong(if0,igrid)
    lat = flat(if0,igrid)
    CALL ll2xyz(long,lat,x3,y3,z3)
    CALL starea2(x1,y1,z1,x2,y2,z2,x3,y3,z3,varea)
    jstar(iv0,1,igrid) = 1.0d0/varea
  ENDDO
ENDDO

! ---------------------------------------------------------------

END SUBROUTINE buildJ

! ===============================================================

SUBROUTINE buildH

! To build the mimetic H operator. This operator is diagonal
! with coefficient given by the ratio ldist/ddist

USE grid

IMPLICIT NONE
INTEGER :: ie0, igrid

! ---------------------------------------------------------------

! Allocate array for size of operator stencil
ALLOCATE(nhsten(nedgex,ngrids))

! Stencil is a single edge
nhsten = 1

! Find maximum values in order to allocate subsequent arrays
nhsmx = MAXVAL(nhsten)

! Allocate arrays for operator stencils and coefficients
ALLOCATE(hsten(nedgex,nhsmx,ngrids))
ALLOCATE(hstar(nedgex,nhsmx,ngrids))

! Initialize to zero
hsten = 0
hstar = 0.0d0

! And define stencil and coefficient
DO igrid = 1, ngrids
  DO ie0 = 1, nedge(igrid)
    hsten(ie0,1,igrid) = ie0
    hstar(ie0,1,igrid) = ldist(ie0,igrid)/ddist(ie0,igrid)
  ENDDO
ENDDO

! ---------------------------------------------------------------

END SUBROUTINE buildH

! ===============================================================

SUBROUTINE buildR

! To build the R operator. The stencil for face if0 is the set
! of vertices around face if0. The coefficient is the ratio of
! the kite-shaped corner area to the full face area.

USE grid

IMPLICIT NONE
INTEGER :: if0, if2, igrid, ix1, ix2, ie1, ie2, iv0
REAL*8 :: long, lat, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          area1, area2, rmag

! ---------------------------------------------------------------

! Allocate array for size of operator stencil
ALLOCATE(nrsten(nfacex,ngrids))

! Stencil is the set of vertices of each face
nrsten = neoff

! Find maximum values in order to allocate subsequent arrays
nrsmx = MAXVAL(nrsten)

! Allocate arrays for operator stencils and coefficients
ALLOCATE(rsten(nfacex,nrsmx,ngrids))
ALLOCATE(rcoeff(nfacex,nrsmx,ngrids))

! Initialize to zero
rsten = 0
rcoeff = 0.0d0

! And define stencil and coefficients
DO igrid = 1, ngrids
  DO if0 = 1, nface(igrid)
    ! Face centre
    long = flong(if0,igrid)
    lat = flat(if0,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    ! Loop over vertices of face
    DO ix1 = 1, nrsten(if0,igrid)
      ix2 = ix1 + 1
      IF (ix2 > nrsten(if0,igrid)) ix2 = 1
      ! ix1'th vertex is the common vertex of the
      ! ix1'th and ix2'th edges
      iv0 = voff(if0,ix1,igrid)
      ie1 = eoff(if0,ix1,igrid)
      ie2 = eoff(if0,ix2,igrid)
      ! Coordinates of vertex
      long = vlong(iv0,igrid)
      lat = vlat(iv0,igrid)
      CALL ll2xyz(long,lat,x1,y1,z1)
      ! Each edge contributes a triangle to
      ! the kite area.
      ! Find crossing point of each edge and its dual
      ! (using Voronoi property). 
      if2 = fnxtf(if0,ix1,igrid)
      long = flong(if2,igrid)
      lat = flat(if2,igrid)
      CALL ll2xyz(long,lat,x2,y2,z2)
      x2 = 0.5d0*(x0 + x2)
      y2 = 0.5d0*(y0 + y2)
      z2 = 0.5d0*(z0 + z2)
      rmag = 1.0d0/SQRT(x2*x2 + y2*y2 + z2*z2)
      x2 = x2*rmag
      y2 = y2*rmag
      z2 = z2*rmag
      CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area1)
      if2 = fnxtf(if0,ix2,igrid)
      long = flong(if2,igrid)
      lat = flat(if2,igrid)
      CALL ll2xyz(long,lat,x2,y2,z2)
      x2 = 0.5d0*(x0 + x2)
      y2 = 0.5d0*(y0 + y2)
      z2 = 0.5d0*(z0 + z2)
      rmag = 1.0d0/SQRT(x2*x2 + y2*y2 + z2*z2)
      x2 = x2*rmag
      y2 = y2*rmag
      z2 = z2*rmag
      CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area2)
      rsten(if0,ix1,igrid) = iv0
      rcoeff(if0,ix1,igrid) = (area1 + area2)/farea(if0,igrid)
    ENDDO
  ENDDO

ENDDO

! ---------------------------------------------------------------

END SUBROUTINE buildR

! ===============================================================

SUBROUTINE buildinj

! To build injection operator for multigrid

USE grid

IMPLICIT NONE
INTEGER :: if0, if1, igrid, igridp, nf, ix, ixp

! ---------------------------------------------------------------

! Allocate array for size of operator stencil
ALLOCATE(ninj(nfacex,ngrids-1))

! Stencil is 6 or 7 cells on standard buckyball grid
ninjmx = 7

! Allocate arrays for operator stencils and coefficients
ALLOCATE(injsten(nfacex,ninjmx,ngrids-1))
ALLOCATE(injwgt(nfacex,ninjmx,ngrids-1))

! Initialize to zero
injsten = 0
injwgt = 0.0d0

! And define stencil and coefficients
DO igrid = 1, ngrids-1
  igridp = igrid + 1
  DO if0 = 1, nface(igrid)
    ! Face if0 exists on grid igrid and grid igrid+1 and is
    ! the centre of the stencil
    injsten(if0,1,igrid) = if0
    injwgt(if0,1,igrid) = 1.0d0
    ! The neighbours of if0 on grid igrid+1 are the other members
    ! of the stencil
    nf = neoff(if0,igridp)
    ninj(if0,igrid) = 1 + nf
    DO ix = 1, nf
      ixp = ix + 1
      if1 = fnxtf(if0,ix,igridp)
      injsten(if0,ixp,igrid) = if1
      injwgt(if0,ixp,igrid) = 0.5d0
    ENDDO
  ENDDO
ENDDO

! ---------------------------------------------------------------

END SUBROUTINE buildinj

! ===============================================================
!
      SUBROUTINE LL2XYZ(LONG,LAT,X,Y,Z)
!
!     To convert longitude and latitude to cartesian coordinates
!     on the unit sphere
!
      IMPLICIT NONE
!
      REAL*8 LONG,LAT,X,Y,Z,CLN,SLN,CLT,SLT
!
!     ------------------------------------------------------------------
!
      SLN=SIN(LONG)
      CLN=COS(LONG)
      SLT=SIN(LAT)
      CLT=COS(LAT)
!
      X=CLN*CLT
      Y=SLN*CLT
      Z=SLT
!
!     ------------------------------------------------------------------
!
      RETURN
      END
!
!     ==================================================================
!
      SUBROUTINE STAREA2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AREA)
!
!     Calculate the area of the spherical triangle whose corners
!     have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
!     The formula below is more robust to roundoff error than the
!     better known sum of angle - PI formula
!
      IMPLICIT NONE
!
      REAL*8 X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2, &
             D0,D1,D2,S,T0,T1,T2,T3,AREA
!
!
!     Distances between pairs of points
      CALL SPDIST(X0,Y0,Z0,X1,Y1,Z1,D2)
      CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,D0)
      CALL SPDIST(X2,Y2,Z2,X0,Y0,Z0,D1)
!
!     Half perimeter
      S=0.5D0*(D0+D1+D2)
!
!     Tangents
      T0 = TAN(0.5D0*(S-D0))
      T1 = TAN(0.5D0*(S-D1))
      T2 = TAN(0.5D0*(S-D2))
      T3 = TAN(0.5D0*S)
!
!     Area
      AREA = 4.0D0*ATAN(SQRT(T0*T1*T2*T3))
!
      RETURN
      END
!
!     ===================================================================
!
      SUBROUTINE SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
!
!     Calculate the spherical distance S between two points with Cartesian
!     coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

      IMPLICIT NONE

      REAL*8 X1, Y1, Z1, X2, Y2, Z2, S, DX, DY, DZ, AD


      DX = X2 - X1
      DY = Y2 - Y1
      DZ = Z2 - Z1
      AD = SQRT(DX*DX + DY*DY + DZ*DZ)
      S = 2.0D0*ASIN(0.5D0*AD)


      RETURN
      END
!
!     ===================================================================

SUBROUTINE readgrid

! To allocate array space for the grid information in module grid
! and to read the information from file

USE grid
USE channels

IMPLICIT NONE

INTEGER :: if0, ie0, iv0, igrid, ix
CHARACTER*26 :: ygridfile

! ----------------------------------------------------------------

! Ask user for resolution - this will determine filename
WRITE(6,*) 'Resolution of input file: number of faces'
READ(5,*) nfacex
WRITE(ygridfile,'(''gridmap_hex_'',I10.10,''.dat'')') nfacex

! Open file for reading
OPEN(changin,FILE=ygridfile,FORM='UNFORMATTED')

! First read ngrids
READ(changin) ngrids

! Allocate nface, nedge, nvert
ALLOCATE(nface(ngrids), nedge(ngrids), nvert(ngrids))

! Read numbers of faces, edges and vertices on each grid
READ(changin) nface
READ(changin) nedge
READ(changin) nvert

! Find maximum values in order to allocated subsequent arrays
nfacex = MAXVAL(nface)
nedgex = MAXVAL(nedge)
nvertx = MAXVAL(nvert)

! Allocate neoff, neofv
ALLOCATE(neoff(nfacex,ngrids), neofv(nvertx,ngrids))
neoff = 0
neofv = 0

! Read the numbers of edges of each face and vertex on each grid
READ(changin) ((neoff(if0,igrid),           &
                    if0 = 1, nface(igrid)), &
                    igrid = 1, ngrids)
READ(changin) ((neofv(iv0,igrid),           &
                    iv0 = 1, nvert(igrid)), &
                    igrid = 1, ngrids)

! Find maximum values in order to allocate subsequent arrays
nefmx = MAXVAL(neoff)
nevmx = MAXVAL(neofv)

! Allocate connectivity arrays arrays
ALLOCATE(fnxtf(nfacex,nefmx,ngrids), eoff(nfacex,nefmx,ngrids), &
         voff(nfacex,nefmx,ngrids),  fnxte(nedgex,2,ngrids),    &
         vofe(nedgex,2,ngrids),      fofv(nvertx,nevmx,ngrids), &
         eofv(nvertx,nevmx,ngrids))

! Read the connectivity arrays
READ(changin) (((fnxtf(if0,ix,igrid),           &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changin) (((eoff(if0,ix,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changin) (((voff(if0,ix,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changin) (((fnxte(ie0,ix,igrid),           &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
READ(changin) (((vofe(ie0,ix,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
READ(changin) (((fofv(iv0,ix,igrid),            &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)
READ(changin) (((eofv(iv0,ix,igrid),            &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)

! Allocate the geometrical information arrays
ALLOCATE(flong(nfacex,ngrids), flat(nfacex,ngrids), &
         vlong(nvertx,ngrids), vlat(nvertx,ngrids), &
         farea(nfacex,ngrids), &
         ldist(nedgex,ngrids), ddist(nedgex,ngrids))

! Read the geometrical information arrays
READ(changin) ((flong(if0,igrid),               &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((flat(if0,igrid),                &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((vlong(iv0,igrid),               &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((vlat(iv0,igrid),                &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((farea(if0,igrid),               &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((ldist(ie0,igrid),               &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((ddist(ie0,igrid),               &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)

!     ---------------------------------------------------------------

END SUBROUTINE readgrid

!     ===============================================================

SUBROUTINE writegrid

! To write the grid and operator information to file

USE grid
USE channels

IMPLICIT NONE

INTEGER :: if0, ie0, iv0, igrid, ix
CHARACTER*30 :: ygridfile

! ----------------------------------------------------------------

! Construct filename
WRITE(ygridfile,'(''gridopermap_hex_'',I10.10,''.dat'')') nfacex

! Open file for writing
OPEN(changout,FILE=ygridfile,FORM='UNFORMATTED')


! First write ngrids
WRITE(changout) ngrids

! Write numbers of faces, edges and vertices on each grid
WRITE(changout) nface
WRITE(changout) nedge
WRITE(changout) nvert

! Write the numbers of edges of each face and vertex on each grid
WRITE(changout) ((neoff(if0,igrid),         &
                    if0 = 1, nface(igrid)), &
                    igrid = 1, ngrids)
WRITE(changout) ((neofv(iv0,igrid),         &
                    iv0 = 1, nvert(igrid)), &
                    igrid = 1, ngrids)


! Write the connectivity arrays
WRITE(changout) (((fnxtf(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((eoff(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((voff(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((fnxte(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
WRITE(changout) (((vofe(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
WRITE(changout) (((fofv(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((eofv(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)


! Write the geometrical information arrays
WRITE(changout) ((flong(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((flat(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((vlong(iv0,igrid),             &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((vlat(iv0,igrid),              &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((farea(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ldist(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ddist(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)


! Write the sizes of the operator stencils on each grid
WRITE(changout) ((nisten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((njsten(iv0,igrid),            &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nhsten(ie0,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nrsten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)


! Write the operator stencils and coefficients
WRITE(changout) (((isten(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nismx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((jsten(iv0,ix,igrid),         &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((hsten(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rsten(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((istar(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nismx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((jstar(iv0,ix,igrid),         &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((hstar(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rcoeff(if0,ix,igrid),        &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)


! Write the size of the injection stencil
WRITE(changout) ((ninj(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids-1)

! Write the injection stencil and coefficients
WRITE(changout) (((injsten(if0,ix,igrid),       &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ninjmx),           &
                     igrid = 1, ngrids-1)
WRITE(changout) (((injwgt(if0,ix,igrid),        &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ninjmx),           &
                     igrid = 1, ngrids-1)


!     ---------------------------------------------------------------

END SUBROUTINE writegrid

!     ===============================================================
