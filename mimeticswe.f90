MODULE grid

IMPLICIT NONE

! COUNTING

! Number of grids in multigrid hierarchy
INTEGER :: ngrids

! nface        Number of faces on each grid
! nedge        Number of edges on each gridcp
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
! 2. The positive normal direction n points from
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
! fareamin     Minimum face areaon each grid
REAL*8, ALLOCATABLE :: flong(:,:), flat(:,:), &
                       vlong(:,:), vlat(:,:), &
                       farea(:,:), ldist(:,:), ddist(:,:), &
                       fareamin(:)

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
! kecoeff       Coefficients for computing KE
INTEGER, ALLOCATABLE :: nisten(:,:), njsten(:,:), nhsten(:,:), nrsten(:,:)
INTEGER, ALLOCATABLE :: isten(:,:,:), jsten(:,:,:), &
                        hsten(:,:,:), rsten(:,:,:)
REAL*8, ALLOCATABLE :: istar(:,:,:), jstar(:,:,:), &
                       hstar(:,:,:), rcoeff(:,:,:), kecoeff(:,:,:)
INTEGER :: nismx, njsmx, nhsmx, nrsmx


! RESTRICTION AND PROLONGATION OPERATORS FOR MULTIGRID

! ninj           Number of faces in stencil for restriction operator
! injsten        Stencil for restriction operator
! injwgt         Weights for restriction operator
INTEGER, ALLOCATABLE :: ninj(:,:), injsten(:,:,:)
REAL*8, ALLOCATABLE :: injwgt(:,:,:)
INTEGER :: ninjmx


END MODULE grid

! ================================================================

MODULE runtype

! Default values are set here; they may be overridden
! via the namelist

IMPLICIT NONE

! File containing grid information
!CHARACTER*31 :: ygridfile = 'gridopermap_hex_0000010242.dat '
CHARACTER*31 :: ygridfile = 'gridopermap_cube_0000013824.dat'

! File to output grid coordinates (for use in generating reference
! solutions).
!CHARACTER*30 :: ygridcoords = 'gridcoords_hex_0000010242.dat '
CHARACTER*30 :: ygridcoords = 'gridcoords_cube_0000013824.dat'

! Run identifier for naming output files (CURRENTLY NOT USED)
CHARACTER*6 :: runid = '000001'

! True for a restart run (CURRENTLY NOT USED)
LOGICAL :: lrestart

! Name of restart file (CURRENTLY NOT USED)
CHARACTER*30 :: yresfile = 'run000001_restart_00000000.dat'


END MODULE runtype

! ================================================================

MODULE constants

! Various physical and geometrical constants

IMPLICIT NONE

! Pi
REAL*8, PARAMETER :: pi = 3.14159265358979323d0, piby2 = 0.5d0*pi

! Earth's radius
REAL*8 :: rearth = 6371220.0d0

! Gravitational acceletration
REAL*8 :: gravity = 9.80616d0

! Rotation rate of the earth
REAL*8 :: rotatn = 7.29212d-5

! Dual cell integrals of planetary vorticity
REAL*8, ALLOCATABLE :: planvort2(:)


END MODULE constants

! ================================================================

MODULE state

! Variables describing the model state.
! Note on naming convention: an appended 2 indicates a variable
! integrated over a cell or dual cell area (a discrete 2-form);
! an appended 1 indicates a variable integrated along a
! primal or dual edge (a discrete 1-form). Names with neither
! indicate point values or length or area averages.

USE grid
IMPLICIT NONE

! OROGRAPHY
! orog2             Primal cell area integral of orography
REAL*8, ALLOCATABLE :: orog2(:)

! PROGNOSTIC VARIABLES
! phi2              Primal cell area integral of geopotential
! v1                Circulation along dual edge
! c2                Area integral of tracer concentration
REAL*8, ALLOCATABLE :: phi2(:), v1(:), c2(:)


! VARIABLES USED IN CHECKING THEORETICAL PROPERTIES
! xphibar2          Prognostic dual cell integral of geopotential
! xzeta2            Prognostic pv-like tracer
REAL*8, ALLOCATABLE :: xphibar2(:), xzeta2(:)


END MODULE state

! ================================================================

MODULE work

! Variables used in various parts of the calculation, 
! particularly the advection scheme, which it is therefore
! useful to store

IMPLICIT NONE

! vbar1             Time averaged circulation along dual edge
! ubar1             Time averaged flux across primal edge
! vperpbar1         Time averaged flux across dual edge
! uperpbar1         Time averaged circulation along primal edge
! mf1               Time integrated mass flux across primal edge
! mfperp1           Time integrated mass flux across dual edge
! divfac            Factor related to divergence used in modifying swept areas
REAL*8, ALLOCATABLE :: vbar1(:), ubar1(:), vperpbar1(:), uperpbar1(:), &
                       mf1(:), mfperp1(:), divfac(:)


REAL*8, ALLOCATABLE :: pva(:), u1_temp(:), &
          b0(:), b0_temp(:), gb1(:), div2(:), zeta2(:), &
          pv(:), qperp1(:), divmf(:), phirefe(:), &
          phi0(:), phi_inc(:), v_inc(:), &
          rphi2(:), rv1(:), rhs(:), &
          phi2_new(:), v1_new(:), &
          phibar2(:)


END MODULE work

! ================================================================

MODULE errdiag

! Variables used in computing error diagnostics

IMPLICIT NONE

! phi2_init           Initial geopotential
! v1_init             Initial v
REAL*8, ALLOCATABLE :: phi2_init(:), v1_init(:), &
                       phiexac(:), pvexac(:), &
                       phierr(:), pverr(:), v1err(:)

END MODULE errdiag

! ================================================================

MODULE helmcoeff

! Reference state values on the grid hierarchy needed for
! multigrid helmholtz solver

USE grid
IMPLICIT NONE

! phiref         Primal grid cell values of reference geopotential
!                on finest grid
! nusq           Cell edge values of reference geopotential
!                times alpha_pg*alpha_v*dt*dt
!                on all grids
! helmdiag       Diagonal coefficient of Helmholtz operator
!                on all grids
! underrel       Under-relaxation parameter on all grids
REAL*8, ALLOCATABLE :: phiref(:), nusq(:,:), helmdiag(:,:), underrel(:)


END MODULE helmcoeff

! ================================================================

MODULE advection

! Fields related to advection scheme

USE grid
IMPLICIT NONE

! degree         Degree of polynomial fit
! monotone       Use limiter if true (NOT YET IMPLEMENTED)
! nmonomial      Number of monomials to build polynomials of
!                given degree
! ngauss         Number of Gauss points needed in each direction
!                to integrate fit over swept area
! nstenadvf      Number of cells in stencil for given face
! nexadvf        Number of cells in stencil to be fitted exactly
! nadvfmx        Maximum number of cells in advection stencil
!                for given face
! stenadvf       Stencil for constructing polynomial fit for
!                given face
! intmonf        Integrals of monomials over all cells of stencil
!                of given face
! nstenadvv      Number of dual cells in stencil for given vertex
! nexadvv        Number of dual cells in stencil to be fitted exactly
! nadvvmx        Maximum number of dual cells in advection stencil
!                for given vertex
! stenadvv       Stencil for constructing polynomial fit for
!                given vertex
! intmonv        Integrals of monomials over all dual cells of stencil
!                of given vertex
! xgauss         Gauss points for the interval [0,1]
! wgauss         Gauss weights for the interval [0,1]

INTEGER :: degree = 2
LOGICAL :: monotone = .false.
INTEGER :: nmonomial, ngauss
INTEGER :: nadvfmx, nadvvmx
INTEGER, ALLOCATABLE :: nstenadvf(:), nexadvf(:), &
                        stenadvf(:,:), &
                        nstenadvv(:), nexadvv(:), &
                        stenadvv(:,:)
REAL*8, ALLOCATABLE :: intmonf(:,:,:), intmonv(:,:,:)
REAL*8 :: xgauss(3), wgauss(3)


END MODULE advection

! ================================================================

MODULE timestep

! Information related to timestepping
! Default values are set here; they may be overridden
! via the namelist

IMPLICIT NONE

! Time step
REAL*8 :: dt = 300.0d0

! Off-centring parameters
! alpha_v, beta_v         For the velocity used for advection
! alpha_pg, beta_pg       For the pressure gradient term               
REAL*8 :: alpha_v = 0.5d0, alpha_pg = 0.5d0
REAL*8 :: beta_v, beta_pg

! Number of iterations in nonlinear solver
INTEGER :: niter = 4

! Length of integration (in steps)
INTEGER :: nstop = 1

! Frequency of output dumps (in steps)
INTEGER :: noutput = 1

! Frequency of restart dumps (in steps)
INTEGER :: nrestart = 0


! Current time step
INTEGER :: istep

! Current time
REAL*8 :: time


END MODULE timestep

! ================================================================

MODULE channels

! Tidy list of all I/O channels in one place to avoid accidental
! overuse of any channel number

IMPLICIT NONE

INTEGER, PARAMETER :: channml = 20          ! For reading namelists
INTEGER, PARAMETER :: changrid = 25         ! Grid information
INTEGER, PARAMETER :: chanerr = 42          ! Time series of basic error measures
INTEGER, PARAMETER :: chandiag = 43         ! Time series of basic diagnostics
INTEGER, PARAMETER :: chanrefgrd = 26       ! To dump grid coordinates for reference solutions
INTEGER, PARAMETER :: chanrefin = 27        ! To read reference solution
INTEGER, PARAMETER :: chanerrout = 28       ! To write difference from reference solution
! INTEGER, PARAMETER :: chanresin = 50        ! Input channel for restart run
! INTEGER, PARAMETER :: chanresout = 60       ! Restart dumps
INTEGER, PARAMETER :: chandumpm1 = 80       ! Quick look dump file for Matlab (primal grid fields)
INTEGER, PARAMETER :: chandumpm2 = 81       ! Quick look dump file for Matlab (dual grid fields)

END MODULE channels

! ================================================================

PROGRAM mimeticswe

IMPLICIT NONE

! ----------------------------------------------------------------

CALL preliminary
print *,'Done preliminary'

! Dump the grid coordinates for use in generating
! reference solutions
! CALL dumpgrid
! print *,'Done dumpgrid'

! Test the various exterior derivative and Hodge star operators
CALL testop
print *,'Done testop'

! Test the multigrid solver
! CALL testmg
! print *,'Done testmg'

! Generate system matrix for linearized SWEs on the f-sphere
! CALL fsphere
! print *,'Done fsphere'

! Generate system matrix for more general basic state on rotating sphere
! CALL nmodes
! print *,'Done nmodes'

! Test advection scheme
! CALL testadv
! print *,'Done testadv'

! Integrate in time
! CALL integrate
! print *,'Done integration'

! ----------------------------------------------------------------

END PROGRAM mimeticswe

! ================================================================

SUBROUTINE preliminary

! Preliminary calculations and setting up

USE constants
USE grid
USE timestep

IMPLICIT NONE
INTEGER :: iv0
REAL*8 :: lat

! ----------------------------------------------------------------

! Read namelist information
CALL readnml
print *,'Namelists read'

! ----------------------------------------------------------------

! Read in the grid data
CALL readgrid
print *,'Done readgrid'

! ----------------------------------------------------------------

! Allocate array space now that resolution is known
CALL allocateall
print *,'Done allocateall'

! ----------------------------------------------------------------

! Set up dual cell integrals of planetary vorticity
! This assumes jstar is diagonal
DO iv0 = 1, nvert(ngrids)
  lat = vlat(iv0,ngrids)
  planvort2(iv0) = 2.0d0*rotatn*SIN(lat)/jstar(iv0,1,ngrids)
ENDDO
! (An alternative would be to set up a planetary stream function
! and compute its Laplacian - this would ensure global integral
! vanishes, if that matters.)

! ----------------------------------------------------------------

! Set up information needed by advection scheme
CALL setupadv

! ----------------------------------------------------------------

! Set up coefficients needed to compute KE
CALL buildkecoeff

! ----------------------------------------------------------------

! Set up orography
CALL setorog
print *,'Done setorog'

! ----------------------------------------------------------------

! Set up initial conditions
CALL setini
print *,'Done setini'

! ----------------------------------------------------------------

! Set initial time
time = 0.0d0

! Set beta values
beta_v = 1.0d0 - alpha_v
beta_pg = 1.0d0 - alpha_pg

! ----------------------------------------------------------------

END SUBROUTINE preliminary

! ================================================================

SUBROUTINE integrate

! Perform the time integration

USE grid
USE state
USE timestep
IMPLICIT NONE

INTEGER :: nf, nv
REAL*8 :: phi(nfacex)

! ----------------------------------------------------------------

nf = nface(ngrids)
nv = nvert(ngrids)

istep = 0
CALL output
CALL diagnostics

DO istep = 1, nstop
  CALL step
  ! For testing impact of grid scale vorticity noise
  ! if (istep == 720) CALL addcompmode
  time = time + dt
  IF (modulo(istep,noutput) == 0) then
    CALL output
  ENDIF
  IF (modulo(istep,1) == 0) then
    CALL diagnostics
  ENDIF
  ! Compute and output difference from reference solution if
  ! required at this time
  ! CALL diffref
  print *,'Done step ',istep,'   Time = ',time
  print *,' '
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE integrate

! ================================================================

SUBROUTINE readnml

! Read namelists

USE runtype
USE advection
USE timestep
USE channels
IMPLICIT NONE


NAMELIST /rundata/ ygridfile, ygridcoords, runid, lrestart, yresfile
NAMELIST /advectdata/ degree, monotone
NAMELIST /timedata/ dt, niter, alpha_v, alpha_pg, nstop, noutput, nrestart

OPEN(channml,FILE='swenml.in',DELIM='APOSTROPHE')
READ(channml,rundata)
READ(channml,advectdata)
READ(channml,timedata)
CLOSE(channml)


! Compute quantities derived from namelist value
nmonomial = (degree + 1)*(degree + 2)/2
ngauss = 1 + degree/2     ! Integer arithmetic!
IF (ngauss > 3) THEN
  PRINT *,'*** Increase the dimension of xgauss and wgauss to ',ngauss
  PRINT *,'in MODULE advection ***'
  STOP
ENDIF


END SUBROUTINE readnml

! ================================================================

SUBROUTINE setorog

! Set orography field

USE constants
USE grid
USE state
IMPLICIT NONE

INTEGER :: if0
REAL*8 :: phis0, latc, longc, rr0, rr, rr2, long, lat

! ----------------------------------------------------------------

! Zero orography
orog2 = 0.0d0*farea(:,ngrids)

! Williamson et al. test case 5
phis0 = 2000.0d0*gravity
rr0 = pi/9.0d0
latc = pi/6.0d0
longc = 0.5d0*pi   ! Should really be 1.5*pi but this
                   ! makes plotting easier

! DO if0 = 1, nface(ngrids)
  ! Orography at dual vertices
  ! long = flong(if0,ngrids)
  ! lat = flat(if0,ngrids)
  ! Or orography at cell centroid
  ! CALL centroid(if0,long,lat)
  ! rr2 = (long - longc)**2 + (lat - latc)**2
  ! rr = SQRT(MIN(rr0*rr0,rr2))
  ! orog2(if0) = phis0*(1.0 - rr/rr0)*farea(if0,ngrids)
! ENDDO


! ----------------------------------------------------------------

END SUBROUTINE setorog

! ================================================================

SUBROUTINE setini

! Initialize prognostic fields

USE constants
USE grid
USE state
USE errdiag
IMPLICIT NONE

INTEGER :: nf, ne, nv, if0, iv0, ie0, j, jy, nygg
REAL*8 :: u00, phi00, psi(nvertx), utemp(nedgex), &
          alpha, ca, sa, clat, slat, slon
REAL*8, ALLOCATABLE :: psigg(:), hgg(:)
REAL*8 :: l1, l2, &
          lat0, lat1, lat2, en, umen, dygg, psi1, psi2, &
          beta, totvol, totarea, den, &
          cc1, cc2, cc3, cc4, u1, u2, e1, e2, hpert, &
          long, lat, x1, y1, z1, x2, y2, z2, xc, yc, zc, &
          xe, ye, ze, mag, dx, dy, dz
INTEGER :: ic


! ----------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)


! ic = 1   Resting, constant phi
! ic = 2   Williamson et al. test case 2
! ic = 3   Williamson et al. test case 5
! ic = 4   Solid body rotation for normal mode calculation
! ic = 5   Galewsky et al. barotropically unstable jet - specify stream fn
! ic = 6   Galewsky et al. barotropically unstable jet - specify V
ic = 5


IF (ic ==1) THEN
  ! Resting state with uniform phi
  phi2 = 1.0d5*farea(:,ngrids)
  v1 = 0.0d0

ELSEIF (ic == 2 .OR. ic == 3 .OR. ic ==4) THEN

  ! Balanced zonal flow at angle alpha
  IF (ic == 2) THEN
    ! Test case 2
    u00 = 2.0d0*pi*rearth/(12.0d0*86400.0d0)
    phi00 = 2.94d4
  ELSEIF (ic == 3) THEN
    ! Test case 5
    u00 = 20.0d0
    phi00 = 5960.0d0*gravity
  ELSEIF (ic == 4) THEN
    u00 = 0.0d0
    phi00 = 1.0d8
  ENDIF
  alpha = 0.0d0*pi    ! Rotation rate should be zero if alpha <> 0
  ca = COS(alpha)
  sa = SIN(alpha)

  DO iv0 = 1, nv
    clat = COS(vlat(iv0,ngrids))
    slat = SIN(vlat(iv0,ngrids))
    slon = SIN(vlong(iv0,ngrids))
    ! Minus the stream function
    psi(iv0) = u00*rearth*(ca*slat + sa*clat*slon)
  ENDDO
  DO if0 = 1, nf
    ! Initialize at dual vertex...
    ! long = flong(if0,ngrids)
    ! lat = flat(if0,ngrids)
    ! ...or initialize at centroid
    CALL centroid(if0,long,lat)
    clat = COS(lat)
    slat = SIN(lat)
    slon = SIN(long)
    ! Geopotential
    phi2(if0) = phi00 - (rearth*2.0d0*rotatn*u00 + u00*u00)*0.5d0* &
                        (ca*slat + sa*clat*slon)**2
    phi2(if0) = phi2(if0)*farea(if0,ngrids)
  ENDDO

  ! Non-divergent velocity field
  ! U = - D_1 (stream function); sign taken care of above
  CALL Dprimal1(psi,utemp,ngrids,nv,ne)
  ! Corresponding V field
  CALL HodgeHinv(utemp,v1,ngrids,ne)
  print *,'u00 = ',u00
  print *,'Max v = ',MAXVAL(v1/ddist(:,ngrids))

ELSEIF ((ic == 5) .OR. (ic == 6)) THEN

  ! Galewsky et al test case
  nygg = 2*FLOOR(SQRT(REAL(nf)))
  ALLOCATE(hgg(nygg+1), psigg(nygg+1))
  u00 = 80.0
  lat0 = pi/7.0
  lat1 = pi/2.0 - lat0
  en = exp(-4/(lat1 - lat0)**2)
  umen = u00/en
  totvol = 0.0D0
  totarea = 0.0D0
  ! Integrate to tabulate h and psi as functions of geographical
  ! latitude
  dygg = pi/nygg
  hgg(1) = 0.0D0
  psigg(1) = 0.0D0
  DO j = 2, nygg
    l1 = (j-2)*dygg - piby2
    den = (l1 - lat0)*(l1 - lat1)
    IF (den .lt. 0.0D0) THEN
      u1 = umen*exp(1.0D0/den)
    ELSE
      u1 = 0.0D0
    ENDIF
    l2 = (j-1)*dygg - piby2
    den = (l2 - lat0)*(l2 - lat1)
    IF (den .lt. 0.0D0) THEN
      u2 = umen*exp(1.0D0/den)
    ELSE
      u2 = 0.0D0
    ENDIF
    psigg(j) = psigg(j-1) - 0.5d0*(u1 + u2)*dygg
    u1 = u1*(2.0d0*rotatn*SIN(l1) + TAN(l1)*u1/rearth)
    u2 = u2*(2.0d0*rotatn*SIN(l2) + TAN(l2)*u2/rearth)
    hgg(j) = hgg(j-1) - rearth*0.5d0*(u1 + u2)*dygg
    totarea = totarea + COS(l2)*dygg
    totvol = totvol + hgg(j)*COS(l2)*dygg
  ENDDO
  psigg(nygg+1) = psigg(nygg)
  hgg(nygg+1) = hgg(nygg)
  totvol = totvol/(totarea*gravity)
  hgg = hgg + (1.0D4 - totvol)*gravity

  ! Now assign h as a function of geographical latitude
  ! using interpolation from tabulated values
  totvol = 0.00
  totarea = 0.0D0
  DO if0 = 1, nf
    ! l1 = flat(if0,ngrids) + piby2
    CALL centroid(if0,long,lat)
    l1 = lat + piby2
    jy = FLOOR(l1/dygg) + 1
    beta = (l1 - (jy - 1)*dygg)/dygg
    IF (jy == 1 .OR. jy == nygg) THEN
      ! Linear interpolation
      cc2 = 1.0D0 - beta
      cc3 = beta
      phi2(if0) = (cc2*hgg(jy) + cc3*hgg(jy+1))*farea(if0,ngrids)
    ELSE
      ! Cubic interpolation
      cc1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
      cc2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
      cc3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
      cc4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
      phi2(if0) = (cc1*hgg(jy-1) + cc2*hgg(jy) + cc3*hgg(jy+1) + cc4*hgg(jy+2))*farea(if0,ngrids)
    ENDIF
    totarea = totarea + farea(if0,ngrids)
    totvol = totvol + phi2(if0)
  ENDDO
  ! Now calculate velocity components by interpolating
  ! stream function to each vertex
  DO iv0 = 1, nv
    l1 = vlat(iv0,ngrids) + piby2
    jy = FLOOR(l1/dygg) + 1
    beta = (l1 - (jy - 1)*dygg)/dygg
    IF (jy == 1 .OR. jy == nygg) THEN
      ! Linear interpolation
      cc2 = 1.0D0 - beta
      cc3 = beta
      psi(iv0) = cc2*psigg(jy) + cc3*psigg(jy+1)
    ELSE
      ! Cubic interpolation
      cc1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
      cc2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
      cc3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
      cc4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
      psi(iv0) = cc1*psigg(jy-1) + cc2*psigg(jy) + cc3*psigg(jy+1) + cc4*psigg(jy+2)
    ENDIF
  ENDDO
  psi = psi*rearth
  DEALLOCATE(hgg, psigg)
  
  ! Geopotential perturbation
  alpha = 1.0D0/3.0D0
  beta = 1.0D0/15.0D0
  hpert = 120.0D0
  lat2 = 0.5D0*piby2
  DO if0 = 1, nf
    ! l2 = flat(if0,ngrids)
    ! l1 = flong(if0,ngrids)
    CALL centroid(if0,long,lat)
    l2 = lat
    l1 = long
    clat = COS(l2)
    IF (l1 > pi) l1 = l1 - 2.0d0*pi
    e1 = EXP(-(l1/alpha)**2)
    e2 = EXP(-((lat2 - l2)/beta)**2)
    phi2(if0) = phi2(if0) + gravity*hpert*clat*e1*e2*farea(if0,ngrids)
  ENDDO

  IF (ic == 5) THEN
    ! Non-divergent velocity field
    ! U = - D_1 (stream function)
    CALL Dprimal1(psi,utemp,ngrids,nv,ne)
    utemp = -utemp
    ! Corresponding V field
    CALL HodgeHinv(utemp,v1,ngrids,ne)
  ELSE
    DO ie0 = 1, ne
      ! Find location and latitude of dual edge midpoint and
      ! vector displacement along dual edge
      if0 = fnxte(ie0,1,ngrids)
      long = flong(if0,ngrids)
      lat = flat(if0,ngrids)
      CALL ll2xyz(long,lat,x1,y1,z1)
      if0 = fnxte(ie0,2,ngrids)
      long = flong(if0,ngrids)
      lat = flat(if0,ngrids)
      CALL ll2xyz(long,lat,x2,y2,z2)
      xc = 0.5d0*(x1 + x2)
      yc = 0.5d0*(y1 + y2)
      zc = 0.5d0*(z1 + z2)
      mag = 1.0d0/(SQRT(xc*xc + yc*yc + zc*zc))
      xc = xc*mag
      yc = yc*mag
      zc = zc*mag
      dx = x2 - x1
      dy = y2 - y1
      dz = z2 - z1
      CALL xyz2ll(xc,yc,zc,long,lat)
      ! Unit eastward vector
      xe = -yc
      ye = xc
      ze = 0.0d0
      mag = 1.0d0/(SQRT(xe*xe + ye*ye + ze*ze))
      xe = xe*mag
      ye = ye*mag
      ze = ze*mag
      den = (lat - lat0)*(lat - lat1)
      IF (den .lt. 0.0D0) THEN
        u2 = umen*exp(1.0D0/den)
      ELSE
        u2 = 0.0D0
      ENDIF
      v1(ie0) = u2*rearth*(xe*dx + ye*dy + ze*dz)
    ENDDO
  ENDIF
  print *,'u00 = ',u00
  print *,'Max v = ',MAXVAL(v1/ddist(:,ngrids))

ENDIF

! Correct depth to allow for orography, if any
phi2 = phi2 - orog2

! Save initial state for error diagnostics
phi2_init = phi2
v1_init = v1

! Initialize dual grid diagnostic tracers
CALL operR(phi2,xphibar2,ngrids,nf,nv)
CALL Ddual2(v1,xzeta2,ngrids,ne,nv)
xzeta2 = xzeta2 + planvort2


! ----------------------------------------------------------------

END SUBROUTINE setini

! ================================================================

SUBROUTINE step

! Calculations for one time step

USE constants
USE state
USE work
USE helmcoeff
USE timestep
IMPLICIT NONE

INTEGER :: nf, ne, nv, iter, ie0, if1, if2
REAL*8 :: temp

! ----------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)

! Set up Coefficients for Helmholtz problem
phiref = phi2
CALL buildhelm
CALL HodgeI(phi2,phiref,ngrids,nf)
CALL cell2edge(phiref,phirefe,ngrids,nf,ne)

! ----------------------------------------------------------------

print *,'*** Tidy up use of temporary arrays ***'

! First guess for new time level fields
phi2_new = phi2
v1_new = v1

! Calculate old values of primal edge normal velocity fluxes
CALL HodgeH(v1,u1_temp,ngrids,ne)

! Compute divergence - needed to modify swept areas in
! routine primaladvflx
CALL Dprimal2(u1_temp,div2,ngrids,ne,nf)
CALL HodgeI(div2,divfac,ngrids,nf)
divfac = 1.0d0/(1.0d0 + beta_v*dt*divfac)

! Point values of orography plus phi
CALL HodgeI(orog2+beta_pg*phi2,phi0,ngrids,nf)

! Add point values of old KE
CALL findke(v1,b0,ne,nf)
b0 = phi0 + beta_pg*b0

! Construct old dual cell integrals of vorticity, and hence PV
! First relative vorticity
CALL Ddual2(v1,zeta2,ngrids,ne,nv)
! Add planetary contribution
zeta2 = zeta2 + planvort2
! Dual cell geopotential
CALL operR(phi2,phibar2,ngrids,nf,nv)
! Vertex PV
pv = zeta2/phibar2
! Finally the dual cell area integrals of PV
CALL HodgeJinv(pv,pva,ngrids,nv)

! -----------------------------------------------------------------

! Begin iterative solution of nonlinear problem
DO iter = 1, niter

  ! Cell centre value of new phi
  CALL HodgeI(phi2_new,phi0,ngrids,nf)

  ! Find new value of cell centre KE
  CALL findke(v1_new,b0_temp,ne,nf)
  ! Hence the time integral of KE + phi + orog
  b0_temp = (b0 + alpha_pg*(b0_temp + phi0))*dt

  ! And compute dual edge integrals of gradient
  CALL Ddual1(b0_temp,gb1,ngrids,nf,ne)

  ! Construct time averaged velocities
  vbar1 = beta_v*v1 + alpha_v*v1_new
  CALL HodgeH(vbar1,ubar1,ngrids,ne)
  CALL operW(ubar1,vperpbar1,ngrids,ne)
  CALL HodgeH(vperpbar1,uperpbar1,ngrids,ne)

  ! Compute primal grid mass fluxes
  CALL primaladvflx(phi2,mf1)

  ! Construct dual grid mass fluxes
  CALL operW(mf1,mfperp1,ngrids,ne)

  ! Compute dual grid PV fluxes
  CALL dualadvflx(pva,qperp1)

  ! Divergence of mass flux
  CALL Dprimal2(mf1,divmf,ngrids,ne,nf)

  ! Residual in mass equation
  rphi2 = phi2_new - phi2 + divmf

  ! Residual in momentum equation
  rv1 = v1_new - v1 - qperp1 + gb1

  print *,'RMS rphi = ',SQRT(SUM((rphi2/farea(:,ngrids))**2)/nf)
  print *,'RMS rv   = ',SQRT(SUM((rv1/ddist(:,ngrids))**2)/ne)

  ! Build RHS of Helmholtz problem
  CALL HodgeH(rv1,u1_temp,ngrids,ne)
  u1_temp = phirefe*u1_temp
  CALL Dprimal2(u1_temp,rhs,ngrids,ne,nf)
  rhs =  rphi2 - rhs*alpha_v*dt

  ! Solve Helmholtz problem
  CALL mgsolve(phi_inc,rhs,ngrids)

  ! Backsubstitute for v_inc
  CALL HodgeI(phi_inc,b0_temp,ngrids,nf)
  CALL Ddual1(b0_temp,gb1,ngrids,nf,ne)
  v_inc = -rv1 - gb1*alpha_pg*dt

  ! Update prognostic variables
  phi2_new = phi2_new + phi_inc
  v1_new = v1_new + v_inc

ENDDO

! Update diagnostics tracers using last estimate of
! dual mass flux
CALL Ddual2(mfperp1,pva,ngrids,ne,nv)
xphibar2 = xphibar2 + pva
pv = xzeta2/phibar2
CALL HodgeJinv(pv,pva,ngrids,nv)
CALL dualadvflx(pva,qperp1)
CALL Ddual2(qperp1,pva,ngrids,ne,nv)
xzeta2 = xzeta2 + pva

! Rename prognostic variables ready for next step
phi2 = phi2_new
v1 = v1_new


! ----------------------------------------------------------------

END SUBROUTINE step

! ================================================================

SUBROUTINE setupadv

! Set up information needed for advection scheme

USE advection
IMPLICIT NONE

! ----------------------------------------------------------------

! Build stencils for advection
CALL buildadvsten
print *,'  Done buildadvsten'

! Build integrals of monomials over advection stencils
CALL buildintmon
print *,'  Done buildintmon'

! Build coefficients used to construct polynomial fits
! for advection scheme
CALL buildadvcoeff
print *,'  Done buildadvcoeff'

! Set up Gauss points for swept area integrals
CALL setupgauss
print *,'  Done setupgauss'

! ----------------------------------------------------------------

END SUBROUTINE setupadv

! ================================================================

SUBROUTINE allocateall

! Allocate array space that will be needed throughout the code

USE constants
USE state
USE work
USE errdiag
USE helmcoeff
USE advection

IMPLICIT NONE

! ----------------------------------------------------------------

! Arrays in module constants
ALLOCATE(planvort2(nvert(ngrids)))

! Arrays in module state
ALLOCATE(orog2(nface(ngrids)))
ALLOCATE(phi2(nface(ngrids)), v1(nedge(ngrids)), c2(nface(ngrids)))
ALLOCATE(xphibar2(nvert(ngrids)), xzeta2(nvert(ngrids)))

! Arrays in module work
ALLOCATE(vbar1(nedge(ngrids)), ubar1(nedge(ngrids)), &
         vperpbar1(nedge(ngrids)), uperpbar1(nedge(ngrids)), &
         mf1(nedge(ngrids)), mfperp1(nedge(ngrids)), &
         divfac(nface(ngrids)))


ALLOCATE(pva(nvertx), u1_temp(nedgex), &
         b0(nfacex), b0_temp(nfacex), gb1(nedgex), div2(nfacex), zeta2(nvertx), &
         pv(nvertx), qperp1(nedgex), divmf(nfacex), phirefe(nedgex), &
         phi0(nfacex), phi_inc(nfacex), v_inc(nedgex), &
         rphi2(nfacex), rv1(nedgex), rhs(nfacex), &
         phi2_new(nfacex), v1_new(nedgex), &
         phibar2(nvertx))


! Arrays in module errdiag
ALLOCATE(phi2_init(nface(ngrids)), v1_init(nedge(ngrids)))
ALLOCATE(phiexac(nface(ngrids)), pvexac(nvert(ngrids)))
ALLOCATE(phierr(nface(ngrids)), pverr(nvert(ngrids)), v1err(nedge(ngrids)))

! Arrays in module helmholtz
ALLOCATE(phiref(nface(ngrids)))
ALLOCATE(nusq(nedgex,ngrids))
ALLOCATE(helmdiag(nfacex,ngrids))
ALLOCATE(underrel(ngrids))

! Arrays in module advection
! We don't know the maximum stencil size in general until
! we have built the stencils. To avoid convoluted code,
! set maximum stencil size here for known grids and degrees
! of fit and check later once stencils are built.
IF (nefmx == 3) THEN
  ! Triangular primal grid
  IF (degree == 0) THEN
    ! At most one cell for piecewise constant
    nadvfmx = 1
  ELSEIF (degree == 1) THEN
    ! At most 4 cells for piecewise linear
    nadvfmx = 4
  ELSEIF (degree == 2) THEN
    ! At most 10 cells for piecewise quadratic
    nadvfmx = 10
  ELSEIF (degree == 3) THEN
    ! At most 10 cells for piecewise cubic
    nadvfmx = 10
  ELSEIF (degree == 4) THEN
    ! At most 22 cells for piecewise quartic
    nadvfmx = 22
  ELSE
    PRINT *,'Configuration not known.'
    PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
    STOP
  ENDIF
ELSEIF (nefmx == 4) THEN
  ! Quadrilateral primal grid
  IF (degree == 0) THEN
    ! At most one cell for piecewise constant
    nadvfmx = 1
  ELSEIF (degree == 1) THEN
    ! At most 5 cells for piecewise linear
    nadvfmx = 5
  ELSEIF (degree == 2) THEN
    ! At most 9 cells for piecewise quadratic
    nadvfmx = 9
  ELSEIF (degree == 3) THEN
    ! At most 21 cells for piecewise cubic
    nadvfmx = 21
  ELSEIF (degree == 4) THEN
    ! At most 23 cells for piecewise quartic
    nadvfmx = 23
  ELSE
    PRINT *,'Configuration not known.'
    PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
    STOP
  ENDIF
ELSEIF (nefmx == 6) THEN
  ! Hexagonal primal grid
  IF (degree == 0) THEN
    ! At most one cell for piecewise constant
    nadvfmx = 1
  ELSEIF (degree == 1) THEN
    ! At most 7 cells for piecewise linear
    nadvfmx = 7
  ELSEIF (degree == 2) THEN
    ! At most 7 cells for piecewise quadratic
    nadvfmx = 7
  ELSEIF (degree == 3) THEN
    ! At most 13 cells for piecewise cubic
    nadvfmx = 13
  ELSEIF (degree == 4) THEN
    ! At most 19 cells for piecewise quartic
    nadvfmx = 19
  ELSE
    PRINT *,'Configuration not known.'
    PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
    STOP
  ENDIF
ELSE
  PRINT *,'Configuration not known.'
  PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
  STOP
ENDIF

IF (nevmx == 3) THEN
  ! Triangular dual grid
  IF (degree == 0) THEN
    ! At most one cell for piecewise constant
    nadvvmx = 1
  ELSEIF (degree == 1) THEN
    ! At most 4 cells for piecewise linear
    nadvvmx = 4
  ELSEIF (degree == 2) THEN
    ! At most 10 cells for piecewise quadratic
    nadvvmx = 10
  ELSEIF (degree == 3) THEN
    ! At most 10 cells for piecewise cubic
    nadvvmx = 10
  ELSEIF (degree == 4) THEN
    ! At most 22 cells for piecewise quartic
    nadvvmx = 22
  ELSE
    PRINT *,'Configuration not known.'
    PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
    STOP
  ENDIF
ELSEIF (nevmx == 4) THEN
  ! Quadrilateral dual grid
  IF (degree == 0) THEN
    ! At most one cell for piecewise constant
    nadvvmx = 1
  ELSEIF (degree == 1) THEN
    ! At most 5 cells for piecewise linear
    nadvvmx = 5
  ELSEIF (degree == 2) THEN
    ! At most 9 cells for piecewise quadratic
    nadvvmx = 9
  ELSEIF (degree == 3) THEN
    ! At most 21 cells for piecewise cubic
    nadvvmx = 21
  ELSEIF (degree == 4) THEN
    ! At most 21 cells for piecewise quartic
    nadvvmx = 21
  ELSE
    PRINT *,'Configuration not known.'
    PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
    STOP
  ENDIF
ELSEIF (nevmx == 6) THEN
  ! Hexagonal grid
  IF (degree == 0) THEN
    ! At most one cell for piecewise constant
    nadvvmx = 1
  ELSEIF (degree == 1) THEN
    ! At most 7 cells for piecewise linear
    nadvvmx = 7
  ELSEIF (degree == 2) THEN
    ! At most 7 cells for piecewise quadratic
    nadvvmx = 7
  ELSEIF (degree == 3) THEN
    ! At most 13 cells for piecewise cubic
    nadvvmx = 13
  ELSEIF (degree == 4) THEN
    ! At most 19 cells for piecewise quartic
    nadvvmx = 19
  ELSE
    PRINT *,'Configuration not known.'
    PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
    STOP
  ENDIF
ELSE
  PRINT *,'Configuration not known.'
  PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
  STOP
ENDIF

ALLOCATE(nstenadvf(nfacex), nexadvf(nfacex), &
         stenadvf(nfacex,nadvfmx), intmonf(nfacex,nadvfmx,nmonomial))
ALLOCATE(nstenadvv(nvertx), nexadvv(nvertx), &
         stenadvv(nvertx,nadvvmx), intmonv(nvertx,nadvvmx,nmonomial))


! ----------------------------------------------------------------

END SUBROUTINE allocateall

! ================================================================

SUBROUTINE buildhelm

! Build the cell edge values related to phiref needed for
! the multigrid Helmholtz solver at all resolutions

USE grid
USE helmcoeff
USE constants
USE timestep

IMPLICIT NONE
INTEGER :: nf2, ne2, nf1, igrid, if1, if2, if3, ie1, ie2, &
           ixd2, ixh, ixd1, ixi
REAL*8 :: const, temp, temp1(nfacex), temp2(nfacex), &
          cd2, ch, cd1, ci
! real*8 :: phi00, mu2

! ---------------------------------------------------------------

! Useful constant
const = alpha_pg*alpha_v*dt*dt

! Under-relaxation parameter. Should be the best compromise
! between damping small-scale error when Laplacian dominates
! and damping large-scale error when undifferentiated term
! dominates. Depends on grid structure and grid resolution.
! *** There is scope to optimize here ***
DO igrid = 1, ngrids
  IF (nefmx == 6) THEN
    ! Hexagonal grid
    ! mu2 = const*phi00/fareamin(igrid)
    ! underrel(igrid) = 1.0d0
    underrel(igrid) = 0.8d0
    ! underrel(igrid) = (4.0d0*mu2 + 1.0d0)/(6.0d0*mu2 + 1.0d0)
  ELSEIF (nefmx == 4) THEN
    ! Cubic grid
    ! mu2 = const*phi00/fareamin(igrid)
    ! underrel(igrid) = 1.0d0
    underrel(igrid) = 0.8d0
    ! underrel(igrid) = (4.0d0*mu2 + 1.0d0)/(6.0d0*mu2 + 1.0d0)
  ELSE
    print *,'Choose a sensible value for underrel in buildhelm'
  ENDIF
ENDDO


! Construct cell edge values of wave Courant number

! Finest grid
nf2 = nface(ngrids)
ne2 = nedge(ngrids)
temp2 = const*phiref
CALL HodgeI(temp2,temp1,ngrids,nf2)
CALL cell2edge(temp1,nusq(1,ngrids),ngrids,nf2,ne2)

! Coarser grids
DO igrid = ngrids-1, 1, -1
  nf1 = nface(igrid+1)
  temp1(1:nf1) = temp2(1:nf1)
  nf2 = nface(igrid)
  ne2 = nedge(igrid)
  CALL restrict(temp1,nf1,temp2,nf2,igrid)
  CALL HodgeI(temp2,temp1,igrid,nf2)
  CALL cell2edge(temp1,nusq(1,igrid),igrid,nf2,ne2)
ENDDO


! Extract diagonal coefficient of Helmholtz operator
DO igrid = 1, ngrids
  ! Loop over cells
  DO if1 = 1, nface(igrid)

    temp = -1.0d0
    ! Loop over edges of if1 involved in Dprimal2 operator
    DO ixd2 = 1, neoff(if1,igrid)
      ie1 = eoff(if1,ixd2,igrid)
      cd2 = -eoffin(if1,ixd2,igrid)
      ! Loop over edges involved in H operator
      DO ixh = 1, nhsten(ie1,igrid)
        ie2 = hsten(ie1,ixh,igrid)
        ch = hstar(ie1,ixh,igrid)*nusq(ie1,igrid)
        ! Loop over cells involved in Ddual1 operator
        cd1 = 1.0d0
        DO ixd1 = 1, 2
          if2 = fnxte(ie2,ixd1,igrid)
          cd1 = -cd1
          ! Loop over cells in I operator
          DO ixi = 1, nisten(if2,igrid)
            if3 = isten(if2,ixi,igrid)
            ci = istar(if2,ixi,igrid)
            IF (if3 == if1) temp = temp + cd2*ch*cd1*ci
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    helmdiag(if1,igrid) = temp

  ENDDO
ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildhelm

! ================================================================

SUBROUTINE buildadvsten

! Build stencils for advection scheme

USE advection
IMPLICIT NONE

INTEGER :: if0, if1, if2, nlist, list(100), ixs, ixs2, ixn, ixl, &
           iv0, iv1, iv2, ie1
LOGICAL :: alreadysten, alreadylist, lfound

! ---------------------------------------------------------------

! PRIMAL GRID

print *,'*** nexadvf frozen at 1 ***'

! Loop over faces
DO if0 = 1, nface(ngrids)

  ! Start with the face itself
  nstenadvf(if0) = 1
  stenadvf(if0,1) = if0
  ! And demand that this cell should be fitted exactly
  nexadvf(if0) = 1

  ! Iteratively expand stencil until the number of cells
  ! is at least the number of monomials
  DO WHILE (nstenadvf(if0) < nmonomial)

    ! Current cells must be fitted exactly
    ! nexadvf(if0) = nstenadvf(if0)

    ! And look for some new ones
    nlist = 0
    list = 0

    ! Find neighbours of cells currently in the stencil
    DO ixs = 1, nstenadvf(if0)
      if1 = stenadvf(if0,ixs)
      DO ixn = 1, neoff(if1,ngrids)
        if2 = fnxtf(if1,ixn,ngrids)

        ! Is it already in the stencil?
        alreadysten = .false.
        DO ixs2 = 1, nstenadvf(if0)
          IF (if2 == stenadvf(if0,ixs2)) alreadysten = .true.
        ENDDO

        IF (.NOT. alreadysten) THEN
          ! If it's already on the temporary list, flag the fact that
          ! we've seen it more than once
          alreadylist = .false.
          DO ixl = 1, nlist
            IF (if2 == ABS(list(ixl))) THEN
              ! It's already on the list; make a note, and store as
              ! a positive number 
              alreadylist = .true.
              list(ixl) = if2
            ENDIF
          ENDDO
          IF (.NOT. alreadylist) THEN
            ! It's not already on the list; add it as a negative number
            ! to indicate this is the first time
            nlist = nlist + 1
            list(nlist) = -if2
          ENDIF
        ENDIF

      ENDDO
    ENDDO

    ! If we found any that are neighbours more than once then take them
    lfound = .false.
    DO ixl = 1, nlist
      IF (list(ixl) > 0) THEN
        nstenadvf(if0) = nstenadvf(if0) + 1
        stenadvf(if0,nstenadvf(if0)) = list(ixl)
        lfound = .true.
      ENDIF
    ENDDO

    ! Otherwise, take those that are neighbours just once
    IF (.NOT. lfound) THEN
      DO ixl = 1, nlist
        nstenadvf(if0) = nstenadvf(if0) + 1
        stenadvf(if0,nstenadvf(if0)) = ABS(list(ixl))
      ENDDO
    ENDIF

  ENDDO

ENDDO

IF (MAXVAL(nstenadvf) .NE. nadvfmx) THEN
  PRINT *,'nadvfmx = ',nadvfmx,' but MAXVAL(nstenadvf) = ',MAXVAL(nstenadvf)
  STOP
ENDIF


! DUAL GRID

print *,'*** nexadvv frozen at 1 ***'

! Loop over dual cells
DO iv0 = 1, nvert(ngrids)

  ! Start with the dual cell itself
  nstenadvv(iv0) = 1
  stenadvv(iv0,1) = iv0
  ! And demand that this cell should be fitted exactly
  nexadvv(iv0) = 1

  ! Iteratively expand stencil until the number of cells
  ! is at least the number of monomials
  DO WHILE (nstenadvv(iv0) < nmonomial)

    ! Current cells must be fitted exactly
    ! nexadvv(iv0) = nstenadvv(iv0)

    ! And look for some new ones
    nlist = 0
    list = 0

    ! Find neighbours of cells currently in the stencil
    DO ixs = 1, nstenadvv(iv0)
      iv1 = stenadvv(iv0,ixs)
      DO ixn = 1, neofv(iv1,ngrids)
        ! We don't store the vertices neighbouring a given vertex,
        ! so look along edges of the vertex to find neighbours
        ie1 = eofv(iv1,ixn,ngrids)
        iv2 = vofe(ie1,1,ngrids)
        IF (iv2 == iv1) iv2 = vofe(ie1,2,ngrids)

        ! Is it already in the stencil?
        alreadysten = .false.
        DO ixs2 = 1, nstenadvv(iv0)
          IF (iv2 == stenadvv(iv0,ixs2)) alreadysten = .true.
        ENDDO

        IF (.NOT. alreadysten) THEN
          ! If it's already on the temporary list, flag the fact that
          ! we've seen it more than once
          alreadylist = .false.
          DO ixl = 1, nlist
            IF (iv2 == ABS(list(ixl))) THEN
              ! It's already on the list; make a note, and store as
              ! a positive number 
              alreadylist = .true.
              list(ixl) = iv2
            ENDIF
          ENDDO
          IF (.NOT. alreadylist) THEN
            ! It's not already on the list; add it as a negative number
            ! to indicate this is the first time
            nlist = nlist + 1
            list(nlist) = -iv2
          ENDIF
        ENDIF

      ENDDO
    ENDDO

    ! If we found any that are neighbours more than once then take them
    lfound = .false.
    DO ixl = 1, nlist
      IF (list(ixl) > 0) THEN
        nstenadvv(iv0) = nstenadvv(iv0) + 1
        stenadvv(iv0,nstenadvv(iv0)) = list(ixl)
        lfound = .true.
      ENDIF
    ENDDO

    ! Otherwise, take those that are neighbours just once
    IF (.NOT. lfound) THEN
      DO ixl = 1, nlist
        nstenadvv(iv0) = nstenadvv(iv0) + 1
        stenadvv(iv0,nstenadvv(iv0)) = ABS(list(ixl))
      ENDDO
    ENDIF

  ENDDO

ENDDO

IF (MAXVAL(nstenadvv) .NE. nadvvmx) THEN
  PRINT *,'nadvvmx = ',nadvvmx,' but MAXVAL(nstenadvv) = ',MAXVAL(nstenadvv)
  STOP
ENDIF


! ---------------------------------------------------------------

END SUBROUTINE buildadvsten

! ================================================================

SUBROUTINE buildintmon

! Build integrals of monomials over all cells in advection stencils

! *** Note: it may be better to build these on the fly rather than
! store and retrieve them. ***

USE constants
USE advection
IMPLICIT NONE

INTEGER :: if0, if1, if2, ixs, ixt, ie1, iv0, iv1, iv2, igauss, &
           px, py, m

REAL*8 :: long, lat, x0, y0, z0, x1, y1, z1, xn1, yn1, zn1, &
          xn2, yn2, zn2, mag, area, wgt, p(3), xg, yg, zg, &
          s, sinth, costh, xx, yy, &
          fn, x2, y2, z2, x3, y3, z3, x4, y4, z4

! ---------------------------------------------------------------

! PRIMAL GRID

! Initialize to zero
intmonf = 0.0d0

! Loop over all faces
DO if0 = 1, nface(ngrids)

  ! Position vector of centre of face if0
  long = flong(if0,ngrids)
  lat = flat(if0,ngrids)
  CALL ll2xyz(long,lat,x0,y0,z0)

  ! Find direction of first neighbour to establish axes of
  ! Local coordinate system
  if1 = fnxtf(if0,1,ngrids)
  ! Position vector of face if1
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  CALL ll2xyz(long,lat,x1,y1,z1)

  ! Unit normal to plane containing points 0 and 1
  xn1 = y0*z1 - z0*y1
  yn1 = z0*x1 - x0*z1
  zn1 = x0*y1 - y0*x1
  mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
  xn1 = xn1/mag
  yn1 = yn1/mag
  zn1 = zn1/mag

  ! Loop over all cells in stencils
  DO ixs = 1, nstenadvf(if0)

    if2 = stenadvf(if0,ixs)

    ! Centre of this cell
    long = flong(if2,ngrids)
    lat = flat(if2,ngrids)
    CALL ll2xyz(long,lat,x2,y2,z2)

    ! Subdivide cell into triangles
    ! (This is overkill if the cell is itself a triangle)
    DO ixt = 1, neoff(if2,ngrids)
      ie1 = eoff(if2,ixt,ngrids)
      iv1 = vofe(ie1,1,ngrids)
      long = vlong(iv1,ngrids)
      lat = vlat(iv1,ngrids)
      CALL ll2xyz(long,lat,x3,y3,z3)
      iv2 = vofe(ie1,2,ngrids)
      long = vlong(iv2,ngrids)
      lat = vlat(iv2,ngrids)
      CALL ll2xyz(long,lat,x4,y4,z4)
      ! Area of this triangle
      CALL starea2(x2,y2,z2,x3,y3,z3,x4,y4,z4,area)
      wgt = area*rearth*rearth/3.0d0
      ! Loop over Gauss points for this triangle
      DO igauss = 1, 3
        p(1:3) = 1.0d0/6.0d0
        p(igauss) = 4.0d0/6.0d0
        xg = p(1)*x2 + p(2)*x3 + p(3)*x4
        yg = p(1)*y2 + p(2)*y3 + p(3)*y4
        zg = p(1)*z2 + p(2)*z3 + p(3)*z4
        mag = SQRT(xg*xg + yg*yg + zg*zg)
        xg = xg/mag
        yg = yg/mag
        zg = zg/mag
        ! Express Gauss point in terms of face if0's local coordinates
        ! Distance
        xn2 = y0*zg - z0*yg
        yn2 = z0*xg - x0*zg
        zn2 = x0*yg - y0*xg
        mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
        s = ASIN(mag)*rearth
        ! Unit normal
        xn2 = xn2/mag
        yn2 = yn2/mag
        zn2 = zn2/mag
        ! Angle relative to local x-axis
        sinth = x0*(yn1*zn2 - zn1*yn2) &
              + y0*(zn1*xn2 - xn1*zn2) &
              + z0*(xn1*yn2 - yn1*xn2)
        costh = xn1*xn2 + yn1*yn2 + zn1*zn2
        ! Finally obtain local coordinates
        xx = s*costh
        yy = s*sinth
        ! Loop over monomials
        px = 0
        py = 0
        DO m = 1, nmonomial
          fn = (xx**px)*(yy**py)
          intmonf(if0,ixs,m) = intmonf(if0,ixs,m) + wgt*fn
          px = px - 1
          py = py + 1
          IF (px < 0) THEN
            px = py
            py = 0
          ENDIF
        ENDDO      ! End loop over monomials
      ENDDO        ! End loop over Gauss points
    ENDDO          ! End loop over sub-triangles

  ENDDO            ! End loop over cells in stencil

ENDDO              ! End loop over faces


! ---------------------------------------------------------------

! DUAL GRID

! Initialize to zero
intmonv = 0.0d0

! Loop over all dual
DO iv0 = 1, nvert(ngrids)

  ! Position vector of centre of dual cell iv0
  long = vlong(iv0,ngrids)
  lat = vlat(iv0,ngrids)
  CALL ll2xyz(long,lat,x0,y0,z0)

  ! Find direction of first neighbour to establish axes of
  ! Local coordinate system
  ie1 = eofv(iv0,1,ngrids)
  iv1 = vofe(ie1,1,ngrids)
  IF (iv1 == iv0) iv1 = vofe(ie1,2,ngrids)
  ! Position vector of dual cell iv1
  long = vlong(iv1,ngrids)
  lat = vlat(iv1,ngrids)
  CALL ll2xyz(long,lat,x1,y1,z1)

  ! Unit normal to plane containing points 0 and 1
  xn1 = y0*z1 - z0*y1
  yn1 = z0*x1 - x0*z1
  zn1 = x0*y1 - y0*x1
  mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
  xn1 = xn1/mag
  yn1 = yn1/mag
  zn1 = zn1/mag

  ! Loop over all dual cells in stencils
  DO ixs = 1, nstenadvv(iv0)

    iv2 = stenadvv(iv0,ixs)

    ! Centre of this dual cell
    long = vlong(iv2,ngrids)
    lat = vlat(iv2,ngrids)
    CALL ll2xyz(long,lat,x2,y2,z2)

    ! Subdivide cell into triangles
    ! (This is overkill if the cell is itself a triangle)
    DO ixt = 1, neofv(iv2,ngrids)
      ie1 = eofv(iv2,ixt,ngrids)
      if1 = fnxte(ie1,1,ngrids)
      long = flong(if1,ngrids)
      lat = flat(if1,ngrids)
      CALL ll2xyz(long,lat,x3,y3,z3)
      if2 = fnxte(ie1,2,ngrids)
      long = flong(if2,ngrids)
      lat = flat(if2,ngrids)
      CALL ll2xyz(long,lat,x4,y4,z4)
      ! Area of this triangle
      CALL starea2(x2,y2,z2,x3,y3,z3,x4,y4,z4,area)
      wgt = area*rearth*rearth/3.0d0
      ! Loop over Gauss points for this triangle
      DO igauss = 1, 3
        p(1:3) = 1.0d0/6.0d0
        p(igauss) = 4.0d0/6.0d0
        xg = p(1)*x2 + p(2)*x3 + p(3)*x4
        yg = p(1)*y2 + p(2)*y3 + p(3)*y4
        zg = p(1)*z2 + p(2)*z3 + p(3)*z4
        mag = SQRT(xg*xg + yg*yg + zg*zg)
        xg = xg/mag
        yg = yg/mag
        zg = zg/mag
        ! Express Gauss point in terms of dual cell iv0's local coordinates
        ! Distance
        xn2 = y0*zg - z0*yg
        yn2 = z0*xg - x0*zg
        zn2 = x0*yg - y0*xg
        mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
        s = ASIN(mag)*rearth
        ! Unit normal
        xn2 = xn2/mag
        yn2 = yn2/mag
        zn2 = zn2/mag
        ! Angle relative to local x-axis
        sinth = x0*(yn1*zn2 - zn1*yn2) &
              + y0*(zn1*xn2 - xn1*zn2) &
              + z0*(xn1*yn2 - yn1*xn2)
        costh = xn1*xn2 + yn1*yn2 + zn1*zn2
        ! Finally obtain local coordinates
        xx = s*costh
        yy = s*sinth
        ! Loop over monomials
        px = 0
        py = 0
        DO m = 1, nmonomial
          fn = (xx**px)*(yy**py)
          intmonv(iv0,ixs,m) = intmonv(iv0,ixs,m) + wgt*fn
          px = px - 1
          py = py + 1
          IF (px < 0) THEN
            px = py
            py = 0
          ENDIF
        ENDDO      ! End loop over monomials
      ENDDO        ! End loop over Gauss points
    ENDDO          ! End loop over sub-triangles

  ENDDO            ! End loop over dual cells in stencil

ENDDO              ! End loop over dual cells



! ---------------------------------------------------------------

END SUBROUTINE buildintmon

! ================================================================

SUBROUTINE buildadvcoeff

! Build the coefficients used to construct polynomial fits for
! advection scheme. The polynomial fit should be exact for
! nexadvf or nexadvv cells in the stencil, and the mean square
! residual over the other cells should be minimized. This leads to
! to linear problem for the coefficients of the polynomial fit
! and some Lagrange multipliers. Finding the inverse of the matrix
! involved allows us to save the matrix that gives the polynomial
! coefficients in terms of the cell integrals of the advected quantity.

USE advection
IMPLICIT NONE

INTEGER :: if0, iv0, i, ns, ne, nm
REAL*8, ALLOCATABLE :: l(:,:), lt(:,:), ltl(:,:), &
                       m(:,:), minv(:,:), r(:,:), q(:,:)

! ---------------------------------------------------------------

! PRIMAL GRID

! Loop over faces
DO if0 = 1, nface(ngrids)

  ! Determine the stencil size and the number of cells fitted
  ! exactly; hence allocate space for matrices
  ns = nstenadvf(if0)
  ne = nexadvf(if0)
  nm = nmonomial + ne
  ALLOCATE(l(ns,nmonomial), lt(nmonomial,ns), ltl(nmonomial,nmonomial))
  ALLOCATE(m(nm,nm), minv(nm,nm), r(nm,ns), q(nm,ns))

  ! Extract the matrix of integrals of monomials
  l = intmonf(if0,1:ns,1:nmonomial)
  lt = TRANSPOSE(l)
  ltl = MATMUL(lt,l)

  ! Build matrix for linear system
  m(1:nmonomial,1:nmonomial) = ltl
  m(1:nmonomial,nmonomial+1:nm) = lt(1:nmonomial,1:ne)
  m(nmonomial+1:nm,1:nmonomial) = l(1:ne,1:nmonomial)
  m(nmonomial+1:nm,nmonomial+1:nm) = 0.0d0

  ! Invert matrix
  CALL matinv(m,minv,nm)

  ! Matrix relating cell integrals to RHS of linear system
  r(1:nmonomial,1:ns) = lt(1:nmonomial,1:ns)
  r(nmonomial+1:nm,1:ns) = 0.0d0
  DO i = 1,ne
    r(nmonomial+i,i) = 1.0d0
  ENDDO

  ! Matrix giving polynomial coefficients in terms of cell integrals
  q = MATMUL(minv,r)
  
  ! Save the part we need
  lt = q(1:nmonomial,1:ns)
  l = TRANSPOSE(lt)
  intmonf(if0,1:ns,1:nmonomial) = l

  DEALLOCATE(l, lt, ltl, m, minv, r, q)

ENDDO


! DUAL GRID

! Loop over dual cells
DO iv0 = 1, nvert(ngrids)

  ! Determine the stencil size and the number of cells fitted
  ! exactly; hence allocate space for matrices
  ns = nstenadvv(iv0)
  ne = nexadvv(iv0)
  nm = nmonomial + ne
  ALLOCATE(l(ns,nmonomial), lt(nmonomial,ns), ltl(nmonomial,nmonomial))
  ALLOCATE(m(nm,nm), minv(nm,nm), r(nm,ns), q(nm,ns))

  ! Extract the matrix of integrals of monomials
  l = intmonv(iv0,1:ns,1:nmonomial)
  lt = TRANSPOSE(l)
  ltl = MATMUL(lt,l)

  ! Build matrix for linear system
  m(1:nmonomial,1:nmonomial) = ltl
  m(1:nmonomial,nmonomial+1:nm) = lt(1:nmonomial,1:ne)
  m(nmonomial+1:nm,1:nmonomial) = l(1:ne,1:nmonomial)
  m(nmonomial+1:nm,nmonomial+1:nm) = 0.0d0

  ! Invert matrix
  CALL matinv(m,minv,nm)

  ! Matrix relating cell integrals to RHS of linear system
  r(1:nmonomial,1:ns) = lt(1:nmonomial,1:ns)
  r(nmonomial+1:nm,1:ns) = 0.0d0
  DO i = 1,ne
    r(nmonomial+i,i) = 1.0d0
  ENDDO

  ! Matrix giving polynomial coefficients in terms of cell integrals
  q = MATMUL(minv,r)
  
  ! Save the part we need
  lt = q(1:nmonomial,1:ns)
  l = TRANSPOSE(lt)
  intmonv(iv0,1:ns,1:nmonomial) = l

  DEALLOCATE(l, lt, ltl, m, minv, r, q)

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildadvcoeff

! ================================================================

SUBROUTINE setupgauss

! Initialize Gauss points and weights for integrating
! over swept areas in advection scheme

USE advection
IMPLICIT NONE
REAL*8 :: rr3, r3by5

! ---------------------------------------------------------------

! Only lower order versions required so simply assign
! explicitly

IF (ngauss == 1) THEN
  xgauss(1) = 0.5d0
  wgauss(1) = 1.0d0
ELSEIF (ngauss == 2) THEN
  rr3 = 1.0d0/SQRT(3.0D0)
  xgauss(1) = 0.5d0*(1.0d0 - rr3)
  xgauss(2) = 0.5d0*(1.0d0 + rr3)
  wgauss(1) = 0.5d0
  wgauss(2) = 0.5d0
ELSEIF (ngauss == 3) THEN
  r3by5 = SQRT(3.0d0/5.0d0)
  xgauss(1) = 0.5d0*(1.0d0 - r3by5)
  xgauss(2) = 0.5d0
  xgauss(3) = 0.5d0*(1.0d0 + r3by5)
  wgauss(1) = 5.0d0/18.0d0
  wgauss(2) = 4.0d0/9.0d0
  wgauss(3) = 5.0d0/18.0d0
ELSE
  PRINT *,'Need to set up Gauss points and weights for degree = ',degree
  STOP
ENDIF


! ---------------------------------------------------------------

END SUBROUTINE setupgauss

! ================================================================

SUBROUTINE primaladvflx(chi2,flx1)

! Compute advective flux of input field.
! chi2 is the area integral over primal cells of the
! density or concentration of the field to be advected.
! flx1 is the net flux across primal cell edges over one time step.

! It is assumed that the fields u1 and uperp1 are already
! available.

USE constants
USE state
USE work
USE advection
USE timestep
IMPLICIT NONE

REAL*8, INTENT(IN) :: chi2(nfacex)
REAL*8, INTENT(OUT) :: flx1(nedgex)
INTEGER :: nf, ne, ie0, if0, iv1, iv2, m, px, py, ixs, &
           if1, i, j, ns
REAL*8 :: long, lat, x0, y0, z0, &
          x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, xn2, yn2, zn2, &
          s, sinth, costh, xx1, yy1, xx2, yy2, mag, dx, dy, &
          tx, ty, nx, ny, wx, wy, udt, vdt, xg, yg, aswept, &
          flx, poly, fn, xgi, xgj, wgi, wgj, a(nmonomial), p, &
          sgn, temp
REAL*8 :: cn, maxc

! ---------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)

maxc = 0.0d0

! print *,'primaladvflx:  cheaper to loop over cells'

! Loop over edges
DO ie0 = 1, ne

  ! Decide which is the upwind cell
  ! and orient the ends of the edge
  IF (ubar1(ie0) > 0.0d0) THEN
    if0 = fnxte(ie0,1,ngrids)
    iv1 = vofe(ie0,1,ngrids)
    iv2 = vofe(ie0,2,ngrids)
    sgn = 1.0d0
  ELSE
    if0 = fnxte(ie0,2,ngrids)
    iv1 = vofe(ie0,2,ngrids)
    iv2 = vofe(ie0,1,ngrids)
    sgn = -1.0d0
  ENDIF

  ! *** just for diagnostics ***
  cn = ubar1(ie0)*dt/farea(if0,ngrids)
  maxc = MAX(ABS(cn),maxc)

  ! Build the polynomial fit for cell if0
  ns = nstenadvf(if0)
  a = 0.0d0
  DO ixs = 1, ns
    if1 = stenadvf(if0,ixs)
    p = chi2(if1)
    DO m = 1, nmonomial
      a(m) = a(m) + intmonf(if0,ixs,m)*p
    ENDDO
  ENDDO

  ! Centre of upwind cell
  long = flong(if0,ngrids)
  lat = flat(if0,ngrids)
  CALL ll2xyz(long,lat,x0,y0,z0)

  ! Find direction of first neighbour to establish axes of
  ! local coordinate system
  if1 = fnxtf(if0,1,ngrids)
  ! Position vector of face if1
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  CALL ll2xyz(long,lat,x1,y1,z1)

  ! Unit normal to plane containing points 0 and 1
  xn1 = y0*z1 - z0*y1
  yn1 = z0*x1 - x0*z1
  zn1 = x0*y1 - y0*x1
  mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
  xn1 = xn1/mag
  yn1 = yn1/mag
  zn1 = zn1/mag

  ! Position vectors of edge ends,
  ! expressed in local coordinates.
  ! Note: these could be precomputed and stored if that
  ! was more efficient.
  ! First one:
  long = vlong(iv1,ngrids)
  lat = vlat(iv1,ngrids)
  CALL ll2xyz(long,lat,x2,y2,z2)
  ! Distance
  xn2 = y0*z2 - z0*y2
  yn2 = z0*x2 - x0*z2
  zn2 = x0*y2 - y0*x2
  mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
  s = ASIN(mag)*rearth
  ! Unit normal
  xn2 = xn2/mag
  yn2 = yn2/mag
  zn2 = zn2/mag
  ! Angle relative to local x-axis
  sinth = x0*(yn1*zn2 - zn1*yn2) &
        + y0*(zn1*xn2 - xn1*zn2) &
        + z0*(xn1*yn2 - yn1*xn2)
  costh = xn1*xn2 + yn1*yn2 + zn1*zn2
  ! Finally obtain local coordinates
  xx1 = s*costh
  yy1 = s*sinth
  ! Second one:
  long = vlong(iv2,ngrids)
  lat = vlat(iv2,ngrids)
  CALL ll2xyz(long,lat,x2,y2,z2)
  ! Distance
  xn2 = y0*z2 - z0*y2
  yn2 = z0*x2 - x0*z2
  zn2 = x0*y2 - y0*x2
  mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
  s = ASIN(mag)*rearth
  ! Unit normal
  xn2 = xn2/mag
  yn2 = yn2/mag
  zn2 = zn2/mag
  ! Angle relative to local x-axis
  sinth = x0*(yn1*zn2 - zn1*yn2) &
        + y0*(zn1*xn2 - xn1*zn2) &
        + z0*(xn1*yn2 - yn1*xn2)
  costh = xn1*xn2 + yn1*yn2 + zn1*zn2
  ! Finally obtain local coordinates
  xx2 = s*costh
  yy2 = s*sinth

  ! Unit normal and tangent expressed in local coordinates
  dx = xx2 - xx1
  dy = yy2 - yy1
  mag = SQRT(dx*dx + dy*dy)
  tx = dx/mag
  ty = dy/mag
  nx = ty
  ny = -tx

  ! Swept area
  aswept = ubar1(ie0)*dt*divfac(if0)

  ! Displacement in normal and tangential directions
  temp = sgn*dt/ldist(ie0,ngrids)
  udt = ubar1(ie0)*temp
  vdt = uperpbar1(ie0)*temp

  ! Displacement parallel to velocity
  wx = udt*nx + vdt*tx
  wy = udt*ny + vdt*ty

  ! Integrate over swept area
  flx = 0.0d0
  ! Loop over gauss points
  DO i = 1, ngauss
    xgi = xgauss(i)
    wgi = wgauss(i)
    DO j = 1, ngauss
      xgj = xgauss(j)
      wgj = wgauss(j)

      ! Gauss point in local coordinates
      xg = xx1 - xgi*wx + xgj*dx
      yg = yy1 - xgi*wy + xgj*dy

      ! Evaluate polynomial fit
      ! Loop over monomials
      poly = 0.0d0
      px = 0
      py = 0
      DO m = 1, nmonomial
        fn = (xg**px)*(yg**py)
        poly = poly + a(m)*fn
        px = px - 1
        py = py + 1
        IF (px < 0) THEN
          px = py
          py = 0
        ENDIF
      ENDDO
      flx = flx + poly*wgi*wgj
    ENDDO
  ENDDO
  flx1(ie0) = flx*aswept

ENDDO

! print *,'Max Courant No. (primal) = ',maxc

! ---------------------------------------------------------------

END SUBROUTINE primaladvflx

! ================================================================

SUBROUTINE dualadvflx(chi2,flx1)

! Compute advective flux of input field.
! chi2 is the area integral over dual cells of the
! mixing ratio of the field to be advected.
! flx1 is the net flux across dual cell edges over one time step.

! It is assumed that the fields v1, vperp1 and mfperp1 are already
! available.

USE constants
USE state
USE work
USE advection
USE timestep
IMPLICIT NONE

REAL*8, INTENT(IN) :: chi2(nvertx)
REAL*8, INTENT(OUT) :: flx1(nedgex)
INTEGER :: nv, ne, ie0, iv0, if1, if2, m, px, py, ixs, &
           iv1, i, j, ns, ie1
REAL*8 :: long, lat, x0, y0, z0, &
          x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, xn2, yn2, zn2, &
          s, sinth, costh, xx1, yy1, xx2, yy2, mag, dx, dy, &
          tx, ty, nx, ny, wx, wy, udt, vdt, xg, yg, mswept, &
          flx, poly, fn, xgi, xgj, wgi, wgj, a(nmonomial), p, &
          sgn, temp
REAL*8 :: cn, maxc

! ---------------------------------------------------------------

nv = nvert(ngrids)
ne = nedge(ngrids)

maxc = 0.0d0

! Loop over edges
DO ie0 = 1, ne

  ! Decide which is the upwind dual cell
  ! and orient the ends of the edge
  IF (mfperp1(ie0) > 0.0d0) THEN
    iv0 = vofe(ie0,1,ngrids)
    if1 = fnxte(ie0,2,ngrids)
    if2 = fnxte(ie0,1,ngrids)
    sgn = 1.0d0
  ELSE
    iv0 = vofe(ie0,2,ngrids)
    if1 = fnxte(ie0,1,ngrids)
    if2 = fnxte(ie0,2,ngrids)
    sgn = -1.0d0
  ENDIF

  ! *** just for diagnostics ***
  cn = vperpbar1(ie0)*dt*jstar(iv0,1,ngrids)
  maxc = MAX(ABS(cn),maxc)

  ! Build the polynomial fit for dual cell iv0
  ns = nstenadvv(iv0)
  a = 0.0d0
  DO ixs = 1, ns
    iv1 = stenadvv(iv0,ixs)
    p = chi2(iv1)
    DO m = 1, nmonomial
      a(m) = a(m) + intmonv(iv0,ixs,m)*p
    ENDDO
  ENDDO

  ! Centre of upwind cell
  long = vlong(iv0,ngrids)
  lat = vlat(iv0,ngrids)
  CALL ll2xyz(long,lat,x0,y0,z0)

  ! Find direction of first neighbour to establish axes of
  ! local coordinate system
  ie1 = eofv(iv0,1,ngrids)
  iv1 = vofe(ie1,1,ngrids)
  IF (iv1 == iv0) iv1 = vofe(ie1,2,ngrids)
  ! Position vector of dual cell iv1
  long = vlong(iv1,ngrids)
  lat = vlat(iv1,ngrids)
  CALL ll2xyz(long,lat,x1,y1,z1)

  ! Unit normal to plane containing points 0 and 1
  xn1 = y0*z1 - z0*y1
  yn1 = z0*x1 - x0*z1
  zn1 = x0*y1 - y0*x1
  mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
  xn1 = xn1/mag
  yn1 = yn1/mag
  zn1 = zn1/mag

  ! Position vectors of edge ends,
  ! expressed in local coordinates.
  ! Note: these could be precomputed and stored if that
  ! was more efficient.
  ! First one:
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  CALL ll2xyz(long,lat,x2,y2,z2)
  ! Distance
  xn2 = y0*z2 - z0*y2
  yn2 = z0*x2 - x0*z2
  zn2 = x0*y2 - y0*x2
  mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
  s = ASIN(mag)*rearth
  ! Unit normal
  xn2 = xn2/mag
  yn2 = yn2/mag
  zn2 = zn2/mag
  ! Angle relative to local x-axis
  sinth = x0*(yn1*zn2 - zn1*yn2) &
        + y0*(zn1*xn2 - xn1*zn2) &
        + z0*(xn1*yn2 - yn1*xn2)
  costh = xn1*xn2 + yn1*yn2 + zn1*zn2
  ! Finally obtain local coordinates
  xx1 = s*costh
  yy1 = s*sinth
  ! Second one:
  long = flong(if2,ngrids)
  lat = flat(if2,ngrids)
  CALL ll2xyz(long,lat,x2,y2,z2)
  ! Distance
  xn2 = y0*z2 - z0*y2
  yn2 = z0*x2 - x0*z2
  zn2 = x0*y2 - y0*x2
  mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
  s = ASIN(mag)*rearth
  ! Unit normal
  xn2 = xn2/mag
  yn2 = yn2/mag
  zn2 = zn2/mag
  ! Angle relative to local x-axis
  sinth = x0*(yn1*zn2 - zn1*yn2) &
        + y0*(zn1*xn2 - xn1*zn2) &
        + z0*(xn1*yn2 - yn1*xn2)
  costh = xn1*xn2 + yn1*yn2 + zn1*zn2
  ! Finally obtain local coordinates
  xx2 = s*costh
  yy2 = s*sinth

  ! Unit normal and tangent expressed in local coordinates
  dx = xx2 - xx1
  dy = yy2 - yy1
  mag = SQRT(dx*dx + dy*dy)
  tx = dx/mag
  ty = dy/mag
  nx = ty
  ny = -tx

  ! Swept mass
  mswept = mfperp1(ie0)

  ! Displacement in normal and tangential directions
  temp = sgn*dt/ddist(ie0,ngrids)
  udt = vperpbar1(ie0)*temp
  vdt = -vbar1(ie0)*temp

  ! Displacement parallel to velocity
  wx = udt*nx + vdt*tx
  wy = udt*ny + vdt*ty

  ! Integrate over swept area
  flx = 0.0d0
  ! Loop over gauss points
  DO i = 1, ngauss
    xgi = xgauss(i)
    wgi = wgauss(i)
    DO j = 1, ngauss
      xgj = xgauss(j)
      wgj = wgauss(j)

      ! Gauss point in local coordinates
      xg = xx1 - xgi*wx + xgj*dx
      yg = yy1 - xgi*wy + xgj*dy

      ! Evaluate polynomial fit
      ! Loop over monomials
      poly = 0.0d0
      px = 0
      py = 0
      DO m = 1, nmonomial
        fn = (xg**px)*(yg**py)
        poly = poly + a(m)*fn
        px = px - 1
        py = py + 1
        IF (px < 0) THEN
          px = py
          py = 0
        ENDIF
      ENDDO
      flx = flx + poly*wgi*wgj
    ENDDO
  ENDDO
  flx1(ie0) = flx*mswept

ENDDO

!print *,'Max Courant No. (dual) =   ',maxc

! ---------------------------------------------------------------

END SUBROUTINE dualadvflx

! ================================================================

SUBROUTINE buildkecoeff

! Build coefficients used to calculate KE

USE constants
USE grid
IMPLICIT NONE

INTEGER :: if0, if1, if2, ixe

REAL*8 :: long, lat, x0, y0, z0, x1, y1, z1, xn1, yn1, zn1, &
          xn2, yn2, zn2, mag, s, sinth, costh, xx, yy, &
          x2, y2, z2, cx, cy, mm(2,2), rdet, sg

! ---------------------------------------------------------------

! Loop over all cells
DO if0 = 1, nface(ngrids)

  ! Position vector of centre of cell if0
  long = flong(if0,ngrids)
  lat = flat(if0,ngrids)
  CALL ll2xyz(long,lat,x0,y0,z0)

  ! Find direction of first neighbour to establish axes of
  ! Local coordinate system
  if1 = fnxtf(if0,1,ngrids)
  ! Position vector of face if1
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  CALL ll2xyz(long,lat,x1,y1,z1)

  ! Unit normal to plane containing points 0 and 1
  xn1 = y0*z1 - z0*y1
  yn1 = z0*x1 - x0*z1
  zn1 = x0*y1 - y0*x1
  mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
  xn1 = xn1/mag
  yn1 = yn1/mag
  zn1 = zn1/mag

  ! Initialize accumulated matrix
  mm = 0.0d0

  ! Loop over edges of cell
  DO ixe = 1, neoff(if0,ngrids)

    if2 = fnxtf(if0,ixe,ngrids)

    ! Centre of this cell
    long = flong(if2,ngrids)
    lat = flat(if2,ngrids)
    CALL ll2xyz(long,lat,x2,y2,z2)

    ! Express this point in terms of face if0's local coordinates
    ! Distance
    xn2 = y0*z2 - z0*y2
    yn2 = z0*x2 - x0*z2
    zn2 = x0*y2 - y0*x2
    mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
    s = ASIN(mag)*rearth
    ! Unit normal
    xn2 = xn2/mag
    yn2 = yn2/mag
    zn2 = zn2/mag
    ! Angle relative to local x-axis
    sinth = x0*(yn1*zn2 - zn1*yn2) &
          + y0*(zn1*xn2 - xn1*zn2) &
          + z0*(xn1*yn2 - yn1*xn2)
    costh = xn1*xn2 + yn1*yn2 + zn1*zn2
    ! Sign convention for this edge
    sg = -eoffin(if0,ixe,ngrids)
    ! Finally obtain local coordinates
    xx = sg*s*costh
    yy = sg*s*sinth
    ! And accumulate in mm (we don't really need to store
    ! all cpts because of symmetry)
    mm(1,1) = mm(1,1) + xx*xx
    mm(2,1) = mm(2,1) + yy*xx
    mm(1,2) = mm(1,2) + xx*yy
    mm(2,2) = mm(2,2) + yy*yy
    ! Save xx and yy for later
    kecoeff(if0,ixe,1) = xx
    kecoeff(if0,ixe,2) = yy

  ENDDO

  ! Now multiply by mm inverse to get final coeffcients
  DO ixe = 1, neoff(if0,ngrids)

    rdet = 1.0d0/(mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1))
    cx = (mm(2,2)*kecoeff(if0,ixe,1) - mm(1,2)*kecoeff(if0,ixe,2))*rdet
    cy = (mm(1,1)*kecoeff(if0,ixe,2) - mm(2,1)*kecoeff(if0,ixe,1))*rdet
    kecoeff(if0,ixe,1) = cx
    kecoeff(if0,ixe,2) = cy

  ENDDO

ENDDO              ! End loop over cells


END SUBROUTINE buildkecoeff

! ================================================================

SUBROUTINE findke(v,k,ne,nf)

! Compute the kinetic energy per unit mass at cell centres

USE grid
IMPLICIT NONE
INTEGER, INTENT(IN) :: ne, nf
REAL*8, INTENT(IN) :: v(ne)
REAL*8, INTENT(OUT) :: k(nf)
INTEGER :: if0, ixe, ie1
REAL*8 :: ux, uy

! ----------------------------------------------------------------

! Loop over cells
DO if0 = 1, nf

  ! Initialize vector velocity
  ux = 0.0d0
  uy = 0.0d0

  ! Loop over edges
  DO ixe = 1, neoff(if0,ngrids)

    ie1 = eoff(if0,ixe,ngrids)
    ux = ux + v(ie1)*kecoeff(if0,ixe,1)
    uy = uy + v(ie1)*kecoeff(if0,ixe,2)

  ENDDO

  ! Finally compute KE
  k(if0) = 0.5d0*(ux*ux + uy*uy)

ENDDO


END SUBROUTINE findke

! ================================================================

SUBROUTINE Dprimal1(f,df,igrid,nv,ne)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises pointwise values
! at vertices; df comprises integrals of the derivative
! of f along primal cell edges.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nv, ne
REAL*8, INTENT(IN) :: f(nv)
REAL*8, INTENT(OUT) :: df(ne)
INTEGER :: ie1, iv1, iv2

! ----------------------------------------------------------------

DO ie1 = 1, ne
  iv1 = vofe(ie1,1,igrid)
  iv2 = vofe(ie1,2,igrid)
  df(ie1) = f(iv2) - f(iv1)
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Dprimal1

! ================================================================

SUBROUTINE Dprimal2(f,df,igrid,ne,nf)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises integrals along
! primal edges; df comprises integrals of the derivative
! over primal cells.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne, nf
REAL*8, INTENT(IN) :: f(ne)
REAL*8, INTENT(OUT) :: df(nf)
INTEGER :: if1, ix, ie1
REAL*8 :: temp

! ----------------------------------------------------------------

DO if1 = 1, nf
  temp = 0.0d0
  DO ix = 1, neoff(if1,igrid)
    ie1 = eoff(if1,ix,igrid)
    temp = temp - f(ie1)*eoffin(if1,ix,igrid)
  ENDDO
  df(if1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Dprimal2

! ================================================================

SUBROUTINE Ddual1(f,df,igrid,nf,ne)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises pointwise values
! at face centres; df comprises integrals of the derivative
! of f along dual cell edges.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
REAL*8, INTENT(IN) :: f(nf)
REAL*8, INTENT(OUT) :: df(ne)
INTEGER :: ie1, if1, if2

! ----------------------------------------------------------------

DO ie1 = 1, ne
  if1 = fnxte(ie1,1,igrid)
  if2 = fnxte(ie1,2,igrid)
  df(ie1) = f(if2) - f(if1)
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Ddual1

! ================================================================

SUBROUTINE Ddual2(f,df,igrid,ne,nv)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises integrals along
! dual edges; df comprises integrals of the derivative
! over dual cells.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne, nv
REAL*8, INTENT(IN) :: f(ne)
REAL*8, INTENT(OUT) :: df(nv)
INTEGER :: iv1, ix, ie1
REAL*8 :: temp

! ----------------------------------------------------------------

DO iv1 = 1, nv
  temp = 0.0d0
  DO ix = 1, neofv(iv1,igrid)
    ie1 = eofv(iv1,ix,igrid)
    temp = temp + f(ie1)*eofvin(iv1,ix,igrid)
  ENDDO
  df(iv1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Ddual2

! ================================================================

SUBROUTINE HodgeI(f1,f2,igrid,nf)

! Apply the Hodge star I operator that converts face
! integrals f1 to face centre values f2 on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(OUT) :: f2(nf)
INTEGER :: if1, if2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO if1 = 1, nf
  temp = 0.0d0
  DO ix = 1, nisten(if1,igrid)
    if2 = isten(if1,ix,igrid)
    temp = temp + f1(if2)*istar(if1,ix,igrid)
  ENDDO
  f2(if1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE HodgeI

! ================================================================

SUBROUTINE HodgeJ(f1,f2,igrid,nv)

! Apply the Hodge star J operator that converts dual face
! integrals f1 to vertex values f2 on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nv
REAL*8, INTENT(IN) :: f1(nv)
REAL*8, INTENT(OUT) :: f2(nv)
INTEGER :: iv1, iv2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO iv1 = 1, nv
  temp = 0.0d0
  DO ix = 1, njsten(iv1,igrid)
    iv2 = jsten(iv1,ix,igrid)
    temp = temp + f1(iv2)*jstar(iv1,ix,igrid)
  ENDDO
  f2(iv1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE HodgeJ

! ================================================================

SUBROUTINE HodgeH(f1,f2,igrid,ne)

! Apply the Hodge star H operator that converts dual edge
! integrals (circulations) f1 to primal edge integrals (fluxes) f2
! on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: ie1, ie2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO ie1 = 1, ne
  temp = 0.0d0
  DO ix = 1, nhsten(ie1,igrid)
    ie2 = hsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*hstar(ie1,ix,igrid)
  ENDDO
  f2(ie1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE HodgeH

! ================================================================

SUBROUTINE HodgeIinv(f1,f2,igrid,nf)

! Apply the inverse Hodge star operator I^{-1} that converts
! primal cell area integrals f1 to primal cell point values f2
! on grid igrid.

USE grid
IMPLICIT NONE

INTEGER :: niter = 10
INTEGER, INTENT(IN) :: igrid, nf
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(OUT) :: f2(nf)
INTEGER :: if1, iter
REAL*8 :: temp(nf)

! ----------------------------------------------------------------

! First guess based on diagonal I
DO if1 = 1, nf
  f2(if1) = f1(if1)/istar(if1,1,igrid)
ENDDO

IF (nismx > 1) THEN
  ! I is not diagonal, so use Jacobi iteration to invert
  DO iter = 1, niter
    CALL HodgeI(f2,temp,igrid,nf)
    temp = f1 - temp
    DO if1 = 1, nf
      f2(if1) = f2(if1) + temp(if1)/istar(if1,1,igrid)
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------

END SUBROUTINE HodgeIinv

! ================================================================

SUBROUTINE HodgeJinv(f1,f2,igrid,nv)

! Apply the inverse Hodge star operator J^{-1} that converts
! dual cell area integrals f1 to dual cell point values f2
! on grid igrid.

USE grid
IMPLICIT NONE

INTEGER :: niter = 10
INTEGER, INTENT(IN) :: igrid, nv
REAL*8, INTENT(IN) :: f1(nv)
REAL*8, INTENT(OUT) :: f2(nv)
INTEGER :: iv1, iter
REAL*8 :: temp(nv)

! ----------------------------------------------------------------

! First guess based on diagonal J
DO iv1 = 1, nv
  f2(iv1) = f1(iv1)/jstar(iv1,1,igrid)
ENDDO

IF (njsmx > 1) THEN
  ! J is not diagonal, so use Jacobi iteration to invert
  DO iter = 1, niter
    CALL HodgeJ(f2,temp,igrid,nv)
    temp = f1 - temp
    DO iv1 = 1, nv
      f2(iv1) = f2(iv1) + temp(iv1)/jstar(iv1,1,igrid)
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------

END SUBROUTINE HodgeJinv

! ================================================================

SUBROUTINE HodgeHinv(f1,f2,igrid,ne)

! Apply the inverse Hodge star operator H^{-1} that converts primal edge
! integrals (fluxes) f1 to dual edge integrals (circulations) f2
! on grid igrid.

USE grid
IMPLICIT NONE

INTEGER :: niter = 20
INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: ie1, iter
REAL*8 :: temp(ne)

! ----------------------------------------------------------------

! First guess based on diagonal H
DO ie1 = 1, ne
  f2(ie1) = f1(ie1)/hstar(ie1,1,igrid)
ENDDO

IF (nhsmx > 1) THEN
  ! H is not diagonal, so use Jacobi iteration to invert
  print *,'*** HodgeHinv:  think about number of iterations needed...'
  DO iter = 1, niter
    CALL HodgeH(f2,temp,igrid,ne)
    temp = f1 - temp
    DO ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/hstar(ie1,1,igrid)
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------

END SUBROUTINE HodgeHinv

! ================================================================

SUBROUTINE operW(f1,f2,igrid,ne)

! Apply the W operator:
! given fluxes f1 across primal edges, construct
! the rotated fluxes across dual edges f2, on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: if1, ie1, ie2, ix, ix1, ix2, ixv, ne1
REAL*8 :: ss, w

! ----------------------------------------------------------------

! Initialize to zero
f2 = 0.0d0

! Loop over faces
DO if1 = 1, nface(igrid)
  ne1 = neoff(if1,igrid)
  ! For each edge of this face
  DO ix1 = 1, ne1
    ss = -0.5
    ie1 = eoff(if1,ix1,igrid)
    ! Find the contribution to f2 from every other
    ! edge of this face
    DO ix = 0, ne1 - 2
      ixv = MOD(ix1 + ix - 1,ne1) + 1
      ix2 = MOD(ix1 + ix,ne1) + 1
      ie2 = eoff(if1,ix2,igrid)
      ss = ss + rcoeff(if1,ixv,igrid)
      w = -ss*eoffin(if1,ix1,igrid)*eoffin(if1,ix2,igrid)
      f2(ie1) = f2(ie1) + w*f1(ie2)
    ENDDO
  ENDDO
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE operW

! ================================================================

SUBROUTINE operR(f1,f2,igrid,nf,nv)

! Apply the R operator:
! given face integrals f1 on primal faces, map to dual cell
! integrals f2, on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, nv
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(OUT) :: f2(nv)
INTEGER :: if1, iv1, ix1, ne1

! ----------------------------------------------------------------

! Initialize to zero
f2 = 0.0d0

! Loop over faces
DO if1 = 1, nface(igrid)
  ne1 = neoff(if1,igrid)
  ! Share out this face's contributions to its surrounding vertices
  DO ix1 = 1, ne1
    iv1 = voff(if1,ix1,igrid)
    f2(iv1) = f2(iv1) + f1(if1)*rcoeff(if1,ix1,igrid)
  ENDDO
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE operR

! ================================================================

SUBROUTINE perp(f1,f2,igrid,ne)

! Apply the perp operator:
! given circulations f1 along dual edges, construct
! the rotated fluxes across dual edges f2, on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
REAL*8 :: temp(ne)

! ----------------------------------------------------------------

! First construct fluxes across primal edges
CALL HodgeH(f1,temp,igrid,ne)

! Then apply W operator to rotate
CALL operW(temp,f2,igrid,ne)

! ----------------------------------------------------------------

END SUBROUTINE perp

! ================================================================

SUBROUTINE graddiv(f1,f2,igrid,ne,nf)

! Calculate the grad-div component of the vector Laplacian
! operator applied to circulations

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne, nf
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
REAL*8 :: temp1(ne), temp2(nf), temp3(nf)

! ----------------------------------------------------------------

! Build up from simpler operators
CALL HodgeH(f1,temp1,igrid,ne)
CALL Dprimal2(temp1,temp2,igrid,ne,nf)
CALL HodgeI(temp2,temp3,igrid,nf)
CALL Ddual1(temp3,f2,igrid,nf,ne)

! ----------------------------------------------------------------

END SUBROUTINE graddiv

! ================================================================

SUBROUTINE curlcurl(f1,f2,igrid,ne,nv)

! Calculate the curlcurl component of the vector Laplacian
! operator applied to circulations

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne, nv
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
REAL*8 :: temp1(nv), temp2(nv), temp3(ne)

! ----------------------------------------------------------------

! Build up from simpler operators
CALL Ddual2(f1,temp1,igrid,ne,nv)
CALL HodgeJ(temp1,temp2,igrid,nv)
CALL Dprimal1(temp2,temp3,igrid,nv,ne)
CALL HodgeHinv(temp3,f2,igrid,ne)
f2 = -f2


! ----------------------------------------------------------------


END SUBROUTINE curlcurl

! ================================================================

SUBROUTINE cell2edge(f1,f2,igrid,nf,ne)

! Simple interpolation from cell centre point values to
! edges

USE grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: ie1, if1, if2


DO ie1 = 1, ne
  if1 = fnxte(ie1,1,igrid)
  if2 = fnxte(ie1,2,igrid)
  f2(ie1) = 0.5d0*(f1(if1) + f1(if2))
ENDDO


END SUBROUTINE cell2edge

! ================================================================

SUBROUTINE restrict(f1,nf1,f2,nf2,igrid)

! To perform the restriction operation needed for a multigrid solver.
! Restrict field f1 from grid igrid + 1 to grid igrid and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms).

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: nf1, nf2, igrid
REAL*8, INTENT(IN) :: f1(nf1)
REAL*8, INTENT(OUT) :: f2(nf2)
INTEGER :: if1, if2, ix
REAL*8 :: wgt

! Safety check
IF (nf2 .ne. nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine restrict'
  STOP
ENDIF

DO if2 = 1, nf2
  f2(if2) = 0.0d0
  DO ix = 1, ninj(if2,igrid)
    if1 = injsten(if2,ix,igrid)
    wgt = injwgt(if2,ix,igrid)
    f2(if2) = f2(if2) + wgt*f1(if1)
  ENDDO
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE restrict

! ================================================================

SUBROUTINE prolong(f2,nf2,f1,nf1,igrid)

! To perform the prolongation operation needed for a multigrid solver.
! Prolong field f2 from grid igrid to grid igrid + 1 and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms); so f2 must be converted to point values
! (zero-form) first and f1 must be converted back to an area-integral
! (2-form) at the end. The prolongation operator is the adjoint of
! the restriction operator, so uses the same stencil and weights.

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: nf1, nf2, igrid
REAL*8, INTENT(IN) :: f2(nf2)
REAL*8, INTENT(OUT) :: f1(nf1)
INTEGER :: if1, if2, ix, igridp
REAL*8 :: wgt, f2if2, temp1(nf1), temp2(nf2)

! Safety check
IF (nf2 .ne. nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine prolong'
  STOP
ENDIF

igridp = igrid + 1
CALL HodgeI(f2,temp2,igrid,nf2)
temp1 = 0.0d0
DO if2 = 1, nf2
  f2if2 = temp2(if2)
  DO ix = 1, ninj(if2,igrid)
    if1 = injsten(if2,ix,igrid)
    wgt = injwgt(if2,ix,igrid)
    temp1(if1) = temp1(if1) + wgt*f2if2
  ENDDO
ENDDO
CALL HodgeIinv(temp1,f1,igridp,nf1)

! ----------------------------------------------------------------

END SUBROUTINE prolong

! ================================================================

SUBROUTINE helmholtz(f,hf,igrid,nf,ne)

! To apply the Helmholtz operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

USE helmcoeff

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
REAL*8, INTENT(IN) :: f(nf)
REAL*8, INTENT(OUT) :: hf(nf)
REAL*8 :: temp1(nf), temp2(ne), temp3(ne)


CALL HodgeI(f,temp1,igrid,nf)
CALL Ddual1(temp1,temp2,igrid,nf,ne)
CALL HodgeH(temp2,temp3,igrid,ne)
temp2 = temp3*nusq(:,igrid)
CALL Dprimal2(temp2,hf,igrid,ne,nf)
hf = hf - f


END SUBROUTINE helmholtz

! ================================================================

SUBROUTINE residual(f,rhs,res,igrid,nf,ne)

! Compute the residual res in the helmholtz equation on grid igrid
! when f is the input field and rhs is the right hand side. Note that
! f, rhs and res are area integrals (2-forms).

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
REAL*8, INTENT(IN) :: f(nf), rhs(nf)
REAL*8, INTENT(OUT) :: res(nf)


CALL helmholtz(f,res,igrid,nf,ne)
res = rhs - res


END SUBROUTINE residual

! ================================================================

SUBROUTINE relax(f,rhs,igrid,nf,ne,niter)

! To carry out niter Jacobi relaxation iterations for the multigrid
! solver on grid igrid

USE helmcoeff

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne, niter
REAL*8, INTENT(IN) :: rhs(nf)
REAL*8, INTENT(INOUT) :: f(nf)
REAL*8, ALLOCATABLE :: res(:), finc(:)
REAL*8 :: u
INTEGER :: iter

ALLOCATE(res(nf), finc(nf))

u = underrel(igrid)
DO iter = 1, niter
  CALL residual(f,rhs,res,igrid,nf,ne)
  finc = res/helmdiag(1:nf,igrid)
  f = f + u*finc
ENDDO

DEALLOCATE(res, finc)

END SUBROUTINE relax

! ================================================================

SUBROUTINE mgsolve(phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 nusq H Ddual1 I phi - phi = RHS
!
! using full multigrid algorithm.
! Coefficients are contained in module helmholtz.

USE grid

IMPLICIT NONE

! Numbers of iterations on coarsest grid and other grids
INTEGER, PARAMETER :: niterc = 10, niter = 2, npass = 1

INTEGER, INTENT(IN) :: ng
REAL*8, INTENT(IN) :: rr(nfacex)
REAL*8, INTENT(OUT) :: phi(nfacex)

INTEGER :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
REAL*8, ALLOCATABLE :: ff(:,:), rf(:,:)
REAL*8 :: temp1(nfacex)

! ------------------------------------------------------------------------

! Allocate space on all grids
ALLOCATE(ff(nfacex,ngrids),rf(nfacex,ngrids))

! One pass should be enough. Warn user if npass is set to
! some other value for testing
IF (npass > 1) PRINT *,'mgsolve: npass = ',npass

! ------------------------------------------------------------------------

! Initialize solution to zero
phi = 0.0d0

! For diagnostics
! nf1 = nface(ngrids)
! ne1 = nedge(ngrids)
! CALL residual(phi,rr,temp1,ngrids,nf1,ne1)
! print *,'Pass ',0,'  RMS residual = ',SQRT(SUM(temp1*temp1)/nf1)

DO ipass = 1, npass

  ! Initialize rhs as residual using latest estimate
  IF (ipass == 1) THEN
    ! No need to do the calculation
    rf(:,ngrids) = rr(:)
  ELSE
    nf1 = nface(ngrids)
    ne1 = nedge(ngrids)
    CALL residual(phi,rr,rf(1,ngrids),ngrids,nf1,ne1)
  ENDIF

  ! Initialize correction to solution on all grids to zero
  ff = 0.0d0

  ! Restrict right hand side to each grid in the hierarchy
  DO igrid = ngrids-1, ngrids-ng+1, -1
    igridp = igrid + 1
    nf1 = nface(igridp)
    nf2 = nface(igrid)
    CALL restrict(rf(1,igridp),nf1,rf(1,igrid),nf2,igrid)
  ENDDO

  ! Iterate to convergence on coarsest grid
  igrid = ngrids-ng+1
  nf1 = nface(igrid)
  ne1 = nedge(igrid)
  ff(1:nf1,igrid) = 0.0d0
  CALL relax(ff(1,igrid),rf(1,igrid),igrid,nf1,ne1,niterc)

  ! Sequence of growing V-cycles
  DO igridp = ngrids-ng+2, ngrids

    igrid = igridp - 1
    nf1 = nface(igridp)
    ne1 = nedge(igridp)
    nf2 = nface(igrid)
    ne2 = nedge(igrid)

    ! Prolong solution to grid igridp
    ! and execute one V-cycle starting from grid igridp

    ! Prolong
    CALL prolong(ff(1,igrid),nf2,ff(1,igridp),nf1,igrid)

    ! Descending part of V-cycle
    DO jgrid = igrid, ngrids-ng+1, -1
    
      jgridp = jgrid + 1
      nf1 = nface(jgridp)
      ne1 = nedge(jgridp)
      nf2 = nface(jgrid)
      ne2 = nedge(jgrid)

      ! Relax on grid jgridp
      CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

      ! Calculate residual on jgridp
      CALL residual(ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)
  
      ! Restrict residual to jgrid
      CALL restrict(temp1,nf1,rf(1,jgrid),nf2,jgrid)

      ! Set correction first guess to zero on grid jgrid-1
      ff(1:nf2,jgrid) = 0.0d0

    ENDDO

    ! Iterate to convergence on coarsest grid
    jgrid = ngrids-ng+1
    nf1 = nface(jgrid)
    ne1 = nedge(jgrid)
    ff(1:nf1,jgrid) = 0.0d0
    CALL relax(ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

    ! Ascending part of V-cycle
    DO jgrid = ngrids-ng+1, igrid

      jgridp = jgrid + 1
      igrid = igrid - 1
      nf1 = nface(jgridp)
      ne1 = nedge(jgridp)
      nf2 = nface(jgrid)
      ne2 = nedge(jgrid)

      ! Prolong correction to jgridp
      CALL prolong(ff(1,jgrid),nf2,temp1,nf1,jgrid)

      ! Add correction to solution on jgridp
      ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

      ! Relax on grid jgridp
      CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

    ENDDO

  ENDDO

  ! Add correction to phi
  phi = phi + ff(:,ngrids)

  ! For diagnostics
  ! nf1 = nface(ngrids)
  ! ne1 = nedge(ngrids)
  ! CALL residual(phi,rr,temp1,ngrids,nf1,ne1)
  ! print *,'Pass ',ipass,'  RMS residual = ',SQRT(SUM(temp1*temp1)/nf1)

ENDDO


DEALLOCATE(ff,rf)

! ----------------------------------------------------------------

END SUBROUTINE mgsolve

! ================================================================

SUBROUTINE readgrid

! To allocate array space for the grid information in module grid
! and to read the information from file

USE runtype
USE constants
USE grid
USE channels

IMPLICIT NONE

INTEGER :: if0, ie0, iv0, igrid, ix, if1, if2, iv1, iv2, ie1

! ----------------------------------------------------------------

! Open file for reading
OPEN(changrid,FILE=ygridfile,FORM='UNFORMATTED')

! First read ngrids
READ(changrid) ngrids


! Allocate nface, nedge, nvert
ALLOCATE(nface(ngrids), nedge(ngrids), nvert(ngrids))

! Read numbers of faces, edges and vertices on each grid
READ(changrid) nface
READ(changrid) nedge
READ(changrid) nvert

! Find maximum values in order to allocated subsequent arrays
nfacex = MAXVAL(nface)
nedgex = MAXVAL(nedge)
nvertx = MAXVAL(nvert)

! Allocate neoff, neofv
ALLOCATE(neoff(nfacex,ngrids), neofv(nvertx,ngrids))

! Read the numbers of edges of each face and vertex on each grid
neoff = 0
neofv = 0
READ(changrid) ((neoff(if0,igrid),          &
                    if0 = 1, nface(igrid)), &
                    igrid = 1, ngrids)
READ(changrid) ((neofv(iv0,igrid),          &
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
READ(changrid) (((fnxtf(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((eoff(if0,ix,igrid),           &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((voff(if0,ix,igrid),           &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((fnxte(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
READ(changrid) (((vofe(ie0,ix,igrid),           &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
READ(changrid) (((fofv(iv0,ix,igrid),           &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((eofv(iv0,ix,igrid),           &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)


! Allocate the geometrical information arrays
ALLOCATE(flong(nfacex,ngrids), flat(nfacex,ngrids), &
         vlong(nvertx,ngrids), vlat(nvertx,ngrids), &
         farea(nfacex,ngrids), &
         ldist(nedgex,ngrids), ddist(nedgex,ngrids), &
         fareamin(ngrids))

! Read the geometrical information arrays
READ(changrid) ((flong(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((flat(if0,igrid),               &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((vlong(iv0,igrid),              &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((vlat(iv0,igrid),               &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((farea(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((ldist(ie0,igrid),              &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((ddist(ie0,igrid),              &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)

! Dimensionalize
farea = farea*rearth*rearth
ldist = ldist*rearth
ddist = ddist*rearth

! Determine smallest face area on each grid
DO igrid = 1, ngrids
  fareamin(igrid) = MINVAL(farea(1:nface(igrid),igrid))
ENDDO


! Allocate arrays for size of operator stencils
ALLOCATE(nisten(nfacex,ngrids), njsten(nvertx,ngrids), &
         nhsten(nedgex,ngrids), nrsten(nfacex,ngrids))

! Read the sizes of the operator stencils on each grid
nisten = 0
njsten = 0
nhsten = 0
nrsten = 0
READ(changrid) ((nisten(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((njsten(iv0,igrid),             &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((nhsten(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((nrsten(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)

! Find maximum values in order to allocate subsequent arrays
nismx = MAXVAL(nisten)
njsmx = MAXVAL(njsten)
nhsmx = MAXVAL(nhsten)
nrsmx = MAXVAL(nrsten)

! Allocate arrays for operator stencils and coefficients
ALLOCATE(isten(nfacex,nismx,ngrids), jsten(nvertx,njsmx,ngrids), &
         hsten(nedgex,nhsmx,ngrids), rsten(nfacex,nrsmx,ngrids))
ALLOCATE(istar(nfacex,nismx,ngrids), jstar(nvertx,njsmx,ngrids), &
         hstar(nedgex,nhsmx,ngrids), rcoeff(nfacex,nrsmx,ngrids), &
         kecoeff(nfacex,nefmx,2))

! Read the operator stencils and coefficients
READ(changrid) (((isten(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nismx),            &
                     igrid = 1, ngrids)
READ(changrid) (((jsten(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((hsten(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((rsten(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((istar(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nismx),            &
                     igrid = 1, ngrids)
READ(changrid) (((jstar(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((hstar(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((rcoeff(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)

! Dimensionalize
istar = istar/(rearth*rearth)
jstar = jstar/(rearth*rearth)


! Construct the tables eoffin and eofvin
ALLOCATE(eoffin(nfacex,nefmx,ngrids), eofvin(nvertx,nevmx,ngrids))
DO igrid = 1, ngrids

  DO if1 = 1, nface(igrid)
    DO ix = 1, neoff(if1,igrid)
      ie1 = eoff(if1,ix,igrid)
      if2 = fnxte(ie1,1,igrid)
      IF (if1 == if2) THEN
        ! This edge points out of face if1
	eoffin(if1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points into face if1
	eoffin(if1,ix,igrid) = 1.0d0
      ENDIF
    ENDDO
  ENDDO

  DO iv1 = 1, nvert(igrid)
    DO ix = 1, neofv(iv1,igrid)
      ie1 = eofv(iv1,ix,igrid)
      iv2 = vofe(ie1,1,igrid)
      IF (iv1 == iv2) THEN
        ! This edge points away from vertex iv1
	eofvin(iv1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points towards vertex iv1
	eofvin(iv1,ix,igrid) = 1.0d0
      ENDIF
    ENDDO
  ENDDO

ENDDO


! Allocate array for size of restriction stencil
ALLOCATE(ninj(nfacex,ngrids-1))

! Read the size of the restriction stencil on each grid
ninj = 0
READ(changrid) ((ninj(if0,igrid),              &
                    if0 = 1, nface(igrid)),    &
                    igrid = 1, ngrids-1)

! Find maximum value in order to allocate subsequent arrays
ninjmx = MAXVAL(ninj)

! Allocate arrays for restriction stencils and weights
ALLOCATE(injsten(nfacex,ninjmx,ngrids-1))
ALLOCATE(injwgt(nfacex,ninjmx,ngrids-1))

! Read the restriction stencil and weights
READ(changrid) (((injsten(if0,ix,igrid),       &
                    if0 = 1, nface(igrid)),    &
                    ix = 1, ninjmx),           &
                    igrid = 1, ngrids-1)
READ(changrid) (((injwgt(if0,ix,igrid),        &
                    if0 = 1, nface(igrid)),    &
                    ix = 1, ninjmx),           &
                    igrid = 1, ngrids-1)

! -------------------------------------------------------------------

END SUBROUTINE readgrid

!     ===============================================================
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
      SUBROUTINE XYZ2LL(X,Y,Z,LONG,LAT)
!
!     To convert cartesian coordinates to longitude and latitude
!
      IMPLICIT NONE
!   
      REAL*8 X,Y,Z,LONG,LAT,PI,TLN,TLT,R
!
!     -------------------------------------------------------------------
!
      PI=4.0D0*ATAN(1.0D0)
!
      IF (X.EQ.0.0D0) THEN
        IF (Y.GE.0.0D0) THEN
          LONG=0.5D0*PI
        ELSE
          LONG=1.5D0*PI
        ENDIF
      ELSE
        TLN=Y/X
        LONG=ATAN(TLN)
        IF (X.LT.0.0D0) THEN
          LONG=LONG+PI
        ENDIF
        IF (LONG.LT.0.0D0) THEN
          LONG=LONG+2.0D0*PI
        ENDIF
      ENDIF
!
      R=SQRT(X*X+Y*Y)
      IF (R.EQ.0.0D0) THEN
        IF (Z.GT.0.0D0) THEN
          LAT=0.5D0*PI
        ELSE
          LAT=-0.5D0*PI
        ENDIF
      ELSE
        TLT=Z/R
        LAT=ATAN(TLT)
      ENDIF
!
!     --------------------------------------------------------------------
!
      RETURN
      END
!
!     ====================================================================
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

SUBROUTINE centroid(if0,long,lat)

! Find the centroid of cell if0 on grid ngrids

USE grid
IMPLICIT NONE
INTEGER, INTENT(IN) :: if0
REAL*8, INTENT(OUT) :: long, lat
INTEGER :: ixe, ie1, iv1, iv2
REAL*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          xc, yc, zc, a, aby3, mag

! -----------------------------------------------------------------------

! Coordinates of centre of face
long1 = flong(if0,ngrids)
lat1 = flat(if0,ngrids)
CALL ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the face
! Hence find area of face and centroid
xc = 0.0d0
yc = 0.0d0
zc = 0.0d0
DO ixe = 1, neoff(if0,ngrids)
  ie1 = eoff(if0,ixe,ngrids)
  iv1 = vofe(ie1,1,ngrids)
  iv2 = vofe(ie1,2,ngrids)
  long1 = vlong(iv1,ngrids)
  lat1 = vlat(iv1,ngrids)
  CALL ll2xyz(long1,lat1,x1,y1,z1)
  long1 = vlong(iv2,ngrids)
  lat1 = vlat(iv2,ngrids)
  CALL ll2xyz(long1,lat1,x2,y2,z2)
  CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
  aby3 = a/3.0d0
  xc = xc + (x0 + x1 + x2)*aby3
  yc = yc + (y0 + y1 + y2)*aby3
  zc = zc + (z0 + z1 + z2)*aby3
ENDDO
mag = SQRT(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag
CALL xyz2ll(xc,yc,zc,long,lat)

! -----------------------------------------------------------------------

END SUBROUTINE centroid

!     ===================================================================

SUBROUTINE matinv(a,ainv,n)

! Invert an n x n matrix a by Gaussian elimination and put the
! result in matrix ainv. No fancy pivoting or checking for
! divide by zero!! Upon exit a should be the n x n identity matrix.

IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL*8, INTENT(INOUT) :: a(n,n)
REAL*8, INTENT(OUT) :: ainv(n,n)
INTEGER :: i, j
REAL*8 :: temp

! -----------------------------------------------------------------------

! Initialize ainv to the identify matrix
ainv = 0.0d0
DO i = 1, n
  ainv(i,i) = 1.0d0
ENDDO

! Eliminate subdiagonal terms
DO j = 1, n
  temp = a(j,j)
  a(j,j:n) = a(j,j:n)/temp
  ainv(j,:) = ainv(j,:)/temp
  DO i = j+1, n
    temp = a(i,j)
    a(i,j:n) = a(i,j:n) - temp*a(j,j:n)
    ainv(i,:) = ainv(i,:) - temp*ainv(j,:)
  ENDDO
ENDDO

! Eliminate above diagonal terms
DO j = n, 2, -1
  DO i = 1, j-1
    temp = a(i,j)
    a(i,j) = a(i,j) - temp*a(j,j)   ! Should be zero if we've done it all right!
    ainv(i,:) = ainv(i,:) - temp*ainv(j,:)
  ENDDO
ENDDO


END SUBROUTINE matinv

! =======================================================================

SUBROUTINE opendumpfiles

! Initialize input output tasks:
! Open files and write headers

USE channels
USE grid
USE timestep
IMPLICIT NONE

REAL*4, ALLOCATABLE :: flongstar4(:), flatstar4(:), &
                       vlongstar4(:), vlatstar4(:)
CHARACTER*8 :: ytime
CHARACTER*16 :: yname

! -----------------------------------------------------------------------

ALLOCATE(flongstar4(nfacex), flatstar4(nfacex), &
         vlongstar4(nvertx), vlatstar4(nvertx))

! Construct timestep element of filename
WRITE(ytime,'(I8.8)') istep

! File for Matlab primal grid output
yname = 'dump1_'//ytime//'.m'
OPEN(chandumpm1,FILE=yname)

! Convert to single precision to reduce size of output
! and improve readability!
flongstar4 = flong(:,ngrids)
flatstar4 = flat(:,ngrids)
! Write header information
WRITE(chandumpm1,*)   'npt = ',nfacex,';'
WRITE(chandumpm1,*)   'nr = ',nvertx,';'
WRITE(chandumpm1,*)   'nprmx = ',nevmx,';'
WRITE(chandumpm1,*)   'rlist = [ ...'
WRITE(chandumpm1,889) fofv(:,:,ngrids)
WRITE(chandumpm1,*)   ' ];'
WRITE(chandumpm1,*)   'rlist = reshape(rlist,nr,nprmx);'
WRITE(chandumpm1,*)   'long = [ ...'
WRITE(chandumpm1,888) flongstar4
WRITE(chandumpm1,*)   ' ];'
WRITE(chandumpm1,*)   'lat = [ ...'
WRITE(chandumpm1,888) flatstar4
WRITE(chandumpm1,*)   ' ];'

! File for Matlab dual grid output
yname = 'dump2_'//ytime//'.m'
OPEN(chandumpm2,FILE=yname)

! Convert to single precision to reduce size of output
! and improve readability!
vlongstar4 = vlong(:,ngrids)
vlatstar4 = vlat(:,ngrids)
! Write header information
WRITE(chandumpm2,*)   'npt = ',nvertx,';'
WRITE(chandumpm2,*)   'nr = ',nfacex,';'
WRITE(chandumpm2,*)   'nprmx = ',nefmx,';'
WRITE(chandumpm2,*)   'rlist = [ ...'
WRITE(chandumpm2,889) voff(:,:,ngrids)
WRITE(chandumpm2,*)   ' ];'
WRITE(chandumpm2,*)   'rlist = reshape(rlist,nr,nprmx);'
WRITE(chandumpm2,*)   'long = [ ...'
WRITE(chandumpm2,888) vlongstar4
WRITE(chandumpm2,*)   ' ];'
WRITE(chandumpm2,*)   'lat = [ ...'
WRITE(chandumpm2,888) vlatstar4
WRITE(chandumpm2,*)   ' ];'

888 FORMAT(E16.4)
889 FORMAT(I8)

DEALLOCATE(flongstar4, flatstar4, vlongstar4, vlatstar4)

! -------------------------------------------------------------------

END SUBROUTINE opendumpfiles

! ===================================================================

SUBROUTINE closedumpfiles

! Close Matlab output files

USE channels
IMPLICIT NONE

! -------------------------------------------------------------------

CLOSE(chandumpm1)
CLOSE(chandumpm2)

! -------------------------------------------------------------------

END SUBROUTINE closedumpfiles

! ===================================================================

SUBROUTINE dumpm(q,ytitle,npt,ygrid)

! Dump a 2D field on an unstructured grid ready for reading and
! contouring by Matlab

! q(npt)                 The data to be output
! ytitle                 Title of data
! ygrid                  'primal' or 'dual  '

USE channels
IMPLICIT NONE

INTEGER, INTENT(IN) :: npt
REAL*8, INTENT(IN) :: q(npt)
CHARACTER*(*), INTENT(IN) :: ytitle
CHARACTER*6, INTENT(IN) :: ygrid
REAL*4 :: qstar4(npt)
INTEGER chanm

! -----------------------------------------------------------------------

! Convert to single precision to reduce size of output
! and improve readability! (Shouldn't matter with formatted output)
qstar4 = q

IF (ygrid == 'primal') THEN
  chanm = chandumpm1
ELSEIF (ygrid == 'dual  ') THEN
  chanm = chandumpm2
ELSE
  PRINT *,'Invalid option: ygrid = ',ygrid,' in subroutine dumpm'
  STOP
ENDIF

! Output
WRITE(chanm,*)   'ytitle = ''',ytitle,''';'
WRITE(chanm,*)   'q = [ ...'
WRITE(chanm,888) qstar4
WRITE(chanm,*)   ' ];'
WRITE(chanm,*)   'jtcontour'

888 FORMAT(E16.4)

! -----------------------------------------------------------------------

END SUBROUTINE dumpm

! =======================================================================

SUBROUTINE output

! Standard output fields

USE constants
USE grid
USE state
USE errdiag
USE timestep
USE work
IMPLICIT NONE

INTEGER :: nf, ne, nv
REAL*8 :: phiv2(nvertx)

! -----------------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)


! Open output files and write header information
CALL opendumpfiles

! Surface height
phi0 = (phi2 + orog2)/(gravity*farea(:,ngrids))
CALL dumpm(phi0,'h',nf,'primal')
! Surface geopotential
!phi0 = (phi2 + orog2)/(farea(:,ngrids))
!CALL dumpm(phi0,'phi',nf,'primal')
!phierr = (phi2 - phi2_init)/farea(:,ngrids)
!CALL dumpm(phierr,'phierr',nf,'primal')

! Vorticity and potential vorticity
CALL Ddual2(v1,zeta2,ngrids,ne,nv)
CALL operR(phi2,phiv2,ngrids,nf,nv)
pv = (zeta2 + planvort2)/phiv2
zeta2 = zeta2*jstar(:,1,ngrids)
CALL dumpm(zeta2,'Vorticity',nv,'dual  ')
! CALL dumpm(pv,'PV',nv,'dual  ')

! Close output files
CALL closedumpfiles


! -----------------------------------------------------------------------

END SUBROUTINE output

! =======================================================================

SUBROUTINE diffref

! Compute and output the surface height difference from a reference solution

USE state
USE timestep
USE constants
USE channels
IMPLICIT NONE

INTEGER, PARAMETER :: nreftime = 5
INTEGER :: reftime(nreftime), ilist, if0
REAL*8 :: h(nfacex), href(nfacex), l1, l2, linf
CHARACTER*12 :: yrefpre
CHARACTER*10 :: yres, ytime
CHARACTER*17 :: yname

! -----------------------------------------------------------------------

! Reference solution file prefixes
IF (nefmx == 4) THEN
  yrefpre = 'TC5ref_cube_'
ELSE
  yrefpre = 'TC5ref_hex__'
ENDIF

! Resolution (for building file name)
WRITE(yres,'(I10.10)') nface(ngrids)

! List of times at which reference solution is required
reftime(1) = NINT(3600.0d0)              ! 1 hour
reftime(2) = NINT(86400.0d0)             ! 1 day
reftime(3) = NINT(5.0d0*86400.0)         ! 5 days
reftime(4) = NINT(10.0D0*86400.0)        ! 10 days
reftime(5) = NINT(15.0d0*86400.0)        ! 15 days

! Surface height
h = (phi2 + orog2)/(gravity*farea(:,ngrids))

DO ilist = 1, nreftime

  IF (NINT(istep*dt) == reftime(ilist)) THEN

    ! Reference solution is required at this time
    WRITE(ytime,'(I10.10)') reftime(ilist)
    ! Read reference solution
    OPEN(chanrefin,FILE=yrefpre//yres//'_'//ytime//'.dat',FORM='UNFORMATTED')
    DO if0 = 1, nface(ngrids)
      READ(chanrefin) href(if0)
    ENDDO
    CLOSE(chanrefin)

    ! File for Matlab primal grid output
    yname = 'err1_'//ytime//'.m'
    OPEN(chanerrout,FILE=yname)

    ! Write header information
    WRITE(chanerrout,*)   'npt = ',nfacex,';'
    WRITE(chanerrout,*)   'nr = ',nvertx,';'
    WRITE(chanerrout,*)   'nprmx = ',nevmx,';'
    WRITE(chanerrout,*)   'rlist = [ ...'
    WRITE(chanerrout,889) fofv(:,:,ngrids)
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'rlist = reshape(rlist,nr,nprmx);'
    WRITE(chanerrout,*)   'long = [ ...'
    WRITE(chanerrout,888) flong(:,ngrids)
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'lat = [ ...'
    WRITE(chanerrout,888) flat(:,ngrids)
    WRITE(chanerrout,*)   ' ];'

    ! Output
    WRITE(chanerrout,*)   'ytitle = ''Reference h time ',ytime,''';'
    WRITE(chanerrout,*)   'q = [ ...'
    WRITE(chanerrout,888) href
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'figure(1)'
    WRITE(chanerrout,*)   'jtcontour'
    WRITE(chanerrout,*)   'ytitle = ''Solution h time ',ytime,''';'
    WRITE(chanerrout,*)   'q = [ ...'
    WRITE(chanerrout,888) h
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'figure(2)'
    WRITE(chanerrout,*)   'jtcontour'
    h = h - href
    WRITE(chanerrout,*)   'ytitle = ''h error time ',ytime,''';'
    WRITE(chanerrout,*)   'q = [ ...'
    WRITE(chanerrout,888) h
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'figure(3)'
    WRITE(chanerrout,*)   'jtcontour'

    CLOSE(chanerrout)

    ! Error norms
    l1 = SUM(ABS(h)*farea(:,ngrids))/SUM(farea(:,ngrids))
    l2 = SQRT((SUM(h*h*farea(:,ngrids)))/SUM(farea(:,ngrids)))
    linf = MAXVAL(ABS(h))
    print *,'h error: L1   = ',l1 ,'  L2   = ',l2 ,'  Linf = ',linf

  ENDIF

ENDDO

888 FORMAT(E16.4)
889 FORMAT(I8)


! -----------------------------------------------------------------------

END SUBROUTINE diffref

! =======================================================================

SUBROUTINE diagnostics

! To compute and output some basic diagnostics

USE channels
USE constants
USE timestep
USE state
USE work
USE errdiag
IMPLICIT NONE

INTEGER :: nf, ne, nv, ie0, if1, if2
REAL*8 :: l1, l2, linf, meanphis, ape, ke, pz, temp, mass, &
          phibarerr, xpverr, l1v, l2v, linfv

! -----------------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)


! For test case 2 or other stready test cases
!phierr = (phi2 - phi2_init)/farea(:,ngrids)
!l1 = SUM(ABS(phierr)*farea(:,ngrids))/SUM(farea(:,ngrids))
!l2 = SQRT((SUM(phierr*phierr*farea(:,ngrids)))/SUM(farea(:,ngrids)))
!linf = MAXVAL(ABS(phierr))
!v1err = (v1 - v1_init)/ddist(:,ngrids)
!l1v = SUM(ABS(v1err)*ldist(:,ngrids)*ddist(:,ngrids)) &
!     /SUM(ldist(:,ngrids)*ddist(:,ngrids))
!l2v = SQRT((SUM(v1err*v1err*ldist(:,ngrids)*ddist(:,ngrids))) &
!     /SUM(ldist(:,ngrids)*ddist(:,ngrids)))
!linfv = MAXVAL(ABS(v1err))
!print *,'Phierr: L1   = ',l1 ,'  L2   = ',l2 ,'  Linf = ',linf
!print *,'Verr:   L1   = ',l1v,'  L2   = ',l2v,'  Linf = ',linfv
!! Write to file
!WRITE(chanerr,*) istep, l1, l2, linf


! Compute dual grid geopotential and compare with tracer
CALL operR(phi2,phibar2,ngrids,nf,nv)
phibarerr = MAXVAL(ABS((phibar2 - xphibar2)*jstar(:,1,ngrids))) &
           /MAXVAL(ABS(phibar2*jstar(:,1,ngrids)))
! print *,'Max difference between dual phi and tracer: ',phibarerr


! Compute PV and compare with tracer
! Relative vorticity
CALL Ddual2(v1,zeta2,ngrids,ne,nv)
! Add planetary contribution
zeta2 = zeta2 + planvort2
xpverr = MAXVAL(ABS((zeta2 - xzeta2)/phibar2)) &
       /MAXVAL(ABS(zeta2/phibar2))
!print *,'Max difference between dual PV and tracer: ', xpverr


! Total mass
mass = SUM(phi2)
!print *,'Total mass = ',mass


! Available energy
! First find mean surface geopotential
meanphis = SUM(phi2 + orog2)/SUM(farea(:,ngrids))
! Available potential energy
b0_temp = phi2 + orog2 - meanphis*farea(:,ngrids)
CALL HodgeI(b0_temp,b0,ngrids,nf)
ape = 0.5d0*SUM(b0*b0_temp)
CALL findke(v1,b0,ne,nf)
ke = SUM(b0*phi2)
!print *,'APE = ',ape,'    KE = ',ke,'    Total = ',ape + ke


! Potential enstrophy
pv = zeta2/phibar2
pz = 0.5d0*SUM(zeta2*pv)
! print *,'Potential enstrophy: ',pz

! Write to file
WRITE(chandiag,888) istep, phibarerr, xpverr, mass, ape, ke, pz
888 FORMAT(I8,6E20.12)

! -----------------------------------------------------------------------

END SUBROUTINE diagnostics

! =======================================================================

SUBROUTINE dumpgrid

! Output the grid coordinates in a simple format for use in generating
! reference solutions

USE runtype
USE channels
USE grid

IMPLICIT NONE
INTEGER :: if0
REAL*8 :: long, lat

! -----------------------------------------------------------------------

OPEN(chanrefgrd,FILE=ygridcoords,FORM='UNFORMATTED')

WRITE(chanrefgrd) nface(ngrids)
DO if0 = 1, nface(ngrids)
  ! Dual vertices
  ! long = flong(if0,ngrids)
  ! lat = flat(if0,ngrids)
  ! Or cell centroids
  CALL centroid(if0,long,lat)
  WRITE(chanrefgrd) long, lat
ENDDO

CLOSE(chanrefgrd)

! -----------------------------------------------------------------------

END SUBROUTINE dumpgrid

! =======================================================================

SUBROUTINE addcompmode

! Construct a grid scale pattern in a stream function
! and use this to perturb the velocity field.

USE state
IMPLICIT NONE
INTEGER :: nf, ne, nv, iter
REAL*8 :: psi(nvertx), du(nedgex), dv(nedgex), dz(nvertx), &
          dpsi(nvertx), amp

! -----------------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)

! Initialize to zero
psi = 0.0d0

! Plant seeds:
! Vertex near (pi,0)
psi(13110) = 1.0d0
! Vertex near (0,pi)
psi(16593) = 1.0d0

! Grow noise by antidiffusion
DO iter = 1, 20
  CALL Dprimal1(psi,du,ngrids,nv,ne)
  CALL HodgeHinv(du,dv,ngrids,ne)
  CALL Ddual2(dv,dz,ngrids,ne,nv)
  dz = -dz
  CALL HodgeJ(dz,dpsi,ngrids,nv)
  print *,'maxval dpsi = ',maxval(dpsi)
  psi = psi - 2.0e9*dpsi
  psi = MIN(psi,1.0d0)
  psi = MAX(psi,-1.0d0)
ENDDO

! Now compute velocity increments
CALL Dprimal1(psi,du,ngrids,nv,ne)
CALL HodgeHinv(du,dv,ngrids,ne)

! Normalize to desired amplitude
amp = MAXVAL(ABS(dv)/ddist(:,ngrids))
dv = dv/amp

! And add to velocity field
v1 = v1 + dv


! -----------------------------------------------------------------------

END SUBROUTINE addcompmode

! =======================================================================

SUBROUTINE testop

! To test various exterior derivative and Hodge star operators

USE grid
USE constants
IMPLICIT NONE

INTEGER :: igrid, nf, ne, nv, if1, iv1
REAL*8, ALLOCATABLE :: ff1(:), ff2(:), ff3(:), &
                       fe1(:), fe2(:), fe3(:), &
                       fv1(:), fv2(:), fv3(:)
REAL*8 :: long, lat

! Which properties to check ?
LOGICAL :: ldata = .false.        ! Test data on primal and dual points
LOGICAL :: ldd = .false.          ! DD f = 0 on primal and dual grids
LOGICAL :: linv = .false.         ! Inverses of I, J and H operators
LOGICAL :: llaplacian = .false.    ! Accuracy of scalar and vector Laplacians
LOGICAL :: ltrisk = .true.       ! TRiSK-related properties

! -------------------------------------------------------------------

! Which grid shall we test?
!igrid = 1
DO igrid = 1, ngrids
! DO igrid = 2, 2

nf = nface(igrid)
ne = nedge(igrid)
nv = nvert(igrid)

print *,' '
print *,'--------------------------'
print *,' '
print *,'GRID ',igrid,'  nf = ',nf,'  ne = ',ne,'  nv = ',nv,'  dof = ',nf + ne
print *,' '

print *,'Maximum ddist = ',MAXVAL(ddist(1:ne,igrid))
print *,'max l / min l = ',MAXVAL(ldist(1:ne,igrid))/MINVAL(ldist(1:ne,igrid)), &
      '  max d / min d = ',MAXVAL(ddist(1:ne,igrid))/MINVAL(ddist(1:ne,igrid)), &
      '  max A / min A = ',MAXVAL(farea(1:nf,igrid))/MINVAL(farea(1:nf,igrid))


ALLOCATE(ff1(nf),ff2(nf),ff3(nf))
ALLOCATE(fe1(ne),fe2(ne),fe3(ne))
ALLOCATE(fv1(nv),fv2(nv),fv3(nv))

! Set up test data
DO if1 = 1, nf
  long = flong(if1,igrid)
  lat = flat(if1,igrid)
  ! ff1(if1) = SIN(lat)
  ff1(if1) = COS(lat)*SIN(long)
ENDDO

DO iv1 = 1, nv
  long = vlong(iv1,igrid)
  lat = vlat(iv1,igrid)
  ! fv1(iv1) = SIN(lat)
  fv1(iv1) = COS(lat)*SIN(long)
ENDDO

IF (ldata) THEN

  ! Check ff1 and fv1 are consistent
  CALL HodgeIinv(ff1,ff2,igrid,nf)
  CALL operR(ff2,fv2,igrid,nf,nv)
  CALL HodgeJ(fv2,fv3,igrid,nv)
  print *,'Check dual grid data consistent with primal'
  print *,'fv1 = ',fv1(1:20)
  print *,'fv3 = ',fv3(1:20)
  print *,'fv3 should approx equal fv1 '
  print *,'Biggest diff ',MAXVAL(ABS(fv3 - fv1))
  print *,'  '

ENDIF

IF (ldd) THEN

  ! Check DD = 0
  CALL Ddual1(ff1,fe1,igrid,nf,ne)
  CALL Ddual2(fe1,fv2,igrid,ne,nv)
  print *,'ff1 = ',ff1(1:40)
  print *,'DD ff1 = ',fv2(1:40)
  print *,'DD ff1 should exactly equal zero'
  print *,' '

  CALL Dprimal1(fv1,fe1,igrid,nv,ne)
  CALL Dprimal2(fe1,ff2,igrid,ne,nf)
  print *,'fv1 = ',fv1(1:40)
  print *,'DD fv1 = ',ff2(1:40)
  print *,'DD fv1 should exactly equal zero'
  print *,' '

ENDIF

IF (linv) THEN

  ! Check inverse of I
  CALL HodgeIinv(ff1,ff2,igrid,nf)
  CALL HodgeI(ff2,ff3,igrid,nf)
  print *,'Inverse of I'
  print *,'ff1 = ',ff1(1:40)
  print *,'ff3 = ',ff3(1:40)
  print *,'ff3 should (almost) exactly equal ff1'
  print *,' '

  ! Check inverse of J
  CALL HodgeJinv(fv1,fv2,igrid,nv)
  CALL HodgeJ(fv2,fv3,igrid,nv)
  print *,'Inverse of J'
  print *,'fv1 = ',fv1(1:40)
  print *,'fv3 = ',fv3(1:40)
  print *,'fv3 should (almost) exactly equal fv1'
  print *,' '

  ! Check inverse of H
  CALL Ddual1(ff1,fe1,igrid,nf,ne)
  CALL HodgeH(fe1,fe2,igrid,ne)
  CALL HodgeHinv(fe2,fe3,igrid,ne)
  print *,'Inverse of H'
  print *,'fe1 = ',fe1(1:40)
  print *,'fe3 = ',fe3(1:40)
  print *,'fe3 should (almost) exactly equal fe1'
  print *,' '

ENDIF

IF (llaplacian) THEN

  ! Check Laplacian of scalar (primal cells)
  CALL Ddual1(ff1,fe1,igrid,nf,ne)
  CALL HodgeH(fe1,fe2,igrid,ne)
  CALL Dprimal2(fe2,ff2,igrid,ne,nf)
  CALL HodgeI(ff2,ff3,igrid,nf)
  ff2 = -2.0d0*ff1
  ff3 = ff3*rearth*rearth
  print *,'Laplacian of scalar (primal)'
!  print *,'ff2 = ',ff2(1:40)
!  print *,'ff3 = ',ff3(1:40)
!  print *,'ff3 should approximately equal ff2'
  ff2 = ff3 - ff2
!  print *,'Err = ',ff2(1:40)
  print *,'Max error = ',MAXVAL(ABS(ff2))
  print *,'RMS error = ',SQRT((SUM(ff2*ff2))/nf) 
  print *,' '

  ! Check Laplacian of a scalar (dual cells)
  CALL Dprimal1(fv1,fe1,igrid,nv,ne)
  CALL HodgeHinv(fe1,fe2,igrid,ne)
  CALL Ddual2(fe2,fv2,igrid,ne,nv)
  CALL HodgeJ(fv2,fv3,igrid,nv)
  fv2 = -2.0d0*fv1
  fv3 = -fv3*rearth*rearth
  print *,'Laplacian of scalar (dual)'
!  print *,'fv2 = ',fv2(1:40)
!  print *,'fv3 = ',fv3(1:40)
!  print *,'fv3 should approximately equal fv2'
  fv2 = fv3 - fv2
!  print *,'Err = ',fv2(1:40)
  print *,'Max error = ',MAXVAL(ABS(fv2))
  print *,'RMS error = ',SQRT((SUM(fv2*fv2))/nv) 
  print *,' '

  ! Check Laplacian of vector
  CALL Ddual1(ff1,fe1,igrid,nf,ne)
  CALL Dprimal1(fv1,fe2,igrid,nv,ne)
  CALL HodgeHinv(fe2,fe3,igrid,ne)
  fe1 = fe1 - fe3
  CALL graddiv(fe1,fe2,igrid,ne,nf)
  CALL curlcurl(fe1,fe3,igrid,ne,nv)
  fe2 = fe2 + fe3
  fe3 = -2.0*fe1
  fe2 = fe2*rearth*rearth
  print *,'Laplacian of vector'
!  print *,'fe3 = ',fe3(1:40)
!  print *,'fe2 = ',fe2(1:40)
!  print *,'fe2 should approximately equal fe3'
  fe2 = fe2 - fe3
!  print *,'Err = ',fe2(1:40)
  print *,'Max error = ',MAXVAL(ABS(fe2))
  print *,'RMS error = ',SQRT((SUM(fe2*fe2))/ne) 
  print *,' '

ENDIF

IF (ltrisk) THEN

  ! Check TRiSK property
  CALL Dprimal1(fv1,fe1,igrid,nv,ne)
  CALL operW(fe1,fe2,igrid,ne)
  CALL Ddual2(fe2,fv2,igrid,ne,nv)
  print *,'fv1 = ',fv1(1:40)
  print *,'DWD fv1 = ',fv2(1:40)
  print *,'DWD fv1 should exactly equal zero'
  print *,' '

  CALL Ddual1(ff1,fe1,igrid,nf,ne)
  CALL HodgeH(fe1,fe2,igrid,ne)
  CALL Dprimal2(fe2,ff2,igrid,ne,nf)
  CALL operR(ff2,fv2,igrid,nf,nv)
  CALL operW(fe2,fe1,igrid,ne)
  CALL Ddual2(fe1,fv3,igrid,ne,nv)
  print *,'fv2 = RDHD ff1 = ',fv2(1:40)
  print *,'fv3 = DWHD ff1 = ',fv3(1:40)
  print *,'fv3 should exactly equal -fv2'
  print *,' '

  ! Check antisymmetry of W
  CALL Dprimal1(fv1,fe1,igrid,nv,ne)
  CALL operW(fe1,fe2,igrid,ne)
  fe2 = fe1*fe2
  print *,'fe1 W fe1 = ',SUM(fe2),'   should be zero'
  print *,' '

  ! Check accuracy of W
  CALL Ddual1(ff1*rearth,fe1,igrid,nf,ne)
  CALL HodgeH(fe1,fe2,igrid,ne)
  CALL operW(fe2,fe1,igrid,ne)
  CALL HodgeH(fe1,fe2,igrid,ne)
  CALL Dprimal1(fv1*rearth,fe3,igrid,nv,ne)
  fe2 = fe2/ldist(:,igrid)
  fe3 = fe3/ldist(:,igrid)
  print *,'Uperp = HWHV = ',fe2(1:40)
  print *,'True Uperp = ',fe3(1:40)
  fe2 = fe2 - fe3
  print *,'Accuracy of W: divergent flow'
  print *,'Max error = ',MAXVAL(ABS(fe2))
  print *,'RMS error = ',SQRT(SUM(fe2*fe2)/ne)
  print *,' '
  CALL Dprimal1(fv1*rearth,fe1,igrid,nv,ne)
  CALL operW(fe1,fe2,igrid,ne)
  CALL Ddual1(ff1*rearth,fe3,igrid,nf,ne)
  fe2 = -fe2/ddist(:,igrid)
  fe3 = fe3/ddist(:,igrid)
  print *,'Vperp = WU = ',fe2(1:40)
  print *,'True Vperp = ',fe3(1:40)
  fe2 = fe2 - fe3
  print *,'Accuracy of W: rotational flow'
  print *,'Max error = ',MAXVAL(ABS(fe2))
  print *,'RMS error = ',SQRT(SUM(fe2*fe2)/ne)
  print *,' '

  ! Check properties of R
  ff1 = farea(1:nf,igrid)
  CALL operR(ff1,fv1,igrid,nf,nv)
  CALL HodgeJ(fv1,fv2,igrid,nv)
  print *,'Property of R'
  print *,'fv2 = ',fv2(1:40)
  print *,'fv2 should exactly equal 1'
  print *,' '

ENDIF

DEALLOCATE(ff1,ff2,ff3,fe1,fe2,fe3,fv1,fv2,fv3)

ENDDO

! -------------------------------------------------------------------

END SUBROUTINE testop

!     ================================================================

SUBROUTINE fsphere

! To construct system matrix for linearized equations on an f-sphere
! for subsequent calculation of eigenvalues

USE grid
IMPLICIT NONE

INTEGER :: igrid, nf, ne, nv, nmat, if1, ie1, i, j
REAL*8, ALLOCATABLE :: phi(:), v(:), u(:), vperp(:), temp(:), a(:,:), &
                       phidot(:), vdot(:)
REAL*8 :: f0 = 1.4584d-4, phi0 = 1.0d5

! -------------------------------------------------------------------

! Which grid shall we test?
igrid = 2

nf = nface(igrid)
ne = nedge(igrid)
nv = nvert(igrid)
nmat = nf + ne
print *,'igrid = ',igrid,' nmat = ',nmat

ALLOCATE(phi(nf),v(ne),u(ne),vperp(ne),temp(nf),a(nmat,nmat),phidot(nf),vdot(ne))

! Set each input variable in turn to 1
DO if1 = 1, nf
  phi = 0.0d0
  v = 0.0d0
  phi(if1) = 1.0d0

  ! Compute u and phi tendencies
  CALL HodgeH(v,u,igrid,ne)
  CALL Dprimal2(u,temp,igrid,ne,nf)
  phidot = -phi0*temp

  CALL HodgeI(phi,temp,igrid,nf)
  CALL Ddual1(temp,vdot,igrid,nf,ne)
  vdot = -vdot

  CALL perp(v,vperp,igrid,ne)
  vdot = vdot - f0*vperp

  ! And save as column of system matrix
  a(1:nf,if1) = phidot
  a(nf+1:nmat,if1) = vdot
ENDDO

! Set each input variable in turn to 1
DO ie1 = 1, ne
  phi = 0.0d0
  v = 0.0d0
  v(ie1) = 1.0d0

  ! Compute u and phi tendencies
  CALL HodgeH(v,u,igrid,ne)
  CALL Dprimal2(u,temp,igrid,ne,nf)
  phidot = -phi0*temp

  CALL HodgeI(phi,temp,igrid,nf)
  CALL Ddual1(temp,vdot,igrid,nf,ne)
  vdot = -vdot

  CALL perp(v,vperp,igrid,ne)
  vdot = vdot - f0*vperp

  ! And save as column of system matrix
  a(1:nf,nf+ie1) = phidot
  a(nf+1:nmat,nf+ie1) = vdot
ENDDO


! Write out matrix for use by Matlab
OPEN(24,file='aaa.m')
WRITE(24,*) 'nmat = ',nmat
WRITE(24,*) 'a = ['
DO j = 1, nmat
  DO i = 1, nmat
    WRITE(24,*) a(i,j)
  ENDDO
ENDDO
WRITE(24,*) '];'


DEALLOCATE(phi,v,u,vperp,temp,a)

! -------------------------------------------------------------------

END SUBROUTINE fsphere

! ====================================================================

SUBROUTINE nmodes

! To construct system matrix for linearized equations on a rotating sphere
! for subsequent calculation of eigenvalues. A solid body rotation basic
! state should be set up in initial.

USE grid
USE state
USE errdiag
IMPLICIT NONE

INTEGER :: nf, ne, nv, nmat, if1, ie1, i, j
REAL*8, ALLOCATABLE :: v1exac(:), a(:,:)
REAL*8 :: pert

! -------------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)
nmat = nf + ne
print *,'grid = ',ngrids,' nmat = ',nmat

ALLOCATE(v1exac(ne),a(nmat,nmat))


! First compute effect of time stepping the basic state (because
! it won't be perfectly steady)
phi2 = phi2_init
v1 = v1_init
CALL step
phiexac = phi2
v1exac = v1


! Perturb each input variable in turn
! First the mass field
DO if1 = 1, nf
  phi2 = phi2_init
  v1 = v1_init
  pert = farea(if1,ngrids)
  phi2(if1) = phi2(if1) + pert

  ! Time step perturbed state
  CALL step

  ! And save as column of system matrix
  a(1:nf,if1) = (phi2 - phiexac)/pert
  a(nf+1:nmat,if1) = (v1 - v1exac)/pert
ENDDO

! Now the velocity field
DO ie1 = 1, ne
  phi2 = phi2_init
  v1 = v1_init
  pert = ddist(ie1,ngrids)
  v1(ie1) = v1(ie1) + pert

  ! Time step perturbed state
  CALL step

  ! And save as column of system matrix
  a(1:nf,nf+ie1) = (phi2 - phiexac)/pert
  a(nf+1:nmat,nf+ie1) = (v1 - v1exac)/pert
ENDDO


! Write out matrix for use by Matlab
OPEN(24,file='aaa.m')
WRITE(24,*) 'nmat = ',nmat
WRITE(24,*) 'a = ['
DO j = 1, nmat
  DO i = 1, nmat
    WRITE(24,*) a(i,j)
  ENDDO
ENDDO
WRITE(24,*) '];'


DEALLOCATE(v1exac,a)

! -------------------------------------------------------------------

END SUBROUTINE nmodes

! ====================================================================

SUBROUTINE testmg

! To test the multigrid Helmholtz solver

USE grid
USE helmcoeff    ! Just for testing
IMPLICIT NONE

INTEGER :: nf, ne, nv, if1, iv1
REAL*8, ALLOCATABLE :: ff1(:), ff2(:), ff3(:)
REAL*8 :: long, lat

! -------------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)

print *,' '
print *,'--------------------------'
print *,' '
print *,'Testing mgsolve '
print *,' '

ALLOCATE(ff1(nf),ff2(nf),ff3(nf))


! Build coefficients used in Helmholtz operator on all grids
phiref = 1.0d5*farea(:,ngrids)
CALL buildhelm


! Set up test data
! Large-scale part
DO if1 = 1, nf
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  ! ff2(if1) = SIN(lat)
  ff2(if1) = COS(lat)*SIN(long)
ENDDO
! Convert to area integrals
CALL HodgeIinv(ff2,ff1,ngrids,nfacex)
! Plus small-scale part
ff1(10) = 2.0d0*ff1(10)
print *,'Original field ff1 =     ',ff1(1:40)
print *,' '

CALL helmholtz(ff1,ff2,ngrids,nf,ne)
print *,'ff2 = Helm(ff1) =        ',ff2(1:40)
print *,' '

CALL mgsolve(ff3,ff2,ngrids)
print *,'Soln of Helmholtz ff3 = ', ff3(1:40)
print *,' '


END SUBROUTINE testmg

! ====================================================================

SUBROUTINE testadv

USE channels
USE constants
USE state
USE work
USE errdiag
USE timestep
USE advection
IMPLICIT NONE

INTEGER :: if0, iv0, nf, ne, nv
REAL*8 :: psi(nvertx), long, lat, flx1(nedgex), dphi2(nfacex), &
          xc, yc, zc, x0, y0, z0, u00, r, r0, h0, alpha, ca, sa, &
          clat, slat, slon, phibar(nvertx), dzeta2(nvertx), &
          lonrot, l1, l2, linf
!real*8 :: div2(nfacex), zeta2(nvertx) ! Now in module work
real*8 :: q, qmax
integer :: ivmx, ifmx

! --------------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)

! Constant geopotential
!phi2 = 1.0d0*farea(:,ngrids)

! Dual grid geopotential
!CALL operR(phi2,xphibar2,ngrids,nf,nv)
!CALL HodgeJ(xphibar2,phibar,ngrids,nv)

! Cosine bell: Williamson et al 1992
long = 1.0d0*pi
lat = 0.0d0
CALL ll2xyz(long,lat,xc,yc,zc)
r0 = 1.0d0/3.0d0
h0 = 1000.0d0
qmax = 0.0d0

DO if0 = 1, nf
  long = flong(if0,ngrids)
  lat = flat(if0,ngrids)
  CALL ll2xyz(long,lat,x0,y0,z0)
  CALL spdist(x0,y0,z0,xc,yc,zc,r)
  r = r/r0
  IF (r < 1.0d0) THEN
    q = 0.5d0*h0*(1.0d0 + COS(pi*r))
  ELSE
    q = 0.0d0
  ENDIF
  phiexac(if0) = q
  phi2(if0) = q*farea(if0,ngrids)
  if (q > qmax) then
    qmax = q
    ifmx = if0
  endif
ENDDO
print *,'Max phi = ',qmax,' in cell ',ifmx

! Cosine bell on dual grid
!DO iv0 = 1, nv
!  long = vlong(iv0,ngrids)
!  lat = vlat(iv0,ngrids)
!  CALL ll2xyz(long,lat,x0,y0,z0)
!  CALL spdist(x0,y0,z0,xc,yc,zc,r)
!  r = r/r0
!  IF (r < 1.0d0) THEN
!    q = 0.5d0*h0*(1.0d0 + COS(pi*r))
!  ELSE
!    q = 0.0d0
!  ENDIF
!  pvexac(iv0) = q
!  xzeta2(iv0) = q
!  if (q > qmax) then
!    qmax = q
!    ivmx = iv0
!  endif
!ENDDO
!print *,'Max pv = ',qmax,' in dual cell ',ivmx

! Construct absolute vorticity
!zeta2 = xzeta2*xphibar2

!print *,'max xphibar2 = ',MAXVAL(xphibar2)
!print *,'Max zeta2 = ',MAXVAL(zeta2)

dphi2 = phi2/farea(:,ngrids)
print *,'min and max phi ',minval(dphi2), maxval(dphi2)
CALL output

! Stream function: solid body rotation at angle alpha
u00 = 2.0d0*pi*rearth/(12.0d0*86400.0d0)
!alpha = 0.0d0
alpha = 0.25d0*pi
!alpha = 0.5*pi
ca = COS(alpha)
sa = SIN(alpha)
DO iv0 = 1, nv
  clat = COS(vlat(iv0,ngrids))
  slat = SIN(vlat(iv0,ngrids))
  slon = SIN(vlong(iv0,ngrids))
  psi(iv0) = u00*rearth*(ca*slat + sa*clat*slon)
ENDDO

! Non-divergent velocity field
! U is - D_1 (psi); sign of psi taken care of above
CALL Dprimal1(psi,ubar1,ngrids,nv,ne)
! Corresponding v field
CALL HodgeHinv(ubar1,vbar1,ngrids,ne)
! Perpendicular components
CALL operW(ubar1,vperpbar1,ngrids,ne)
CALL HodgeH(vperpbar1,uperpbar1,ngrids,ne)

! Compute divergence - needed to modify swept areas in
! routine primaladvflx
! (Should really use old velocity, but since div = 0 here
! it shouldn't matter.)
CALL Dprimal2(ubar1,div2,ngrids,ne,nf)
CALL HodgeI(div2,divfac,ngrids,nf)
divfac = 1.0d0/(1.0d0 + beta_v*dt*divfac)


! Primal grid mass flux
! CALL primaladvflx(phi2,mf1)

! Dual grid mass flux
! CALL operW(mf1,mfperp1,ngrids,ne)


! Loop over time steps
DO istep = 1, 12*48

  ! Compute advective fluxes
  CALL primaladvflx(phi2,flx1)

  ! Divergence of flux
  CALL Dprimal2(flx1,dphi2,ngrids,ne,nf)

  ! Update phi2
  phi2 = phi2 - dphi2

  dphi2 = phi2/farea(:,ngrids)
  print *,'Step ',istep,'  Range of phi ',MINVAL(dphi2), MAXVAL(dphi2)
  print *,'Total mass = ',SUM(phi2)

  ! Construct area integral of PV
  !pv = zeta2/phibar

  ! Compute advective fluxes
  !CALL dualadvflx(pv,flx1)

  ! Minus the divergence of flux
  !CALL Ddual2(flx1,dzeta2,ngrids,ne,nv)

  ! Update zeta2
  !zeta2 = zeta2 + dzeta2

  !pv = zeta2/xphibar2
  !print *,'Step ',istep,'  Range of pv ',MINVAL(pv), MAXVAL(pv)
  !print *,'Total pv = ',SUM(zeta2)

  time = time + dt

  ! Compute errors

  ! Centre of bell (initial bell should be at (pi,0))
  lonrot = pi + u00*time/rearth
  xc =  COS(lonrot)
  yc =  SIN(lonrot)*ca
  zc = -SIN(lonrot)*sa
  CALL xyz2ll(xc,yc,zc,long,lat)

  DO if0 = 1, nf
    long = flong(if0,ngrids)
    lat = flat(if0,ngrids)
    CALL ll2xyz(long,lat,x0,y0,z0)
    CALL spdist(x0,y0,z0,xc,yc,zc,r)
    r = r/r0
    IF (r < 1.0d0) THEN
      q = 0.5d0*h0*(1.0d0 + COS(pi*r))
    ELSE
      q = 0.0d0
    ENDIF
    phiexac(if0) = q
    if (q > qmax) then
      qmax = q
      ifmx = if0
    endif
  ENDDO

  ! Cosine bell on dual grid
  !DO iv0 = 1, nv
  !  long = vlong(iv0,ngrids)
  !  lat = vlat(iv0,ngrids)
  !  CALL ll2xyz(long,lat,x0,y0,z0)
  !  CALL spdist(x0,y0,z0,xc,yc,zc,r)
  !  r = r/r0
  !  IF (r < 1.0d0) THEN
  !    q = 0.5d0*h0*(1.0d0 + COS(pi*r))
  !  ELSE
  !    q = 0.0d0
  !  ENDIF
  !  pvexac(iv0) = q
  !  if (q > qmax) then
  !    qmax = q
  !    ivmx = iv0
  !  endif
  !ENDDO

  phierr = dphi2 - phiexac
  l1 = SUM(ABS(phierr)*farea(:,ngrids))/SUM(farea(:,ngrids))
  l2 = SQRT((SUM(phierr*phierr*farea(:,ngrids)))/SUM(farea(:,ngrids)))
  linf = MAXVAL(ABS(phierr))
  print *,'L1   = ',l1
  print *,'L2   = ',l2
  print *,'Linf = ',linf
  WRITE(chanerr,*) istep,l1,l2,linf

  !pverr = pv - pvexac
  !print *,'L1   = ',SUM(ABS(pverr)/jstar(:,1,ngrids))/SUM(1.0d0/jstar(:,1,ngrids))
  !print *,'L2   = ',SQRT((SUM(pverr*pverr/jstar(:,1,ngrids)))/SUM(1.0d0/jstar(:,1,ngrids)))
  !print *,'Linf = ',MAXVAL(ABS(pverr))

  if (modulo(istep,144) == 0) then
    CALL output
  endif

ENDDO

! Compute errors
print *,'L1   = ',SUM(ABS(phierr)*farea(:,ngrids))/SUM(farea(:,ngrids))
print *,'L2   = ',SQRT((SUM(phierr*phierr*farea(:,ngrids)))/SUM(farea(:,ngrids)))
print *,'Linf = ',MAXVAL(ABS(phierr))

!pverr = pv - pvexac
!print *,'L1   = ',SUM(ABS(pverr)/jstar(:,1,ngrids))/SUM(1.0d0/jstar(:,1,ngrids))
!print *,'L2   = ',SQRT((SUM(pverr*pverr/jstar(:,1,ngrids)))/SUM(1.0d0/jstar(:,1,ngrids)))
!print *,'Linf = ',MAXVAL(ABS(pverr))



END SUBROUTINE testadv

! ====================================================================
