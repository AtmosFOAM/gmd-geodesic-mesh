      PROGRAM GENGRID
C
C     Program to generate a nested set of hexagonal-
C     icosahedral grids on the sphere, including cross-
C     reference tables of adjacent faces, edges and vertices,
C     coordinates of faces and vertices, and lengths of edges
C     and areas of faces.
C
C     John Thuburn 1/November/1994
C
C     Modified to write GRIDMAP file unformatted
C
C     John Thuburn 2/May/1997
C
C     Adapted to optimise on both Heikes-Randall criterion
C     and centroidal criterion
C
C     John Thuburn October 2011
C
C     ---------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
C
      INTEGER NFACE(NGRIDS),NEDGE(NGRIDS),NVERT(NGRIDS)
      INTEGER IF1,IF2,IE1,IV1,IF3,NF2,IGRID,IEDGE,IV2,IV0,
     :        IGRIDM,IE2,IE3,IF11,IF12,IF21,IF22,IF31,IF32,
     :        IVN,IEN,IFN1,IFN2,IFN3,IX,IE0,IF0,IFN,IE11,IE12,
     :        IE21,IE22,IE31,IE32,NE,NF,IFAC,ITER,EMAX,FMAX,
     :        IX1,IX2,IV11,IV12,IV21,IV22
      INTEGER NEOFF(NFACEX,NGRIDS), NEOFV(NVERTX,NGRIDS)
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS)
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
      REAL*8 WEIGHT1,WEIGHT2
      REAL*8 PI,ZZ,PHI0,LONG,LAT,X,Y,Z,PX(3),PY(3),PZ(3),
     :       D1X,D1Y,D1Z,D2X,D2Y,D2Z,MAG,VX,VY,VZ,
     :       X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0,ATOT,AFACE,DA,S12,S22,S32,
     :       AMN,AMX,SFAC,EMN,EMX,DMN,DMX,CS,SN,GMN,GMX,X3,Y3,Z3,
     :       W1,W2,W3,TOT,DAV,COST,TOTCST,FTOL,DD,MAXDD,DISP,
     :       MXDISP,BDISP,S
      CHARACTER*26 YNAME
      LOGICAL LSWAP,lprint
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
      COMMON /COMWGT/ WEIGHT1,WEIGHT2
C
C     ---------------------------------------------------------------
C
C     Initialize all tables to zero
      DO 10 IGRID=1,NGRIDS
        DO 20 IF1=1,NFACEX
          DO 30 IF2=1,6
            FNXTF(IF1,IF2,IGRID)=0
            EOFF(IF1,IF2,IGRID)=0
            VOFF(IF1,IF2,IGRID)=0
   30     CONTINUE
          NEOFF(IF1,IGRID)=0
   20   CONTINUE
        DO 40 IE1=1,NEDGEX
          DO 50 IF2=1,2
            FNXTE(IE1,IF2,IGRID)=0
            FEOFE(IE1,IF2,IGRID)=0
            VOFE(IE1,IF2,IGRID)=0
   50     CONTINUE
   40   CONTINUE
        DO 60 IV1=1,NVERTX
          DO 70 IF2=1,3
            FOFV(IV1,IF2,IGRID)=0
            EOFV(IV1,IF2,IGRID)=0
   70     CONTINUE
          NEOFV(IV1,IGRID)=0
   60   CONTINUE
   10 CONTINUE
C
C     ----------------------------------------------------------------
C
C     Cross-reference tables for dodecahedron
C
      IGRID=1
C
      OPEN(80,FILE='dodecahedron.xref')
C
      NFACE(1)=12
      NEDGE(1)=30
      NVERT(1)=20
C
C     Faces next to each face
      READ(80,*)
      READ(80,*) ((FNXTF(IF1,IF2,1),IF2=1,5),IF1=1,12)
C
C     Edges of each face
      READ(80,*)
      READ(80,*) ((EOFF(IF1,IE1,1),IE1=1,5),IF1=1,12)
C
C     Faces next to each edge
      READ(80,*)
      READ(80,*) ((FNXTE(IE1,IF1,1),IF1=1,2),IE1=1,30)
C
C     Faces at the ends of each edge
      READ(80,*)
      READ(80,*) ((FEOFE(IE1,IF1,1),IF1=1,2),IE1=1,30)
C
C     Vertices of each edge
      READ(80,*)
      READ(80,*) ((VOFE(IE1,IV1,1),IV1=1,2),IE1=1,30)
C
C     Vertices of each face
      READ(80,*)
      READ(80,*) ((VOFF(IF1,IV1,1),IV1=1,5),IF1=1,12)
C
C     Faces around each vertex
      READ(80,*)
      READ(80,*) ((FOFV(IV1,IF1,1),IF1=1,3),IV1=1,20)
C
C     Construct table for edges around each vertex
      DO 80 IE1=1,NEDGE(1)
        IV1=VOFE(IE1,1,1)
        IV2=VOFE(IE1,2,1)
        CALL ADDTAB(EOFV(1,1,1),IV1,IE1,NVERTX,3)
        CALL ADDTAB(EOFV(1,1,1),IV2,IE1,NVERTX,3)
   80 CONTINUE
C
C     -----------------------------------------------------------------
C
C     Set up coordinates of face centres
      PI=4.0D0*ATAN(1.0D0)
      ZZ=2.0D0*COS(0.3D0*PI)
      PHI0=2.0D0*ASIN(1.0D0/ZZ)-0.5D0*PI
C
      FLONG(1,1)=0.0D0
      FLAT(1,1)=0.5D0*PI
      DO 200 IF1=1,5
        FLONG(1+IF1,1)=(IF1-1)*0.4D0*PI+0.2D0*PI
        FLAT(1+IF1,1)=PHI0
        FLONG(6+IF1,1)=(IF1-1)*0.4D0*PI
        FLAT(6+IF1,1)=-PHI0
  200 CONTINUE
      FLONG(12,1)=0.0D0
      FLAT(12,1)=-0.5D0*PI
C
C     ---------------------------------------------------------------
C
C     Areas of faces
      DO 300 IF1=1,12
        FAREA(IF1,1)=PI/3.0D0
  300 CONTINUE
C
C     ---------------------------------------------------------------
C
C     Coordinates of vertices
      DO 400 IV1=1,20
C       First find cartesian coords of surrounding faces
        DO 410 IF1=1,3
          IF2=FOFV(IV1,IF1,1)
          LONG=FLONG(IF2,1)
          LAT=FLAT(IF2,1)
          CALL LL2XYZ(LONG,LAT,X,Y,Z)
          PX(IF1)=X
          PY(IF1)=Y
          PZ(IF1)=Z
  410   CONTINUE
C       Two sides of the triangle joining the face centres
        D1X=PX(2)-PX(1)
        D1Y=PY(2)-PY(1)
        D1Z=PZ(2)-PZ(1)
        D2X=PX(3)-PX(1)
        D2Y=PY(3)-PY(1)
        D2Z=PZ(3)-PZ(1)
C       Find the middle (it's an equilateral triangle)
        VX=(PX(1)+PX(2)+PX(3))/3.0D0
        VY=(PY(1)+PY(2)+PY(3))/3.0D0
        VZ=(PZ(1)+PZ(2)+PZ(3))/3.0D0
C       Project back onto the sphere
        MAG=SQRT(VX*VX+VY*VY+VZ*VZ)
        VX=VX/MAG
        VY=VY/MAG
        VZ=VZ/MAG
C       Convert back to latitude/longitude
        CALL XYZ2LL(VX,VY,VZ,LONG,LAT)
        VLONG(IV1,1)=LONG
        VLAT(IV1,1)=LAT
  400 CONTINUE
C
C     Tabulate lengths of edges and distances between face centres
      EMN=1.0D0
      EMX=0.0D0
      DMN=2.0D0
      DMX=0.0D0
      DO 420 IE0=1,NEDGE(IGRID)
C       Vertices at ends of this edge
        IV1=VOFE(IE0,1,IGRID)
        IV2=VOFE(IE0,2,IGRID)
        LONG=VLONG(IV1,IGRID)
        LAT=VLAT(IV1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=VLONG(IV2,IGRID)
        LAT=VLAT(IV2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
        LDIST(IE0,IGRID)=S
        EMN=MIN(EMN,LDIST(IE0,IGRID))
        EMX=MAX(EMX,LDIST(IE0,IGRID))
C       Faces either side of this edge
        IF1=FNXTE(IE0,1,IGRID)
        IF2=FNXTE(IE0,2,IGRID)
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=FLONG(IF2,IGRID)
        LAT=FLAT(IF2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
        DDIST(IE0,IGRID)=S
        DMN=MIN(DMN,DDIST(IE0,IGRID))
        DMX=MAX(DMX,DDIST(IE0,IGRID))
  420 CONTINUE
C
C     ------------------------------------------------------------------
C
C     Generate the next grid in the series
C     ====================================
C
      IF (NGRIDS.EQ.1) GOTO 699
C
  700 CONTINUE
C
      IGRIDM=IGRID
      IGRID=IGRID+1
      NFACE(IGRID)=NFACE(IGRIDM)+NEDGE(IGRIDM)
      NEDGE(IGRID)=4*NEDGE(IGRIDM)
      NVERT(IGRID)=4*NVERT(IGRIDM)
C
C     Old faces keep their coordinates
      DO 602 IF1=1,NFACE(IGRIDM)
        FLONG(IF1,IGRID)=FLONG(IF1,IGRIDM)
        FLAT(IF1,IGRID)=FLAT(IF1,IGRIDM)
  602 CONTINUE
C
C     Loop over old edges
      DO 600 IE0=1,NEDGE(IGRIDM)
C
C       Each old edge generates a new face
        IFN=NFACE(IGRIDM)+IE0
C
C       This is next to the two old faces either side of the old edge
        IF1=FNXTE(IE0,1,IGRIDM)
        IF2=FNXTE(IE0,2,IGRIDM)
        CALL ADDTAB(FNXTF(1,1,IGRID),IF1,IFN,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IF2,IFN,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN,IF1,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN,IF2,NFACEX,6)
C
C       There is a new edge between IF1 and IFN and similarly
C       between IF2 and IFN
        IE1=2*IE0-1
        IE2=2*IE0
        CALL ADDTAB(FNXTE(1,1,IGRID),IE1,IF1,NEDGEX,2)
        CALL ADDTAB(FNXTE(1,1,IGRID),IE1,IFN,NEDGEX,2)
        CALL ADDTAB(EOFF(1,1,IGRID),IF1,IE1,NFACEX,6)
        CALL ADDTAB(EOFF(1,1,IGRID),IFN,IE1,NFACEX,6)
        CALL ADDTAB(FNXTE(1,1,IGRID),IE2,IF2,NEDGEX,2)
        CALL ADDTAB(FNXTE(1,1,IGRID),IE2,IFN,NEDGEX,2)
        CALL ADDTAB(EOFF(1,1,IGRID),IF2,IE2,NFACEX,6)
        CALL ADDTAB(EOFF(1,1,IGRID),IFN,IE2,NFACEX,6)
C
C       Coordinates of the new face
C         Midpoint of ends of edge
c          IV1=VOFE(IE0,1,IGRIDM)
c          IV2=VOFE(IE0,2,IGRIDM)
c          LONG=VLONG(IV1,IGRID)
c          LAT=VLAT(IV1,IGRID)
c          CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
c          LONG=VLONG(IV2,IGRID)
c          LAT=VLAT(IV2,IGRID)
c          CALL LL2XYZ(LONG,LAT,X,Y,Z)
C       Midpoint of old faces
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=FLONG(IF2,IGRID)
        LAT=FLAT(IF2,IGRID)
        CALL LL2XYZ(LONG,LAT,X,Y,Z)
C       Find the middle
        VX=0.5*(X1+X)
        VY=0.5*(Y1+Y)
        VZ=0.5*(Z1+Z)
C       Project back onto the sphere
        MAG=SQRT(VX*VX+VY*VY+VZ*VZ)
        VX=VX/MAG
        VY=VY/MAG
        VZ=VZ/MAG
C       Convert back to latitude/longitude
        CALL XYZ2LL(VX,VY,VZ,LONG,LAT)
        FLONG(IFN,IGRID)=LONG
        FLAT(IFN,IGRID)=LAT
C
  600 CONTINUE
C
C     Loop over old vertices
      DO 650 IV0=1,NVERT(IGRIDM)
C
C       Each vertex has three new faces
        IE1=EOFV(IV0,1,IGRIDM)
        IE2=EOFV(IV0,2,IGRIDM)
        IE3=EOFV(IV0,3,IGRIDM)
        IFN1=NFACE(IGRIDM)+IE1
        IFN2=NFACE(IGRIDM)+IE2
        IFN3=NFACE(IGRIDM)+IE3
        CALL ADDTAB(FOFV(1,1,IGRID),IV0,IFN1,NVERTX,3)
        CALL ADDTAB(FOFV(1,1,IGRID),IV0,IFN2,NVERTX,3)
        CALL ADDTAB(FOFV(1,1,IGRID),IV0,IFN3,NVERTX,3)
        CALL ADDTAB(VOFF(1,1,IGRID),IFN1,IV0,NFACEX,6)
        CALL ADDTAB(VOFF(1,1,IGRID),IFN2,IV0,NFACEX,6)
        CALL ADDTAB(VOFF(1,1,IGRID),IFN3,IV0,NFACEX,6)
C
C       These faces are mutual neighbours
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN1,IFN2,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN1,IFN3,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN2,IFN1,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN2,IFN3,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN3,IFN1,NFACEX,6)
        CALL ADDTAB(FNXTF(1,1,IGRID),IFN3,IFN2,NFACEX,6)
C
C       Note old faces next to new faces and corresponding edges
        IF11=FNXTF(IFN1,1,IGRID)
        IE11=EOFF(IFN1,1,IGRID)
        IF12=FNXTF(IFN1,2,IGRID)
        IE12=EOFF(IFN1,2,IGRID)
        IF21=FNXTF(IFN2,1,IGRID)
        IE21=EOFF(IFN2,1,IGRID)
        IF22=FNXTF(IFN2,2,IGRID)
        IE22=EOFF(IFN2,2,IGRID)
        IF31=FNXTF(IFN3,1,IGRID)
        IE31=EOFF(IFN3,1,IGRID)
        IF32=FNXTF(IFN3,2,IGRID)
        IE32=EOFF(IFN3,2,IGRID)
C
C       Each old vertex generates three new edges and three
C       new vertices
        DO 660 IX=1,3
C
          IF0=FOFV(IV0,IX,IGRIDM)
          IEN=2*NEDGE(IGRIDM)+3*(IV0-1)+IX
          IVN=NVERT(IGRIDM)+3*(IV0-1)+IX
          CALL ADDTAB(EOFV(1,1,IGRID),IV0,IEN,NVERTX,3)
          CALL ADDTAB(EOFV(1,1,IGRID),IVN,IEN,NVERTX,3)
          CALL ADDTAB(VOFE(1,1,IGRID),IEN,IV0,NEDGEX,2)
          CALL ADDTAB(VOFE(1,1,IGRID),IEN,IVN,NEDGEX,2)
C
C         New vertex is a vertex of the old face
          CALL ADDTAB(VOFF(1,1,IGRID),IF0,IVN,NFACEX,6)
          CALL ADDTAB(FOFV(1,1,IGRID),IVN,IF0,NVERTX,3)
C
C         If the old face is a neighbour of a new face then
C         the new vertex and edge also pertain to the new face
          IF ((IF0.EQ.IF11).OR.(IF0.EQ.IF12)) THEN
            CALL ADDTAB(VOFF(1,1,IGRID),IFN1,IVN,NFACEX,6)
            CALL ADDTAB(FOFV(1,1,IGRID),IVN,IFN1,NVERTX,3)
            CALL ADDTAB(EOFF(1,1,IGRID),IFN1,IEN,NFACEX,6)
            CALL ADDTAB(FNXTE(1,1,IGRID),IEN,IFN1,NEDGEX,2)
          ENDIF
          IF ((IF0.EQ.IF21).OR.(IF0.EQ.IF22)) THEN
            CALL ADDTAB(VOFF(1,1,IGRID),IFN2,IVN,NFACEX,6)
            CALL ADDTAB(FOFV(1,1,IGRID),IVN,IFN2,NVERTX,3)
            CALL ADDTAB(EOFF(1,1,IGRID),IFN2,IEN,NFACEX,6)
            CALL ADDTAB(FNXTE(1,1,IGRID),IEN,IFN2,NEDGEX,2)
          ENDIF
          IF ((IF0.EQ.IF31).OR.(IF0.EQ.IF32)) THEN
            CALL ADDTAB(VOFF(1,1,IGRID),IFN3,IVN,NFACEX,6)
            CALL ADDTAB(FOFV(1,1,IGRID),IVN,IFN3,NVERTX,3)
            CALL ADDTAB(EOFF(1,1,IGRID),IFN3,IEN,NFACEX,6)
            CALL ADDTAB(FNXTE(1,1,IGRID),IEN,IFN3,NEDGEX,2)
          ENDIF
C
C         If the old face is a neighbour of a new face then
C         the edge between the old and new faces pertains
C         to the new vertex
          IF (IF0.EQ.IF11) THEN
            CALL ADDTAB(EOFV(1,1,IGRID),IVN,IE11,NVERTX,3)
            CALL ADDTAB(VOFE(1,1,IGRID),IE11,IVN,NEDGEX,2)
          ENDIF
          IF (IF0.EQ.IF12) THEN
            CALL ADDTAB(EOFV(1,1,IGRID),IVN,IE12,NVERTX,3)
            CALL ADDTAB(VOFE(1,1,IGRID),IE12,IVN,NEDGEX,2)
          ENDIF
          IF (IF0.EQ.IF21) THEN
            CALL ADDTAB(EOFV(1,1,IGRID),IVN,IE21,NVERTX,3)
            CALL ADDTAB(VOFE(1,1,IGRID),IE21,IVN,NEDGEX,2)
          ENDIF
          IF (IF0.EQ.IF22) THEN
            CALL ADDTAB(EOFV(1,1,IGRID),IVN,IE22,NVERTX,3)
            CALL ADDTAB(VOFE(1,1,IGRID),IE22,IVN,NEDGEX,2)
          ENDIF
          IF (IF0.EQ.IF31) THEN
            CALL ADDTAB(EOFV(1,1,IGRID),IVN,IE31,NVERTX,3)
            CALL ADDTAB(VOFE(1,1,IGRID),IE31,IVN,NEDGEX,2)
          ENDIF
          IF (IF0.EQ.IF32) THEN
            CALL ADDTAB(EOFV(1,1,IGRID),IVN,IE32,NVERTX,3)
            CALL ADDTAB(VOFE(1,1,IGRID),IE32,IVN,NEDGEX,2)
          ENDIF
C
  660   CONTINUE
C
  650 CONTINUE
C
C     Calculate coordinates of vertices - Voronoi grid
      DO 667 IV1=1,NVERT(IGRID)
        CALL VRTCRD(IV1,IGRID)
  667 CONTINUE
C
C     Now tweak face and vertex positions a la Heikes and Randall
C     to minimise cost function related to inaccuracy in Laplacian
C     Outer loop over iterations
c      FTOL=1.0D-5
      FTOL=1.0D-1
      TOTCST=0.0D0
      MAXDD=0.0D0
      DO IE1=1,NEDGE(IGRID)
        CALL PENALT(IE1,IGRID,COST,DD)
        TOTCST=TOTCST+COST
	IF (DD.GT.MAXDD) THEN
          MAXDD=DD
	  EMAX=IE1
	ENDIF
      ENDDO
      print *,'Init. Heikes-Randall cost func ',totcst
      print *,'MAXDD=                         ',MAXDD,' Edge ',EMAX
      TOTCST=0.0D0
      MAXDD=0.0D0
      DO IF1=1,NFACE(IGRID)
        CALL PENALC(IF1,IGRID,COST,DD)
        TOTCST=TOTCST+COST
	IF (DD.GT.MAXDD) THEN
	  MAXDD=DD
	  FMAX=IF1
	ENDIF
      ENDDO
      print *,'Init. centroidal cost func     ',totcst
      print *,'MAXDD=                         ',MAXDD,' Face ',FMAX
C
C     Weights for HR (weight1) and centroidal (weight2) contributions
C     to cost function
      WEIGHT1=1.0D0
      WEIGHT2=0.0D0
C     Initial interval size for BRENT
      BDISP=2.0D0**(-IGRID)
C
      DO ITER=1,40
        lprint = .true.
C       Loop over faces of grid
        MXDISP=0.0D0
        DO IF2=13,NFACE(IGRID)
C         Minimize cost associated with face IF1
          CALL POWELL(IF2,IGRID,FTOL,BDISP,DISP,NGRIDS)
          MXDISP=MAX(MXDISP,DISP)
        ENDDO
        if (lprint) then
          print *,'Done iteration ',ITER
          TOTCST=0.0D0
          MAXDD=0.0D0
          DO IE1=1,NEDGE(IGRID)
            CALL PENALT(IE1,IGRID,COST,DD)
            TOTCST=TOTCST+COST
	    IF (DD.GT.MAXDD) THEN
              MAXDD=DD
	      EMAX=IE1
            ENDIF
          ENDDO
          print *,'Heikes-Randall cost function ',totcst
	  print *,'MAXDD=                       ',MAXDD,
     :                      ' Edge ',EMAX
          TOTCST=0.0D0
	  MAXDD=0.0D0
          DO IF1=1,NFACE(IGRID)
            CALL PENALC(IF1,IGRID,COST,DD)
            TOTCST=TOTCST+COST
	    IF (DD.GT.MAXDD) THEN
	      MAXDD=DD
	      FMAX=IF1
	    ENDIF
          ENDDO
          print *,'Centroidal cost function     ',totcst
	  print *,'MAXDD=                       ',MAXDD,
     :                      ' Face ',FMAX
          print *,'Max disp in BRENT = ',MXDISP
        ENDIF

C       BRENT interval size for next iteration is max displacement
C       from this iteration
        BDISP=MXDISP

      ENDDO
C
C     Tabulate lengths of edges and distances between face centres
C     across each edge
      EMN=5.0D0
      EMX=0.0D0
      DMN=5.0D0
      DMX=0.0D0
      DAV=0.0D0
      DO 670 IE0=1,NEDGE(IGRID)
C       Vertices at ends of this edge
        IV1=VOFE(IE0,1,IGRID)
        IV2=VOFE(IE0,2,IGRID)
        LONG=VLONG(IV1,IGRID)
        LAT=VLAT(IV1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=VLONG(IV2,IGRID)
        LAT=VLAT(IV2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
        LDIST(IE0,IGRID)=S
        EMN=MIN(EMN,LDIST(IE0,IGRID))
        EMX=MAX(EMX,LDIST(IE0,IGRID))
C       Faces either side of this edge
        IF1=FNXTE(IE0,1,IGRID)
        IF2=FNXTE(IE0,2,IGRID)
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=FLONG(IF2,IGRID)
        LAT=FLAT(IF2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
        DDIST(IE0,IGRID)=S
        DMN=MIN(DMN,DDIST(IE0,IGRID))
        DMX=MAX(DMX,DDIST(IE0,IGRID))
        DAV=DAV+DDIST(IE0,IGRID)/NEDGE(IGRID)
  670 CONTINUE
C
      PRINT *,' Done grid ',IGRID
C
      print *,'min side: ',emn,' max side: ',emx,' ratio: ',emn/emx
      print *,'min dist: ',dmn,' max dist: ',dmx,' ratio: ',dmn/dmx
      print *,'average side * rearth: ',dav*6371220.0D0
C
C
C
      IF (IGRID.LT.NGRIDS) GOTO 700
  699 CONTINUE
C
C     ------------------------------------------------------------------
C
C     Calculate areas of grid cells on each grid
C
      DO 850 IGRID=1,NGRIDS
C
      AMN=4.0D0*PI
      AMX=0.0D0
      ATOT=0.0D0
      DO 800 IF1=1,NFACE(IGRID)
        AFACE=0.0D0
C       Coordinates of centre of face
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X0,Y0,Z0)
C       Loop over edges in turn and calculate area of triangle
C       formed by the edge and the centre of the face
        IF (IF1.LE.12) THEN
          NE=5
        ELSE
          NE=6
        ENDIF
        DO 810 IE1=1,NE
          IE2=EOFF(IF1,IE1,IGRID)
          IV1=VOFE(IE2,1,IGRID)
          IV2=VOFE(IE2,2,IGRID)
          LONG=VLONG(IV1,IGRID)
          LAT=VLAT(IV1,IGRID)
          CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
          LONG=VLONG(IV2,IGRID)
          LAT=VLAT(IV2,IGRID)
          CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
          CALL STAREA2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,DA)
          AFACE=AFACE+DA
  810   CONTINUE
        FAREA(IF1,IGRID)=AFACE
        ATOT=ATOT+AFACE
  800 CONTINUE
C
      print *,'Grid ',igrid,' total area=',atot
C     Now scale up areas so total is 4*pi
C     (Not really necessary now that we use the correct spherical
C     triangle area)
      SFAC=4.0D0*PI/ATOT
      print *,'fraction of sphere: ',1.0D0/sfac
      DO 820 IF1=1,NFACE(IGRID)
        FAREA(IF1,IGRID)=FAREA(IF1,IGRID)*SFAC
        AMN=MIN(AMN,FAREA(IF1,IGRID))
        AMX=MAX(AMX,FAREA(IF1,IGRID))
  820 CONTINUE
C
      print *,'min area: ',amn,' max area: ',amx,' ratio: ',amn/amx
      print *,'average area =',5.101D14/NFACE(IGRID)
C
  850 CONTINUE
C
C     ------------------------------------------------------------------
C
C     Every vertex has three edges
      DO 902 IGRID=1,NGRIDS
        DO 901 IV0=1,NVERT(IGRID)
          NEOFV(IV0,IGRID)=3
  901   CONTINUE
  902 CONTINUE
C
C     Sort FNXTF into anticlockwise order on each grid
C     and sort EOFF to correspond to FNXTF
      DO 945 IGRID=1,NGRIDS
      DO 900 IF0=1,NFACE(IGRID)
        IF (IF0.LE.12) THEN
          NF=5
        ELSE
          NF=6
        ENDIF
        NEOFF(IF0,IGRID)=NF
C       Coordinates of face IF0
        LONG=FLONG(IF0,IGRID)
        LAT=FLAT(IF0,IGRID)
        CALL LL2XYZ(LONG,LAT,X,Y,Z)
        DO 910 IX1=1,NF-2
C         Coordinates of IX1'th neighbour
          IF1=FNXTF(IF0,IX1,IGRID)
          LONG=FLONG(IF1,IGRID)
          LAT=FLAT(IF1,IGRID)
          CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
          D1X=X1-X
          D1Y=Y1-Y
          D1Z=Z1-Z
          IX2=IX1
  920     CONTINUE
            IX2=IX2+1
C           Coordinates of IX2'th neighbour
            IF2=FNXTF(IF0,IX2,IGRID)
            LONG=FLONG(IF2,IGRID)
            LAT=FLAT(IF2,IGRID)
            CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
            D2X=X2-X
            D2Y=Y2-Y
            D2Z=Z2-Z
            CS=D1X*D2X+D1Y*D2Y+D1Z*D2Z
            IF (CS.LT.0.0D0) GOTO 920
            SN=X*(D1Y*D2Z-D1Z*D2Y)
     :        +Y*(D1Z*D2X-D1X*D2Z)
     :        +Z*(D1X*D2Y-D1Y*D2X)
            IF (SN.LT.0.0D0) GOTO 920
C         IF2 belongs in position IX1+1 so swap them
          IF3=FNXTF(IF0,IX1+1,IGRID)
          FNXTF(IF0,IX1+1,IGRID)=IF2
          FNXTF(IF0,IX2,IGRID)=IF3
  910   CONTINUE
        DO 930 IX1=1,NF
          IF1=FNXTF(IF0,IX1,IGRID)
          IX2=IX1-1
  940     CONTINUE
            IX2=IX2+1
            IE2=EOFF(IF0,IX2,IGRID)
            IF21=FNXTE(IE2,1,IGRID)
            IF22=FNXTE(IE2,2,IGRID)
            IF ((IF21+IF22).NE.(IF0+IF1)) GOTO 940
C         Edge IE2 corresponds to face IF1
          EOFF(IF0,IX2,IGRID)=EOFF(IF0,IX1,IGRID)
          EOFF(IF0,IX1,IGRID)=IE2
  930   CONTINUE
  900 CONTINUE
  945 CONTINUE
C
C     Order VOFF so that the k'th vertex is between the
C     k'th and (k+1)'th edges in EOFF
      DO 946 IGRID=1,NGRIDS
        DO 947 IF0=1,NFACE(IGRID)
          DO 948 IX1=1,NEOFF(IF0,IGRID)
            IX2=IX1+1
            IF (IX2.GT.NEOFF(IF0,IGRID)) IX2=1
            IE1=EOFF(IF0,IX1,IGRID)
            IE2=EOFF(IF0,IX2,IGRID)
            ! Find the common vertex of IE1 and IE2
            IV11=VOFE(IE1,1,IGRID)
            IV12=VOFE(IE1,2,IGRID)
            IV21=VOFE(IE2,1,IGRID)
            IV22=VOFE(IE2,2,IGRID)
            IF ((IV11.EQ.IV21).OR.(IV11.EQ.IV22)) THEN
              IV0=IV11
            ELSEIF ((IV12.EQ.IV21).OR.(IV12.EQ.IV22)) THEN
              IV0=IV12
            ELSE
              PRINT *,'Common vertex not found'
              STOP
            ENDIF
            VOFF(IF0,IX1,IGRID) = IV0
  948     CONTINUE
  947   CONTINUE
  946 CONTINUE

C     Find faces at ends of edges for output grid and
C     sort FEOFE so that FEOFE(1) -> FEOFE(2) is 90 deg anticlockwise
C     of FNXTE(1) -> FNXTE(2) and sort VOFE so that VOFE(1) -> VOFE(2)
C     is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2)
      DO 955 IGRID=1,NGRIDS
      DO 950 IE0=1,NEDGE(IGRID)
        IF1=FNXTE(IE0,1,IGRID)
        IF2=FNXTE(IE0,2,IGRID)
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X,Y,Z)
        LONG=FLONG(IF2,IGRID)
        LAT=FLAT(IF2,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        D1X=X1-X
        D1Y=Y1-Y
        D1Z=Z1-Z
        LSWAP=.FALSE.
        DO 960 IX1=1,2
          IV1=VOFE(IE0,IX1,IGRID)
          DO 970 IX2=1,3
            IF3=FOFV(IV1,IX2,IGRID)
            IF ((IF3.NE.IF1).AND.(IF3.NE.IF2)) THEN
              LONG=FLONG(IF3,IGRID)
              LAT=FLAT(IF3,IGRID)
              CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
              D2X=X2-X
              D2Y=Y2-Y
              D2Z=Z2-Z
              SN=X*(D1Y*D2Z-D1Z*D2Y)
     :          +Y*(D1Z*D2X-D1X*D2Z)
     :          +Z*(D1X*D2Y-D1Y*D2X)
              IF (SN.GT.0.0D0) THEN
                FEOFE(IE0,2,IGRID)=IF3
                IF (IX1.EQ.1) LSWAP=.TRUE.
              ELSE
                FEOFE(IE0,1,IGRID)=IF3
              ENDIF
            ENDIF
  970     CONTINUE
  960   CONTINUE
        IF (LSWAP) THEN
          IV1=VOFE(IE0,1,IGRID)
          VOFE(IE0,1,IGRID)=VOFE(IE0,2,IGRID)
          VOFE(IE0,2,IGRID)=IV1
        ENDIF
  950 CONTINUE
  955 CONTINUE
C
C
C     Find distance between faces along each edge
      IGRID=NGRIDS
      GMN=5.0D0
      GMX=0.0D0
      DO 980 IE0=1,NEDGE(NGRIDS)
C       Faces either end of this edge
        IF1=FEOFE(IE0,1,IGRID)
        IF2=FEOFE(IE0,2,IGRID)
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=FLONG(IF2,IGRID)
        LAT=FLAT(IF2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
        GDIST(IE0,IGRID)=S
        GMN=MIN(GMN,GDIST(IE0,IGRID))
        GMX=MAX(GMX,GDIST(IE0,IGRID))
  980 CONTINUE
C
      print *,'min area: ',amn,' max area: ',amx,' ratio: ',amn/amx
      print *,'min side: ',emn,' max side: ',emx,' ratio: ',emn/emx
      print *,'min dist: ',dmn,' max dist: ',dmx,' ratio: ',dmn/dmx
      print *,'min dist2: ',gmn,' max dist2: ',gmx,' ratio: ',gmn/gmx
C
C     ------------------------------------------------------------------
C
C     Calculate coefficients for interpolating stream function and
C     velocity potential to vertices.
C
      DO 1010 IGRID=1,NGRIDS
C
      DO 1000 IV1=1,NVERT(IGRID)
        IF1=FOFV(IV1,1,IGRID)
        IF2=FOFV(IV1,2,IGRID)
        IF3=FOFV(IV1,3,IGRID)
        LONG=VLONG(IV1,IGRID)
        LAT=VLAT(IV1,IGRID)
        CALL LL2XYZ(LONG,LAT,X,Y,Z)
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=FLONG(IF2,IGRID)
        LAT=FLAT(IF2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        LONG=FLONG(IF3,IGRID)
        LAT=FLAT(IF3,IGRID)
        CALL LL2XYZ(LONG,LAT,X3,Y3,Z3)
        X1=X1-X
        Y1=Y1-Y
        Z1=Z1-Z
        X2=X2-X
        Y2=Y2-Y
        Z2=Z2-Z
        X3=X3-X
        Y3=Y3-Y
        Z3=Z3-Z
        W1=X*(Y2*Z3-Z2*Y3)
     :    +Y*(Z2*X3-X2*Z3)
     :    +Z*(X2*Y3-Y2*X3)
        W2=X*(Y3*Z1-Z3*Y1)
     :    +Y*(Z3*X1-X3*Z1)
     :    +Z*(X3*Y1-Y3*X1)
        W3=X*(Y1*Z2-Z1*Y2)
     :    +Y*(Z1*X2-X1*Z2)
     :    +Z*(X1*Y2-Y1*X2)
        W1=ABS(W1)
        W2=ABS(W2)
        W3=ABS(W3)
        TOT=W1+W2+W3
        W1=W1/TOT
        W2=W2/TOT
        W3=W3/TOT
        COEFF(IV1,1,IGRID)=W1
        COEFF(IV1,2,IGRID)=W2
        COEFF(IV1,3,IGRID)=W3
 1000 CONTINUE
      DO 1005 IV1=NVERT(IGRID)+1,NVERTX
        COEFF(IV1,1,IGRID)=0.0D0
        COEFF(IV1,2,IGRID)=0.0D0
        COEFF(IV1,3,IGRID)=0.0D0
 1005 CONTINUE        
C
 1010 CONTINUE
C
C     ------------------------------------------------------------------
C
C     Create a file of cross reference tables and coordinates etc
C     for use by the model
C
      PRINT *,'Create a GRIDMAP file (0 or 1) ?'
      READ (5,*) IF0
      IF (IF0.EQ.1) THEN
C
        WRITE(YNAME,'(''gridmap_hex_'',I10.10,''.dat'')') NFACEX
        OPEN(22,FILE=YNAME,FORM='UNFORMATTED')
C
        PRINT *,'NFACE = ',NFACE
        PRINT *,'NEDGE = ',NEDGE
        PRINT *,'NVERT = ',NVERT
C
C        WRITE(22,*) 'GRIDMAP for NGRIDS=',NGRIDS
        WRITE(22) NGRIDS
        WRITE(22) NFACE
        WRITE(22) NEDGE
        WRITE(22) NVERT
C        WRITE(22,*) 'Number of edges of each face - all grids'
        WRITE(22) ((NEOFF(IF1,IGRID),
     :                 IF1=1,NFACE(IGRID)),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Number of edges of each vertex - all grids'
        WRITE(22) ((NEOFV(IV1,IGRID),
     :                 IV1=1,NVERT(IGRID)),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Faces next to each face - all grids'
        WRITE(22) (((FNXTF(IF1,IF2,IGRID),
     :                 IF1=1,NFACE(IGRID)),
     :                 IF2=1,6),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Edges of each face - all grids'
        WRITE(22) (((EOFF(IF1,IE2,IGRID),
     :                 IF1=1,NFACE(IGRID)),
     :                 IE2=1,6),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Vertices of each face - all grids'
        WRITE(22) (((VOFF(IF1,IV1,IGRID),
     :                 IF1=1,NFACE(IGRID)),
     :                 IV1=1,6),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Faces next to each edge - all grids'
        WRITE(22) (((FNXTE(IE1,IF2,IGRID),
     :                 IE1=1,NEDGE(IGRID)),
     :                 IF2=1,2),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Vertices of each edge - all grids'
        WRITE(22) (((VOFE(IE1,IV2,IGRID),
     :                 IE1=1,NEDGE(IGRID)),
     :                 IV2=1,2),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Faces around each vertex - all grids'
        WRITE(22) (((FOFV(IV1,IF2,IGRID),
     :                 IV1=1,NVERT(IGRID)),
     :                 IF2=1,3),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Edges around each vertex - all grids'
        WRITE(22) (((EOFV(IV1,IE1,IGRID),
     :                 IV1=1,NVERT(IGRID)),
     :                 IE1=1,3),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Coefficients for interpolation - all grids'
C        WRITE(22) (((COEFF(IV1,IF2,IGRID),
C     :                 IV1=1,NVERT(IGRID)),
C     :                 IGRID=1,NGRIDS),
C     :                 IF2=1,3)
C        WRITE(22,*) 'Longitudes of faces - all grids'
        WRITE(22) ((FLONG(IF1,IGRID),
     :                 IF1=1,NFACE(IGRID)),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Latitudes of faces - all grids'
        WRITE(22) ((FLAT(IF1,IGRID),
     :                 IF1=1,NFACE(IGRID)),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Longitudes of vertices - all grids'
        WRITE(22) ((VLONG(IV1,IGRID),
     :                 IV1=1,NVERT(IGRID)),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Latitudes of vertices - all grids'
        WRITE(22) ((VLAT(IV1,IGRID),
     :                 IV1=1,NVERT(IGRID)),
     :                 IGRID=1,NGRIDS)
C        WRITE(22,*) 'Areas of faces - all grids'
        WRITE(22) ((FAREA(IF1,IGRID),
     :                IF1=1,NFACE(IGRID)),
     :                IGRID=1,NGRIDS)
C        WRITE(22,*) 'Lengths of edges - all grids'
        WRITE(22) ((LDIST(IE1,IGRID),
     :                IE1=1,NEDGE(IGRID)),
     :                IGRID=1,NGRIDS)
C        WRITE(22,*) 'Distance between faces across edges - all grids'
        WRITE(22) ((DDIST(IE1,IGRID),
     :                IE1=1,NEDGE(IGRID)),
     :                IGRID=1,NGRIDS)


C
        CLOSE(22)
C
      ENDIF
C
C
C     Create file of grid coordinates only, for plotting
      OPEN(44,FILE='primalgrid.dat',FORM='FORMATTED')
      DO 200 1 ie1 = 1, nedgex
        iv1 = vofe(ie1,1,ngrids)
        long = vlong(iv1,ngrids)
        lat = vlat(iv1,ngrids)
        CALL ll2xyz(long,lat,x1,y1,z1)
        iv1 = vofe(ie1,2,ngrids)
        long = vlong(iv1,ngrids)
        lat = vlat(iv1,ngrids)
        CALL ll2xyz(long,lat,x2,y2,z2)
        WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
 2001 CONTINUE
      CLOSE(44)
      OPEN(44,FILE='dualgrid.dat',FORM='FORMATTED')
      DO 2002 ie1 = 1, nedgex
        if1 = fnxte(ie1,1,ngrids)
        long = flong(if1,ngrids)
        lat = flat(if1,ngrids)
        CALL ll2xyz(long,lat,x1,y1,z1)
        if1 = fnxte(ie1,2,ngrids)
        long = flong(if1,ngrids)
        lat = flat(if1,ngrids)
        CALL ll2xyz(long,lat,x2,y2,z2)
        WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
 2002 CONTINUE
      CLOSE(44)
C
C
C     -------------------------------------------------------------------
C
C     Check grid cross-references by writing out coords for plotting
C
      goto 503
      OPEN(82,FILE='grid.coords')
C
      DO 501 IGRID=1,NGRIDS
C
C     VLONG(IV1,IGRID) and VLAT(IV1,IGRID) contain the latitude
C     and longitude of the IV1'th vertex on grid IGRID
C
C     1) plot edges of dodecahedron
      WRITE(82,*) NEDGE(IGRID)
      DO 500 IEDGE=1,NEDGE(IGRID)
        IV1=VOFE(IEDGE,1,IGRID)
        IV2=VOFE(IEDGE,2,IGRID)
        LONG=VLONG(IV1,IGRID)
        LAT=VLAT(IV1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=VLONG(IV2,IGRID)
        LAT=VLAT(IV2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        WRITE (82,598) X1,Y1,Z1,X2,Y2,Z2
  500 CONTINUE
C
  501 CONTINUE
C
C     Plot pairs of lines crossing each edge
      IGRID=NGRIDS
      WRITE(82,*) NEDGE(IGRID)
      DO 502 IEDGE=1,NEDGE(IGRID) 
c        IV1=FNXTE(IEDGE,1,IGRID)
c        IV2=FNXTE(IEDGE,2,IGRID)
c        LONG=FLONG(IV1,IGRID)
c        LAT=FLAT(IV1,IGRID)
c        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
c        LONG=FLONG(IV2,IGRID)
c        LAT=FLAT(IV2,IGRID)
c        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
c        WRITE (82,598) X1,Y1,Z1,X2,Y2,Z2
        IV1=FEOFE(IEDGE,1,IGRID)
        IV2=FEOFE(IEDGE,2,IGRID)
        LONG=FLONG(IV1,IGRID)
        LAT=FLAT(IV1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=FLONG(IV2,IGRID)
        LAT=FLAT(IV2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        WRITE (82,598) X1,Y1,Z1,X2,Y2,Z2
  502 CONTINUE
C
      DO 504 IGRID=1,NGRIDS
C
C     FLONG(IF1,IGRID) and FLAT(IF1,IGRID) contain the latitude
C     and longitude of the IF1'th face centre on grid IGRID
C
C     2) plot edges of icosahedron
      WRITE(82,*) NEDGE(IGRID)
      DO 510 IF1=1,NFACE(IGRID)
        LONG=FLONG(IF1,IGRID)
        LAT=FLAT(IF1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        IF (IF1.LE.12) THEN
          NF2=5
        ELSE
          NF2=6
        ENDIF
        DO 520 IF2=1,NF2
          IF3=FNXTF(IF1,IF2,IGRID)
          IF (IF3.GT.IF1) THEN
            LONG=FLONG(IF3,IGRID)
            LAT=FLAT(IF3,IGRID)
            CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
            WRITE(82,598) X1,Y1,Z1,X2,Y2,Z2
          ENDIF
  520   CONTINUE
  510 CONTINUE
  504 CONTINUE
C
  598 FORMAT(6F9.4)
C
      CLOSE(82)
  503 continue
C
C     ------------------------------------------------------------------
C
      STOP
      END
C
C     ==================================================================
C
      SUBROUTINE LL2XYZ(LONG,LAT,X,Y,Z)
C
C     To convert longitude and latitude to cartesian coordinates
C     on the unit sphere
C
      IMPLICIT NONE
C
      REAL*8 LONG,LAT,X,Y,Z,CLN,SLN,CLT,SLT
C
C     ------------------------------------------------------------------
C
      SLN=SIN(LONG)
      CLN=COS(LONG)
      SLT=SIN(LAT)
      CLT=COS(LAT)
C
      X=CLN*CLT
      Y=SLN*CLT
      Z=SLT
C
C     ------------------------------------------------------------------
C
      RETURN
      END
C
C     ==================================================================
C
      SUBROUTINE XYZ2LL(X,Y,Z,LONG,LAT)
C
C     To convert cartesian coordinates to longitude and latitude
C
      IMPLICIT NONE
C   
      REAL*8 X,Y,Z,LONG,LAT,PI,TLN,TLT,R
C
C     -------------------------------------------------------------------
C
      PI=4.0D0*ATAN(1.0D0)
C
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
C
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
C
C     --------------------------------------------------------------------
C
      RETURN
      END
C
C     ====================================================================
C
      SUBROUTINE ADDTAB(TAB,INDEX,ENTRY,DIM1,DIM2)
C
C     Subroutine to add an entry to a table
C
      INTEGER DIM1,DIM2,TAB(DIM1,DIM2),INDEX,ENTRY
C
C     --------------------------------------------------------------------
C
      I=0
C
  100 CONTINUE
      I=I+1
      IF (I.GT.DIM2) THEN
        PRINT *,'**********'
        PRINT *,'TABLE FULL'
        PRINT *,'**********'
        STOP
      ENDIF
      IF (TAB(INDEX,I).NE.0) GOTO 100
      TAB(INDEX,I)=ENTRY
C
C     ---------------------------------------------------------------------
C
      RETURN
      END
C
C     =====================================================================
C
      SUBROUTINE BRENT(AX,BX,CX,TOL,N,P,XIT,XMIN,FMIN,IF1,IGRID,NGCHCK)
C
C     Perform line minimization of the function F using Brent's method.
C     Based on Numerical Recipes
C     AX,BX,CX is assumed to bracket the minimum, i.e. F(BX)<F(AX),F(CX)
C
      IMPLICIT NONE
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
C
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS),NGCHCK
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
C
      INTEGER ITER,ITMAX,N,J,IF1,IGRID
      REAL*8 CGOLD,AX,BX,CX,TOL,XMIN,A,B,X,U,V,W,FX,FU,FV,FW,TOL1,TOL2,
     R       FMIN,E,D,ZEPS,PP,Q,R,ETEMP,XM,P(N),XIT(N)
C
C     Check resolution
      IF (NGRIDS.NE.NGCHCK) THEN
        PRINT *,'ERROR'
        PRINT *,'NGRIDS=',NGRIDS,' in routine BRENT but'
        PRINT *,'NGRIDS=',NGCHCK,' in the calling routine.'
        STOP
      ENDIF
C
C     Maximum number of iterations
      ITMAX=100
C
C     Golden ratio
      CGOLD=0.3819660D0
C
C     Protect against divide by zero
      ZEPS=1.0D-10
C
C     Sort A, B, into ascending order
      IF (AX.LT.BX) THEN
        A=AX
        B=CX
      ELSE
        A=CX
        B=AX
      ENDIF
C
C     Initialize search points and function values
      X=BX
      W=BX
      V=BX
      FLONG(IF1,IGRID)=P(1)+X*XIT(1)
      FLAT(IF1,IGRID)=P(2)+X*XIT(2)
      CALL PEN(IF1,IGRID,FX)
      FW=FX
      FV=FX
C
      E=0.0D0
      D=0.0D0
C
C     MAIN LOOP
      ITER=0
  100 CONTINUE
C
      ITER=ITER+1
      XM=0.5D0*(A+B)
C
C     Check for convergence
C      TOL1=TOL*ABS(X)+ZEPS
C     In the case where the min happens to fall at x=0 but P is
C     non-zero, better to put a typical P value (1.0) instead of x
      TOL1=TOL+ZEPS
      TOL2=2.0D0*TOL1
      IF (ABS(X-XM).LE.TOL2-0.5D0*(B-A)) GOTO 900
C
C     Construct a trial parabolic fit
      IF (ABS(E).GT.TOL1) THEN
        R=(X-W)*(FX-FV)
        Q=(X-V)*(FX-FW)
        PP=(X-V)*Q-(X-W)*R
        Q=2.0D0*(Q-R)
        IF (Q.GT.0.0D0) PP=-PP
        Q=ABS(Q)
        ETEMP=E
        E=D
        IF (ABS(PP).GE.ABS(0.5D0*Q*ETEMP).OR.
     :      PP.LE.Q*(A-X).OR.
     :      PP.GE.Q*(B-X)) THEN
          IF (X.GE.XM) THEN
            E=A-X
          ELSE
            E=B-X
          ENDIF
          D=CGOLD*E
        ELSE
          D=PP/Q
          U=X+D
          IF (U-A.LT.TOL2.OR.B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
        ENDIF
      ELSE
        IF (X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
      ENDIF
C
      IF (ABS(D).GE.TOL1) THEN
        U=X+D
      ELSE
        U=X+SIGN(TOL1,D)
      ENDIF
      FLONG(IF1,IGRID)=P(1)+U*XIT(1)
      FLAT(IF1,IGRID)=P(2)+U*XIT(2)
      CALL PEN(IF1,IGRID,FU)
C
      IF (FU.LE.FX) THEN
        IF (U.GE.X) THEN
          A=X
        ELSE
          B=X
        ENDIF
        V=W
        W=X
        X=U
        FV=FW
        FW=FX
        FX=FU
      ELSE
        IF (U.LT.X) THEN
          A=U
        ELSE
          B=U
        ENDIF
        IF (FU.LE.FW.OR.W.EQ.X) THEN
          V=W
          W=U
          FV=FW
          FW=FU
        ELSEIF (FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
          V=U
          FV=FU
        ENDIF
      ENDIF
C
      IF (ITER.LT.ITMAX) GOTO 100
C
c      PRINT *,'Maximum iterations exceeded in BRENT'
c      STOP
C
  900 CONTINUE
      XMIN=X
      FMIN=FX
C
C     Save vector displacement and new position of min
      DO J=1,N
        XIT(J)=XMIN*XIT(J)
        P(J)=P(J)+XIT(J)
      ENDDO
C
c      print *,'Brent: ',iter,' iterations '
C
      RETURN
      END
C
C     =================================================================
C
      SUBROUTINE POWELL(IF1,IGRID,FTOL,BDISP,MXDISP,NGCHCK)
C
C     Minimize a function in 2 dimensions using Powell's method.
C     Following Numerical Recipes
C
      IMPLICIT NONE
C
      INTEGER N
      PARAMETER(N=2)
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
C
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS),NGCHCK
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
C
      INTEGER J,ITER,ITMAX,IBIG,I,IF1,IGRID,IV1
      REAL*8 P(2),XI(2,2),PT(2),PTT(2),XIT(2),FTOL,FP,FPTT,AX,BX,CX,TOL,
     R       DEL,FRET,T,FAC,CGOLD,XMIN,MXDISP,BDISP
C
C
C     Check resolution
      IF (NGRIDS.NE.NGCHCK) THEN
        PRINT *,'ERROR'
        PRINT *,'NGRIDS=',NGRIDS,' in routine POWELL but'
        PRINT *,'NGRIDS=',NGCHCK,' in the calling routine.'
        STOP
      ENDIF
C
C     Golden ratio
      CGOLD=0.3819660D0
C
      MXDISP=0.0D0
C
C     Brackets for 1D search
c      FAC=2.0D0**(-IGRID)
c      DEL=1.18D0*FAC
cc      DEL=0.01D0*FAC
c      AX=-DEL
c      BX=0.0
cC      BX=(2.0D0*CGOLD-1.0D0)*DEL
c      CX=DEL
      AX=-BDISP
      BX=0.0D0
      CX=BDISP
C
C     Tolerance for 1d search
C      TOL=1.0D-8
      TOL=0.01D0*BDISP
C
      P(1)=FLONG(IF1,IGRID)
      P(2)=FLAT(IF1,IGRID)
      CALL PEN(IF1,IGRID,FRET)
C
C     Initialize search directions
      DO I=1,N
        DO J=1,N
          XI(J,I)=0.0D0
        ENDDO
      ENDDO
      XI(1,1)=1.0D0/COS(FLAT(IF1,IGRID))
      XI(2,2)=1.0D0
C
C     Save initial point
      DO J=1,N
        PT(J)=P(J)
      ENDDO
C
      ITMAX=6
C
C     MAIN LOOP
      ITER=0
  100 CONTINUE
      ITER=ITER+1
      FP=FRET
      IBIG=0
      DEL=0.0D0
C
C     Loop over all directions
      DO I=1,N
C
C       Copy the direction
        DO J=1,N
          XIT(J)=XI(J,I)
        ENDDO
        FPTT=FRET
C
C       Line minimization along direction XIT from P
        CALL BRENT(AX,BX,CX,TOL,N,P,XIT,XMIN,FRET,IF1,IGRID,NGRIDS)
        MXDISP=MAX(MXDISP,XMIN)
C
        IF (ABS(FPTT-FRET).GT.DEL) THEN
          DEL=ABS(FPTT-FRET)
          IBIG=I
        ENDIF
C
      ENDDO
C
C     Are we finished?
      IF (2.0D0*ABS(FP-FRET).LE.FTOL*(ABS(FP)+ABS(FRET))) GOTO 900
C
C     Construct extrapolated point and average direction moved
      DO J=1,N
        PTT(J)=2.0D0*P(J)-PT(J)
        XIT(J)=P(J)-PT(J)
        PT(J)=P(J)
      ENDDO
C
C     Function value at extrapolated point
      FLONG(IF1,IGRID)=PTT(1)
      FLAT(IF1,IGRID)=PTT(2)
      CALL PEN(IF1,IGRID,FPTT)
C
      IF (FPTT.LT.FP) THEN
        T=2.0D0*(FP-2.0D0*FRET+FPTT)*SQRT(FP-FRET-DEL)-DEL*SQRT(FP-FPTT)
        IF (T.LT.0.0D0) THEN
          CALL BRENT(AX,BX,CX,TOL,N,P,XIT,XMIN,FRET,IF1,IGRID,NGRIDS)
          MXDISP=MAX(MXDISP,XMIN)
          DO J=1,N
            XI(J,IBIG)=XI(J,N)
            XI(J,N)=XIT(J)
          ENDDO
        ENDIF
      ENDIF
C
C
      IF (ITER.LT.ITMAX) GOTO 100
C
C
  900 CONTINUE
C
C     Reset long and lat to best value
      FLONG(IF1,IGRID)=P(1)
      FLAT(IF1,IGRID)=P(2)
C
C     And reset positions of affected vertices
      DO I=1,6
        IV1=VOFF(IF1,I,IGRID)
        CALL VRTCRD(IV1,IGRID)
      ENDDO
C
      CALL PEN(IF1,IGRID,FRET)


      RETURN
      END
C
C     =================================================================
C
      SUBROUTINE VRTCRD(IV1,IGRID)
C
C     Find the coordinates of a vertex
C
      IMPLICIT NONE
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX,N
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
      PARAMETER(N=2)
C
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS)
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
C
      INTEGER IV1,IGRID,IF1
      REAL*8 LONG,LAT,X,Y,Z,X1,Y1,Z1,D1X,D1Y,D1Z,D2X,D2Y,D2Z,
     :       VX,VY,VZ,MAG
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
C
C
      IF1=FOFV(IV1,1,IGRID)
      LONG=FLONG(IF1,IGRID)
      LAT=FLAT(IF1,IGRID)
      CALL LL2XYZ(LONG,LAT,X,Y,Z)
      IF1=FOFV(IV1,2,IGRID)
      LONG=FLONG(IF1,IGRID)
      LAT=FLAT(IF1,IGRID)
      CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
      D1X=X-X1
      D1Y=Y-Y1
      D1Z=Z-Z1
      IF1=FOFV(IV1,3,IGRID)
      LONG=FLONG(IF1,IGRID)
      LAT=FLAT(IF1,IGRID)
      CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
      D2X=X-X1
      D2Y=Y-Y1
      D2Z=Z-Z1
      VX=(D1Y*D2Z-D2Y*D1Z)
      VY=(D1Z*D2X-D2Z*D1X)
      VZ=(D1X*D2Y-D2X*D1Y)
C     Make sure it's in the right hemisphere
      IF ((VX*X+VY*Y+VZ*Z).LT.0.0D0) THEN
        VX=-VX
        VY=-VY
        VZ=-VZ
      ENDIF
C     Project back onto the sphere
      MAG=SQRT(VX*VX+VY*VY+VZ*VZ)
      VX=VX/MAG
      VY=VY/MAG
      VZ=VZ/MAG
C     Convert back to latitude/longitude
      CALL XYZ2LL(VX,VY,VZ,LONG,LAT)
      VLONG(IV1,IGRID)=LONG
      VLAT(IV1,IGRID)=LAT
C
      RETURN
      END
C
C     =======================================================================
C
      SUBROUTINE PENALT(IE1,IGRID,COST,DD)
C
C     Find the contribution to the Heikes+Randall penalty function
C     from edge IE1
C
      IMPLICIT NONE
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX,N
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
      PARAMETER(N=2)
C
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS)
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
C
      INTEGER IE1,IGRID,IV1,IF1
      REAL*8 LONG,LAT,X1,Y1,Z1,X2,Y2,Z2,XM,YM,ZM,X,Y,Z,DX,DY,DZ,L2,
     :       R2,COST,MAG,D2,DD
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
C
C
C     First find the ends of the edge
      IV1=VOFE(IE1,1,IGRID)
      LONG=VLONG(IV1,IGRID)
      LAT=VLAT(IV1,IGRID)
      CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
      IV1=VOFE(IE1,2,IGRID)
      LONG=VLONG(IV1,IGRID)
      LAT=VLAT(IV1,IGRID)
      CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
C
C     Calculate length
      DX=X1-X2
      DY=Y1-Y2
      DZ=Z1-Z2
      L2=DX*DX+DY*DY+DZ*DZ
C
C     And midpoint
      XM=0.5D0*(X1+X2)
      YM=0.5D0*(Y1+Y2)
      ZM=0.5D0*(Z1+Z2)
      MAG=1.0D0/SQRT(XM*XM+YM*YM+ZM*ZM)
      XM=XM*MAG
      YM=YM*MAG
      ZM=ZM*MAG
C
C     Find faces either side of edge
      IF1=FNXTE(IE1,1,IGRID)
      LONG=FLONG(IF1,IGRID)
      LAT=FLAT(IF1,IGRID)
      CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
      IF1=FNXTE(IE1,2,IGRID)
      LONG=FLONG(IF1,IGRID)
      LAT=FLAT(IF1,IGRID)
      CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
C
C     Find midpoint
      X=0.5D0*(X1+X2)
      Y=0.5D0*(Y1+Y2)
      Z=0.5D0*(Z1+Z2)
      MAG=1.0D0/SQRT(X*X+Y*Y+Z*Z)
      X=X*MAG
      Y=Y*MAG
      Z=Z*MAG
C
C     Contribution to penalty function
      DX=X-XM
      DY=Y-YM
      DZ=Z-ZM
      D2=DX*DX+DY*DY+DZ*DZ
      DD=SQRT(D2)
      R2=D2/L2
      COST=R2*R2
c      print *,'Edge ',IE1,'  displacement ',DD
C
      RETURN
      END
C
C     =======================================================================
C
      SUBROUTINE PENALF1(IF1,IGRID,TOTCST)
C
C     Find the contribution to the Heikes+Randall penalty function
C     affected by face IF1
C
      IMPLICIT NONE
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX,N
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
      PARAMETER(N=2)
C
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS)
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
C
      INTEGER IGRID,IF1,IV1,IE1,JE1,JE2,IELIST(12),I,J
      REAL*8 COST,TOTCST,DD
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
C
C
      IF (IF1 .LE. 12) THEN
        PRINT *,'Attempt to find cost function PENALF1 for pentagon'
        STOP
      ENDIF
C
C     Find vertices that are affected and update their coordinates
C     (IF1 will never be a pentagon, so only consider hexagons)
C     List edges that are affected by face IF1
      JE1=0
      DO I=1,6
        IV1=VOFF(IF1,I,IGRID)
        CALL VRTCRD(IV1,IGRID)
        DO J=1,3
          IE1=EOFV(IV1,J,IGRID)
C         Check we haven't got this one already
          JE2=1
  300     CONTINUE
          IF (JE2.GT.JE1) GOTO 100
          IF (IELIST(JE2).EQ.IE1) GOTO 200
          JE2=JE2+1
          GOTO 300
  100     CONTINUE
C         Haven't got it so add it to the list
          JE1=JE1+1
          IELIST(JE1)=IE1
  200     CONTINUE
        ENDDO
      ENDDO
C
C     Now add up contributions to penalty function
      TOTCST=0.0D0
      DO JE1=1,12
        IE1=IELIST(JE1)
        CALL PENALT(IE1,IGRID,COST,DD)
        TOTCST=TOTCST+COST
      ENDDO
C
      RETURN
      END
C
C     ===================================================================
C
      SUBROUTINE PENALC(IF1,IGRID,COST,DD)
C
C     Find the contribution to the penalty function that measures departures
C     from centroidal for face IF1
C     PENALF1 should be called first to ensure that vertices are up to date
C
      IMPLICIT NONE
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX,N
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
      PARAMETER(N=2)
C
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS)
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
C
      INTEGER IGRID,IF1,IV1,IV2,IE1,IE2,NE
      REAL*8 COST,LONG,LAT,X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AFACE,
     R       XC,YC,ZC,R2,DX,DY,DZ,MAG,DA,D2,DD
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
C
C
      IF (IF1 .LE. 12) THEN
        NE = 5
      ELSE
        NE = 6
      ENDIF
C
C     Coordinates of centre of face
      LONG=FLONG(IF1,IGRID)
      LAT=FLAT(IF1,IGRID)
      CALL LL2XYZ(LONG,LAT,X0,Y0,Z0)
C     Loop over edges in turn and calculate area of triangle
C     formed by the edge and the centre of the face
C     Hence find area of face and centroid
      XC=0.0D0
      YC=0.0D0
      ZC=0.0D0
      AFACE=0.0D0
      DO 810 IE1=1,NE
        IE2=EOFF(IF1,IE1,IGRID)
        IV1=VOFE(IE2,1,IGRID)
        IV2=VOFE(IE2,2,IGRID)
        LONG=VLONG(IV1,IGRID)
        LAT=VLAT(IV1,IGRID)
        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
        LONG=VLONG(IV2,IGRID)
        LAT=VLAT(IV2,IGRID)
        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
        CALL STAREA2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,DA)
        AFACE=AFACE+DA
        XC=XC+(X0+X1+X2)*DA/3.0D0
        YC=YC+(Y0+Y1+Y2)*DA/3.0D0
        ZC=ZC+(Z0+Z1+Z2)*DA/3.0D0
  810 CONTINUE
      MAG=SQRT(XC*XC+YC*YC+ZC*ZC)
      XC=XC/MAG
      YC=YC/MAG
      ZC=ZC/MAG
C
C     Contribution to penalty function
      DX=X0-XC
      DY=Y0-YC
      DZ=Z0-ZC
      D2=DX*DX+DY*DY+DZ*DZ
      DD=SQRT(D2)
      R2=D2/AFACE
      COST=R2*R2
C
      RETURN
      END
C
C     ===================================================================
C
      SUBROUTINE PENALF2(IF1,IGRID,TOTCST)
C
C     Find the contribution to the centroidal penalty function
C     affected by face IF1
C
      IMPLICIT NONE
C
      INTEGER NGRIDS,NFACEX,NEDGEX,NVERTX,N
      PARAMETER(NGRIDS=3)
      PARAMETER(
     :          NFACEX=5*2**(2*NGRIDS-1)+2,
     :          NEDGEX=15*2**(2*NGRIDS-1),
     :          NVERTX=5*2**(2*NGRIDS))
      PARAMETER(N=2)
C
      INTEGER FNXTF(NFACEX,6,NGRIDS),EOFF(NFACEX,6,NGRIDS),
     :        FNXTE(NEDGEX,2,NGRIDS),FEOFE(NEDGEX,2,NGRIDS),
     :        VOFE(NEDGEX,2,NGRIDS),VOFF(NFACEX,6,NGRIDS),
     :        FOFV(NVERTX,3,NGRIDS),EOFV(NVERTX,3,NGRIDS)
      REAL*8 FLONG(NFACEX,NGRIDS),FLAT(NFACEX,NGRIDS),
     :       FAREA(NFACEX,NGRIDS),
     :       VLONG(NVERTX,NGRIDS),VLAT(NVERTX,NGRIDS),
     :       LDIST(NEDGEX,NGRIDS),DDIST(NEDGEX,NGRIDS),
     :       GDIST(NEDGEX,NGRIDS),COEFF(NVERTX,3,NGRIDS)
C
      INTEGER IGRID,IF1,IF2,IFLIST(7),I
      REAL*8 COST,TOTCST,DD
C
      COMMON /COMGRD/ FNXTF,EOFF,FNXTE,FEOFE,VOFE,VOFF,FOFV,EOFV,
     :                FLONG,FLAT,FAREA,VLONG,VLAT,LDIST,DDIST,GDIST,
     :                COEFF
C
C
      IF (IF1 .LE. 12) THEN
        PRINT *,'Attempt to find cost function PENALF2 for pentagon'
        STOP
      ENDIF
C
C     Find neighbouring faces, i.e. those whose contribution is
C     affected by IF1
C     (IF1 will never be a pentagon, so only consider hexagons)
C     List faces that are affected by face IF1
      DO I=1,6
        IF2=FNXTF(IF1,I,IGRID)
        IFLIST(I)=IF2
      ENDDO
      IFLIST(7)=IF1
C
C     Now add up contributions to penalty function
      TOTCST=0.0D0
      DO I=1,7
        IF2=IFLIST(I)
        CALL PENALC(IF2,IGRID,COST,DD)
        TOTCST=TOTCST+COST
      ENDDO
C
      RETURN
      END
C
C     ===================================================================
CC
      SUBROUTINE PEN(IF1,IGRID,TOTCST)
C
C     Find the contribution to the penalty function affected by face IF1
C     Weighted contributions from the Heikes-Randall penalty
C     function and the centroidal penalty function
C
      IMPLICIT NONE
      INTEGER IF1,IGRID
      REAL*8 C1,C2,TOTCST,DD
      REAL*8 WEIGHT1,WEIGHT2
C
      COMMON /COMWGT/ WEIGHT1,WEIGHT2
C
C
      CALL PENALF1(IF1,IGRID,C1)
      CALL PENALF2(IF1,IGRID,C2)
      TOTCST=WEIGHT1*C1 + WEIGHT2*C2
C
      RETURN
      END
C
C     ===================================================================
C
      SUBROUTINE STAREA2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AREA)
C
C     Calculate the area of the spherical triangle whose corners
C     have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
C     The formula below is more robust to roundoff error than the
C     better known sum of angle - PI formula
C
      IMPLICIT NONE
C
      REAL*8 X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,
     R       D0,D1,D2,S,T0,T1,T2,T3,AREA
C
C
C     Distances between pairs of points
      CALL SPDIST(X0,Y0,Z0,X1,Y1,Z1,D2)
      CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,D0)
      CALL SPDIST(X2,Y2,Z2,X0,Y0,Z0,D1)
C
C     Half perimeter
      S=0.5D0*(D0+D1+D2)
C
C     Tangents
      T0 = TAN(0.5D0*(S-D0))
      T1 = TAN(0.5D0*(S-D1))
      T2 = TAN(0.5D0*(S-D2))
      T3 = TAN(0.5D0*S)
C
C     Area
      AREA = 4.0D0*ATAN(SQRT(T0*T1*T2*T3))
C
      RETURN
      END
C
C     ===================================================================
C
      SUBROUTINE SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
C
C     Calculate the spherical distance S between two points with Cartesian
C     coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

      IMPLICIT NONE

      REAL*8 X1, Y1, Z1, X2, Y2, Z2, S, DX, DY, DZ, AD


      DX = X2 - X1
      DY = Y2 - Y1
      DZ = Z2 - Z1
      AD = SQRT(DX*DX + DY*DY + DZ*DZ)
      S = 2.0D0*ASIN(0.5D0*AD)


      RETURN
      END
C
C     ===================================================================


