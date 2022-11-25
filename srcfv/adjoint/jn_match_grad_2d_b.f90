!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of jn_match_grad_2d in reverse (adjoint) mode (with options with!SliceDeadControl with!SliceDeadInstrs with!St
!aticTaping):
!   gradient     of useful results: wr
!   with respect to varying inputs: wr wd
!   RW status of diff variables: wr:in-out wd:out
!===============================================================================
! JOIN Match for Gradients at centers
!===============================================================================
! =============================================================================
! Join match 2-D
! =============================================================================
SUBROUTINE JN_MATCH_GRAD_2D_B(wr, wrb, prr, gh1r, gh2r, gh3r, gh4r, imr&
& , jmr, wd, wdb, prd, gh1d, gh2d, gh3d, gh4d, imd, jmd, tr)
  IMPLICIT NONE
!
! Variables for dimension -----------------------------------------
  INTEGER, INTENT(IN) :: imr, jmr, imd, jmd
  INTEGER, INTENT(IN) :: gh1r, gh2r, gh3r, gh4r
  INTEGER, INTENT(IN) :: gh1d, gh2d, gh3d, gh4d
! Input variables -------------------------------------------------
  REAL*8, DIMENSION(1-gh1d:imd+gh2d, 1-gh3d:jmd+gh4d, 2), INTENT(IN) :: &
& wd
  REAL*8, DIMENSION(1-gh1d:imd+gh2d, 1-gh3d:jmd+gh4d, 2) :: wdb
  INTEGER, DIMENSION(2, 2), INTENT(IN) :: prr, prd
  INTEGER, DIMENSION(2), INTENT(IN) :: tr
! Output variables ------------------------------------------------
  REAL*8, DIMENSION(1-gh1r:imr+gh2r, 1-gh3r:jmr+gh4r, 2), INTENT(INOUT) &
& :: wr
  REAL*8, DIMENSION(1-gh1r:imr+gh2r, 1-gh3r:jmr+gh4r, 2), INTENT(INOUT) &
& :: wrb
! Local variables -------------------------------------------------
  INTEGER, POINTER :: ind1, ind2
  INTEGER :: ir, jr
  INTEGER, TARGET :: id, jd
  INTEGER :: istep, jstep, idir, jdir, idd, jdd, i, j
  INTEGER :: i1, i2, j1, j2
  INTRINSIC SIGN
  INTRINSIC ABS
  INTEGER :: branch
! -----------------------------------------------------------------
  istep = SIGN(1, tr(1))
  jstep = SIGN(1, tr(2))
  IF (tr(1) .GE. 0.) THEN
    idir = tr(1)
  ELSE
    idir = -tr(1)
  END IF
  IF (tr(2) .GE. 0.) THEN
    jdir = tr(2)
  ELSE
    jdir = -tr(2)
  END IF
!
  IF (istep .EQ. 1) THEN
    i1 = 1
    i2 = prd(2, idir) - prd(1, idir) + 1
  ELSE
    i1 = prd(2, idir) - prd(1, idir) + 1
    i2 = 1
  END IF
!
  IF (jstep .EQ. 1) THEN
    j1 = 1
    j2 = prd(2, jdir) - prd(1, jdir) + 1
  ELSE
    j1 = prd(2, jdir) - prd(1, jdir) + 1
    j2 = 1
  END IF
!
  IF (idir .EQ. 1) THEN
    CALL PUSHCONTROL1B(0)
    ind1 => id
    ind2 => jd
  ELSE IF (idir .EQ. 2) THEN
    CALL PUSHCONTROL1B(1)
    ind1 => jd
    ind2 => id
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
!
  idd = 1
  jdd = 1
!
  DO j=j1,j2,jstep
    DO i=i1,i2,istep
      CALL PUSHINTEGER4(ir)
      ir = prr(1, 1) + idd - 1
      idd = idd + 1
    END DO
    idd = 1
    CALL PUSHINTEGER4(jdd)
    jdd = jdd + 1
  END DO
  wdb = 0.0_8
  DO j=j2-MOD(j2-j1, jstep),j1,-jstep
    CALL POPINTEGER4(jdd)
    DO i=i2-MOD(i2-i1, istep),i1,-istep
      jr = prr(1, 2) + jdd - 1
      wdb(ind1, ind2, :) = wdb(ind1, ind2, :) + wrb(ir, jr, :)
      wrb(ir, jr, :) = 0.0_8
      CALL POPINTEGER4(ir)
    END DO
  END DO
  CALL POPCONTROL1B(branch)
END SUBROUTINE JN_MATCH_GRAD_2D_B
