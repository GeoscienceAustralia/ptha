! Modified from https://github.com/cphyc/Fortran-parallel-sort
!
! It is released under the very permissive licence:
!    Creative Commons Zero v1.0 Universal
! https://creativecommons.org/publicdomain/zero/1.0/


!Subroutine D_mrgrnk (array, inds)
    ! __________________________________________________________
    !   MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
    SORT_INDEX_TYPE , Intent (In) :: array(n)
    Integer(c_int), Intent (Out) :: inds(n)
    integer(c_int), intent(in) :: n
    ! __________________________________________________________
    SORT_INDEX_TYPE2 :: XVALA, XVALB
    !
    Integer(c_int), allocatable :: JWRKT(:)
    Integer(c_int) :: LMTNA, LMTNC, IRNG1, IRNG2
    Integer(c_int) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

    !
    NVAL = Min (SIZE(array), SIZE(inds))
    Select Case (NVAL)
    Case (:0)
       Return
    Case (1)
       inds (1) = 1
       Return
    Case Default
       Continue
    End Select
    allocate(JWRKT(n))
    !
    !  Fill-in the index array, creating ordered couples
    !
    Do IIND = 2, NVAL, 2
       If (array(IIND-1) <= array(IIND)) Then
          inds (IIND-1) = IIND - 1
          inds (IIND) = IIND
       Else
          inds (IIND-1) = IIND
          inds (IIND) = IIND - 1
       End If
    End Do
    If (Modulo(NVAL, 2) /= 0) Then
       inds (NVAL) = NVAL
    End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    LMTNA = 2
    LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
    Do
       If (NVAL <= 2) Exit
       !
       !   Loop on merges of A and B into C
       !
       Do IWRKD = 0, NVAL - 1, 4
          If ((IWRKD+4) > NVAL) Then
             If ((IWRKD+2) >= NVAL) Exit
             !
             !   1 2 3
             !
             If (array(inds(IWRKD+2)) <= array(inds(IWRKD+3))) Exit
             !
             !   1 3 2
             !
             If (array(inds(IWRKD+1)) <= array(inds(IWRKD+3))) Then
                IRNG2 = inds (IWRKD+2)
                inds (IWRKD+2) = inds (IWRKD+3)
                inds (IWRKD+3) = IRNG2
                !
                !   3 1 2
                !
             Else
                IRNG1 = inds (IWRKD+1)
                inds (IWRKD+1) = inds (IWRKD+3)
                inds (IWRKD+3) = inds (IWRKD+2)
                inds (IWRKD+2) = IRNG1
             End If
             Exit
          End If
          !
          !   1 2 3 4
          !
          If (array(inds(IWRKD+2)) <= array(inds(IWRKD+3))) Cycle
          !
          !   1 3 x x
          !
          If (array(inds(IWRKD+1)) <= array(inds(IWRKD+3))) Then
             IRNG2 = inds (IWRKD+2)
             inds (IWRKD+2) = inds (IWRKD+3)
             If (array(IRNG2) <= array(inds(IWRKD+4))) Then
                !   1 3 2 4
                inds (IWRKD+3) = IRNG2
             Else
                !   1 3 4 2
                inds (IWRKD+3) = inds (IWRKD+4)
                inds (IWRKD+4) = IRNG2
             End If
             !
             !   3 x x x
             !
          Else
             IRNG1 = inds (IWRKD+1)
             IRNG2 = inds (IWRKD+2)
             inds (IWRKD+1) = inds (IWRKD+3)
             If (array(IRNG1) <= array(inds(IWRKD+4))) Then
                inds (IWRKD+2) = IRNG1
                If (array(IRNG2) <= array(inds(IWRKD+4))) Then
                   !   3 1 2 4
                   inds (IWRKD+3) = IRNG2
                Else
                   !   3 1 4 2
                   inds (IWRKD+3) = inds (IWRKD+4)
                   inds (IWRKD+4) = IRNG2
                End If
             Else
                !   3 4 1 2
                inds (IWRKD+2) = inds (IWRKD+4)
                inds (IWRKD+3) = IRNG1
                inds (IWRKD+4) = IRNG2
             End If
          End If
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 4
       Exit
    End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
    Do
       If (LMTNA >= NVAL) Exit
       IWRKF = 0
       LMTNC = 2 * LMTNC
       !
       !   Loop on merges of A and B into C
       !
       Do
          IWRK = IWRKF
          IWRKD = IWRKF + 1
          JINDA = IWRKF + LMTNA
          IWRKF = IWRKF + LMTNC
          If (IWRKF >= NVAL) Then
             If (JINDA >= NVAL) Exit
             IWRKF = NVAL
          End If
          IINDA = 1
          IINDB = JINDA + 1
          !
          !   Shortcut for the case when the max of A is smaller
          !   than the min of B. This line may be activated when the
          !   initial set is already close to sorted.
          !
          !          IF (array(inds(JINDA)) <= array(inds(IINDB))) CYCLE
          !
          !  One steps in the C subset, that we build in the final rank array
          !
          !  Make a copy of the rank array for the merge iteration
          !
          JWRKT (1:LMTNA) = inds (IWRKD:JINDA)
          !
          XVALA = array (JWRKT(IINDA))
          XVALB = array (inds(IINDB))
          !
          Do
             IWRK = IWRK + 1
             !
             !  We still have unprocessed values in both A and B
             !
             If (XVALA > XVALB) Then
                inds (IWRK) = inds (IINDB)
                IINDB = IINDB + 1
                If (IINDB > IWRKF) Then
                   !  Only A still with unprocessed values
                   inds (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                   Exit
                End If
                XVALB = array (inds(IINDB))
             Else
                inds (IWRK) = JWRKT (IINDA)
                IINDA = IINDA + 1
                If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                XVALA = array (JWRKT(IINDA))
             End If
             !
          End Do
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 2 * LMTNA
    End Do
    !
    deallocate(JWRKT)
    Return
    !
  !End Subroutine D_mrgrnk
