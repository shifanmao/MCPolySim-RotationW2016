!---------------------------------------------------------------*

! Calculate the change in the twist energy for a MC move
! Kevin J Hou, Winter 2016
!
! Adapted from MC_elas code
! 
     
      SUBROUTINE MC_twist(DETWIST,R,U,W,RP,UP,WP,NT,N,NP,IT1,IT2,TAU)
      
        DOUBLE PRECISION R(NT,3)   ! Bead positions
        DOUBLE PRECISION U(NT,3)   ! Tangent vectors
        DOUBLE PRECISION W(NT,3)   ! Twist Vectors
        DOUBLE PRECISION RP(NT,3)  ! Bead positions
        DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
        DOUBLE PRECISION WP(NT,3)  ! Twist vectors
        INTEGER N,NP,NT           ! Number of beads

        INTEGER IT1               ! Index of test bead 1
        INTEGER IT2               ! Index of test bead 2
      
        DOUBLE PRECISION DETWIST  ! Change in ETWIST
        DOUBLE PRECISION TAU      ! Parameter describing twist energy

        DOUBLE PRECISION OLDOMEGA    ! Old Twist Angle
        DOUBLE PRECISION NEWOMEGA    ! New Twist Angle
        DOUBLE PRECISION OLDNUM      ! Numerator of Old Twist Angle Expression
        DOUBLE PRECISION OLDDEN      ! Denominator of Old Twist Angle Expression
        DOUBLE PRECISION NEWNUM      ! Numerator of New Twist Angle Expression
        DOUBLE PRECISION NEWDEN      ! Denominator of New Twist Angle Expression
        DOUBLE PRECISION F(NT,3)     ! F = U x W, third orthogonal vector
        DOUBLE PRECISION FP(NT,3)    ! F = U x W, third orthogonal vector

 ! Setup

        DETWIST = 0.

! Calculate change in energy

        ! Calculate F = U X W for changed beads
        DO 10 I = MAX(1,IT1-1),IT2
             ! Obtain third vector by crossing U and W
             F(I,1) =  U(I,2)*W(I,3) - U(I,3)-W(I,2)
             F(I,2) = -U(I,1)*W(I,3) + U(I,3)*W(I,1)
             F(I,3) =  U(I,1)*W(I,2) - U(I,2)*W(I,1)
             FP(I,1) =  UP(I,2)*WP(I,3) - UP(I,3)-WP(I,2)
             FP(I,2) = -UP(I,1)*WP(I,3) + UP(I,3)*WP(I,1)
             FP(I,3) =  UP(I,1)*WP(I,2) - UP(I,2)*WP(I,1)
 10     CONTINUE


        DO 20 I = MAX(1,IT1-1),IT2-1
             OLDNUM = (W(I,1)*F(I+1,1)+W(I,2)*F(I+1,2)+W(I,3)*F(I+1,3)) &
             - (W(I+1,1)*F(I,1)+W(I+1,2)*F(I,2)+W(I+1,3)*F(I,3))
             OLDDEN = (F(I,1)*F(I+1,1)+F(I,2)*F(I+1,2)+F(I,3)*F(I+1,3)) &
             + (W(I+1,1)*W(I,1)+W(I+1,2)*W(I,2)+W(I+1,3)*W(I,3))
             NEWNUM = (WP(I,1)*FP(I+1,1)+WP(I,2)*FP(I+1,2)+WP(I,3)*FP(I+1,3)) &
             - (WP(I+1,1)*FP(I,1)+WP(I+1,2)*FP(I,2)+WP(I+1,3)*FP(I,3))
             NEWDEN = (FP(I,1)*FP(I+1,1)+FP(I,2)*FP(I+1,2)+FP(I,3)*FP(I+1,3)) &
             + (WP(I+1,1)*WP(I,1)+WP(I+1,2)*WP(I,2)+WP(I+1,3)*WP(I,3))

             OLDOMEGA = ATAN2(OLDNUM,OLDDEN)
             NEWOMEGA = ATAN2(NEWNUM,NEWDEN)

             DETWIST = DETWIST + (TAU/2)*(NEWOMEGA**2 - OLDOMEGA**2)
 20     CONTINUE

      !Check Loop Bounds (Change from IT-1 IT2 to whole chain and see if result is the same)
      
      RETURN      
      END
      
!---------------------------------------------------------------*
