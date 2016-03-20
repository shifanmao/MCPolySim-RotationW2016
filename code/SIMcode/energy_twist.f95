!---------------------------------------------------------------*
      
!     
!     This subroutine calculates the twist energy given the tangent (U)
!     and one of the twist (W) vectors.
!
!     Kevin J Hou
!     Written 2-17-16
      
      SUBROUTINE energy_twist(ETWIST,R,U,W,NT,N,NP,TAU)
      
      DOUBLE PRECISION ETWIST   ! Twist Energy
      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION W(NT,3)  ! Twist Vector
      DOUBLE PRECISION TAU      ! Twist Parameter
      INTEGER N,NT,NP           ! Number of bead

      DOUBLE PRECISION F(NT,3)  ! Other Twist Vector (U X W)
      DOUBLE PRECISION NUM      ! Numerator of Omega Expression
      DOUBLE PRECISION DEN      ! Denominator of Omega Expression
      DOUBLE PRECISION OMEGA    ! Twist Angle
      
      ! Initialize
      ETWIST = 0;

      ! Calculate third vector
      DO 10 I = 1,NT
             ! Obtain third vector by crossing U and W
             F(I,1) =  U(I,2)*W(I,3) - U(I,3)-W(I,2)
             F(I,2) = -U(I,1)*W(I,3) + U(I,3)*W(I,1)
             F(I,3) =  U(I,1)*W(I,2) - U(I,2)*W(I,1)
 10   CONTINUE

      ! Now Calculate Energy
      DO 20 I = 1,NT-1
             NUM = (W(I,1)*F(I+1,1)+W(I,2)*F(I+1,2)+W(I,3)*F(I+1,3)) &
             - (W(I+1,1)*F(I,1)+W(I+1,2)*F(I,2)+W(I+1,3)*F(I,3))
             DEN = (F(I,1)*F(I+1,1)+F(I,2)*F(I+1,2)+F(I,3)*F(I+1,3)) &
             + (W(I+1,1)*W(I,1)+W(I+1,2)*W(I,2)+W(I+1,3)*W(I,3))
             OMEGA = ATAN2(NUM,DEN)

             ETWIST = ETWIST + (TAU/2)*(OMEGA**2)
 20   CONTINUE
      
      RETURN
      END
      
!---------------------------------------------------------------*
