SUBROUTINE GET_SCALAR_FORCES_TURB
  
   !Get the "forces" (i.e., the RHS of THE Eulerian (scalar) advection
   !equations) accounting for turbulent fluxes
  
   USE ENVIRONMENT
   USE THERMODYNAMIC_VAR
   USE STEPPER 
   USE GRID

   IMPLICIT NONE

   REAL*8, ALLOCATABLE ::  JX(:,:), JZ(:,:), DIV(:,:)

   IF ( .NOT.ALLOCATED(FTH_TURB_DT) ) THEN
      ALLOCATE ( FTH_TURB_DT(NX,NZ), FRV_TURB_DT(NX,NZ) )
   END IF
  
   ALLOCATE ( JX(NX,NZ), JZ(NX,NZ), DIV(NX,NZ) )

   CALL GET_TURBULENT_SGS_PARAMETERS !New version (includes all relevant param)

   !-------------------------------------------------------------------------------
   !Turbulent flux of potential temperature THETA(:,:)
  
   CALL EVAL_GRADIENT_SCALAR ( THETA, JX, JZ )
   
   JX = EDDY_DIFF*JX 
   JZ = EDDY_DIFF*JZ
   
   JZ(1:NX, 1) = 0.D0  !ZERO vertical turbulent flux boundary condition
   JZ(1:NX,NZ) = 0.D0  !ZERO vertical turbulent flux boundary condition
   
   CALL EVAL_DIVERGENCE ( DIV, JX, JZ )
   
   FTH_TURB_DT = DIV*DT
   
   !-------------------------------------------------------------------------------
   !Turbulent flux of vapor mixing ratio QV(:,:)
  
   CALL EVAL_GRADIENT_SCALAR ( RV, JX, JZ )
  
   JX = EDDY_DIFF*JX
   JZ = EDDY_DIFF*JZ
  
   JZ(1:NX, 1) = 0.D0  !ZERO vertical turbulent flux boundary condition
   JZ(1:NX,NZ) = 0.D0  !ZERO vertical turbulent flux boundary condition
   
   CALL EVAL_DIVERGENCE ( DIV, JX, JZ )
  
   FRV_TURB_DT = DIV*DT

   DEALLOCATE (JX,JZ,DIV)
  
END SUBROUTINE GET_SCALAR_FORCES_TURB

SUBROUTINE EVAL_GRADIENT_SCALAR (S,DSX,DSZ)
  
   !Evaluates the gradient of the scalar "S" at every point of the
   !primary grid. Periodic boundary conditions in the horizontal (x)
   !direction and no-flux BC in the vertical (z) direction.

   USE GRID
  
   IMPLICIT NONE
  
   REAL*8    :: S(NX,NZ), DSX(NX,NZ), DSZ(NX,NZ)
   INTEGER*8 :: IX, IZ, PER

   DO IX = 1,NX
      DO IZ = 1,NZ
         !Horizontal gradient (periodic boundary conditions)
         DSX(IX,IZ) = 0.5D0/DX*( S(PER(IX+1,NX),IZ) - S(PER(IX-1,NX),IZ) )
        
         !Vertical gradient
         IF (IZ.EQ.1) THEN
            DSZ(IX,IZ) = 0.5D0/DZ*( -3.D0*S(IX,IZ) + 4.D0*S(IX,IZ+1) - S(IX,IZ+2) )
         ELSE IF (IZ.EQ.NZ) THEN
            DSZ(IX,IZ) = 0.5D0/DZ*(  3.D0*S(IX,IZ) - 4.D0*S(IX,IZ-1) + S(IX,IZ-2) )   
         ELSE
            !Interior points
            DSZ(IX,IZ) = 0.5D0/DZ*( S(IX,IZ+1) - S(IX,IZ-1) )
         END IF
      END DO
   END DO
END SUBROUTINE EVAL_GRADIENT_SCALAR

SUBROUTINE EVAL_DIVERGENCE (DIV,JX,JZ)
  
   !Evaluates the divergence of the vector (JX,JZ) at every point of
   !the primary grid. Periodic boundary conditions in the horizontal
   !(x) direction and no-flux in the vertical (z) direction.
  
   USE GRID
  
   IMPLICIT NONE

   REAL*8    :: DIV(NX,NZ), JX(NX,NZ), JZ(NX,NZ)
   REAL*8    :: AX, AZ
   INTEGER*8 :: IX, IZ, PER
  
   DO IX = 1,NX
      DO IZ = 1,NZ
         !Horizontal contribution (periodic boundary conditions)
         AX = 0.5D0/DX*( JX(PER(IX+1,NX),IZ) - JX(PER(IX-1,NX),IZ) )
        
         !Vertical contribution
         IF (IZ.EQ.1) THEN
            AZ = 0.5D0/DZ*( -3.D0*JZ(IX,IZ) + 4.D0*JZ(IX,IZ+1) - JZ(IX,IZ+2) )
         ELSE IF (IZ.EQ.NZ) THEN
            AZ = 0.5D0/DZ*(  3.D0*JZ(IX,IZ) - 4.D0*JZ(IX,IZ-1) + JZ(IX,IZ-2) )
         ELSE
            !Interior points
            AZ = 0.5D0/DZ*( JZ(IX,IZ+1) - JZ(IX,IZ-1) )
         END IF

         DIV(IX,IZ) = AX + AZ
      END DO
   END DO
END SUBROUTINE EVAL_DIVERGENCE

FUNCTION PER(I,N)
  
   IMPLICIT NONE
   INTEGER*8 :: PER, I, N
  
   IF (I.LT.1) THEN
      PER = I + N
   ELSEIF (I.GT.N) THEN
      PER = I - N
   ELSE
      PER = I
   END IF  
END FUNCTION PER

