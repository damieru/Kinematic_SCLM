!Subroutines for position and velocity fluctuation correction

SUBROUTINE POSITION_CORRECTION (POS)
  
   !The positions are corrected to force consistency
  
   use SD_VARIABLES, only: N_sd,X
   use SD_FUNCTIONS, only: PARTICLE_BC
   USE GRID        , only: NX, NZ, DX, DZ, LX, LZ
   use OMP_LIB
  
   IMPLICIT NONE
  
   REAL*8, intent(inout) :: POS(N_sd,2)
  
   REAL*8 :: P(NX,NZ)
   REAL*8 :: DPX(NX+1,NZ),DPZ(NX,NZ+1) !Gradient at the cell faces
   REAL*8 :: DPX_K,DPZ_K
 
   REAL*8 :: X0,Z0,C,D !Variables for linear interpolation
  
   INTEGER :: IX,IZ, k
 
   !Poisson problem for mean particle kinematic pressure (position correction)

   CALL EVAL_KINEMATIC_PARTICLE_PRESSURE_POS (POS, P, DPX, DPZ )
   
   !$OMP parallel do private(IX, IZ, X0, Z0, C, D, DPX_K, DPZ_K)
   do k = 1,N_sd
      !Linear interpolation of the pressure gradient at the particle
      !position - the pressure gradient is stored at the cell faces
      IX = INT(POS(k,1)/DX) + 1
      IZ = INT(POS(k,2)/DZ) + 1
      
      !NOTA: Alterar para usar as posições dos nós (já calculadas) - vale a pena?
      X0 = DX*(IX-1)
      Z0 = DZ*(IZ-1)

      C = (POS(k,1)-X0)/DX
      D = (POS(k,2)-Z0)/DZ
      
      DPX_K = (1.D0-C)*DPX(IX,IZ) + C*DPX(IX+1,IZ  )
      DPZ_K = (1.D0-D)*DPZ(IX,IZ) + D*DPZ(IX  ,IZ+1)
      
      !Correction
      POS(k,:) = PARTICLE_BC(POS(k,:) - [DPX_K, DPZ_K], LX, LZ)

      !NOTAS PARA AMANHÃ 21.03: 
      ! - Paralelizar este loop (private DPX_K e DPZ_K) - DONE :)
      ! - Verificar se P pode ser opcional (intent out)
      
      !NO REFLECTION of the vertical velocity fluctuation when particle
      !hit no-flux boundaries     
   end do
   !$OMP end parallel do
END SUBROUTINE POSITION_CORRECTION

SUBROUTINE EVAL_KINEMATIC_PARTICLE_PRESSURE_POS (XP,P,PX,PZ)

   USE GRID
   USE STEPPER

   IMPLICIT NONE

   REAL*8 :: XP(2,NN)
   REAL*8 :: P(NX,NZ),PX(NXP,NZ),PZ(NX,NZP)

   REAL*8 :: DENS_0,DENS_LOCAL

   REAL*8 :: S(0:NX-1,1:NZ),F(0:NX-1,1:NZ)
   REAL*8 :: GB(0:NX-1),GT(0:NX-1)

   INTEGER :: H,I,J
   INTEGER :: IX,IZ

   real*8  :: lapl,s1,s2

   H(IX,IZ) = (IZ-1)*NX + IX

   !*************************************************************
   !Right-hand side (RHS) of the Poisson equation for "P"
   !-------------------------------------------------------------

   DENS_0 = REAL(NN) / REAL(NX*NZ)

   CALL UPDATE_LINKED_LIST_P ( XP )

   DO IX = 1,NX
      DO IZ = 1,NZ
         I = H(IX,IZ)
         DENS_LOCAL = REAL( POINT_P(I+1) - POINT_P(I) ) !Local density
         F(IX-1,IZ) = 1.D0 - DENS_LOCAL/DENS_0 !Right-hand side (source)
      END DO
   END DO

   !Bottom and top Neumann boundary conditions 
   DO IX = 1,NX
      GB(IX-1) = 0.D0 
      GT(IX-1) = 0.D0 
   END DO

   !********************************************************************
   !Poisson solver
   call poisson_2d_mixed_bc_xper_yneu_stagg (nx,nz,dx,dz,f,gb,gt,s)

   DO IX = 1,NX
      DO IZ = 1,NZ
         P(IX,IZ) = S(IX-1,IZ)
      END DO
   END DO
   
   !********************************************************************
   !Pressure gradient in the staggered grid

   !1) Vertical gradient
   DO IX = 1,NX
      !Interior points
      DO IZ = 2,NZ
         PZ(IX,IZ) = ( P(IX,IZ) - P(IX,IZ-1) ) / DZ
      END DO
      
      !Boundary points
      PZ(IX,  1) = GB(IX-1) !According to the Neumann b. cond. (bottom)
      PZ(IX,NZP) = GT(IX-1) !According to the Neumann b. cond. (top)
   END DO
   
   !2) Horizontal gradient
   DO IZ = 1,NZ
      !Interior points
      DO IX = 2,NX
         PX(IX,IZ) = ( P(IX,IZ) - P(IX-1,IZ) ) / DX
      END DO
      
      !Boundary points - periodicity
      PX(  1,IZ) = ( P(1,IZ) - P(NX,IZ) ) / DX
      PX(NXP,IZ) = PX(1,IZ)
   END DO
   
   !********************************************************************
   !check solution   
   ! s1 = 0.d0
   ! s2 = 0.d0
   ! do i = 1,nx-2
   !    do j = 2,nz-1
   !       lapl = ( s(i+1,j) + s(i-1,j) - 2.d0*s(i,j) ) / dx**2 &
   !          & + ( s(i,j+1) + s(i,j-1) - 2.d0*s(i,j) ) / dz**2      
   !       s1 = s1 + ( lapl - f(i,j) )**2
   !       s2 = s2 + f(i,j)**2
   !    end do
   ! end do

   ! >>> FOR DEBUGGING PURPOSES <<<
   !write(*,*)
   !write(*,*) "Check FFT-based Poisson solver - POSITIONS"
   !write(*,*) "Relative error: ", sqrt(s1/s2)
END SUBROUTINE EVAL_KINEMATIC_PARTICLE_PRESSURE_POS