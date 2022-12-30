SUBROUTINE MPDATA2D(U1,U2,X,H,N,M,IORD,ISOR,NONOS,IDIV,IFL)

! THIS SUBROUTINE SOLVES 2-D ADVECTIVE TRANSPORT IN CARTESIAN GEOMETRY
! ON STAGGERRED GRID (X,Y VELOCITIES SHIFTED HALF GRID IN X, Y DIR, RESP)
!*************************************************************************
! ADVECTION ALGORITHM: IORD - NUMBER OF ITERATIONS (IORD=1 OVERWRITES
! CALLS TO OTHER OPTIONS AND GIVES SIMPLE UPSTREAM SCHEME); ISOR=1 2ND ORDER
! COMPUTATIONS WHEREAS ISOR=3 AND IORD=3 3D ORDER SCHEME; IDIV=1 ACTIVATES
! CORRECTION FOR DIVERGENT FLOW; NONOS=1 STRICTLY MONOTONE ADVECTION
!   N O T E: idiv MUST be 0 for a nondivergent flow
! A GOOD POINT TO START WOULD BE:
!      PARAMETER(IORD0=2,ISOR0=1,IDIV0=0,NONO=0)
! IFL IS THE FLAG TO DISTINGUISH BETWEEN FIELDS THAT ARE BEING ADVECTED;
! THAT IS NECESSARY TO PROVIDE LATERAL BOUNDARY CONDITIONS ON INLFOW
! AND LOWER BC FOR SEDIMENTING FIELDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ** NOTE THAT THIS ROUTINE WILL WORK FOR FIELD WITH VARIABLE SIGN **
! ** (AS MOMENTUM) SINCE ABSOLUTE VALUES ARE USE IN THE DEFINITION **
! **             OF ANTIDIFFUSIVE VELOCITIES                       **
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! AUTHOR: Piotr Smolarkiewicz (smolar@ncar.ucar.edu), (303)-497-8972
!  modified for this test by Wojciech Grabowski (grabow@ncar.ucr.edu)
!*************************************************************************
!      PARAMETER(N1=251,N2=126) !N1 = NX+1 / N2 = NZ + 1
!      PARAMETER(N1M=N1-1,N2M=N2-1)
!      DIMENSION U1(N+1,M),U2(N,M+1),X(N,M),H(N,M)
!      DIMENSION V1(N1,N2M),V2(N1M,N2),F1(N1,N2M),F2(N1M,N2)
!     *         ,CP(N1M,N2M),CN(N1M,N2M)
!      REAL MX(N1M,N2M),MN(N1M,N2M)

DATA EP/1.E-10/

integer :: N1,N2,N1M,N2M
real :: U1(N+1,M),U2(N,M+1),X(N,M),H(N,M)
real, allocatable :: V1(:,:),V2(:,:),F1(:,:),F2(:,:)
real, allocatable :: CP(:,:),CN(:,:),MX(:,:),MN(:,:)

! for use NOT on a CRAY computer you have to replace CVMGM
! with a substitute which gives the same result: below is
! an example of something that will work:
! CVMGM(a,b,c)= a if c.lt.0 or b if c.ge.0

aneg(d)=.5*(1.-sign(1.,d))
apos(d)=.5*(1.+sign(1.,d))
cvmgm(a,b,c)=aneg(c)*a + apos(c)*b

DONOR(Y1,Y2,A)=CVMGM(Y2,Y1,A)*A
VDYF(X1,X2,A,R)=(ABS(A)-A**2/R)*(ABS(X2)-ABS(X1))/(ABS(X2)+ABS(X1)+EP)
VCORR(A,B,Y1,Y2,R)=-0.125*A*B*Y1/(Y2*R)
VCOR31(A,X0,X1,X2,X3,R)= -(A -3.*ABS(A)*A/R+2.*A**3/R**2)/3.*(ABS(X0)+ABS(X3)-ABS(X1)-ABS(X2))/(ABS(X0)+ABS(X3)+ABS(X1)+ABS(X2)+EP)
VCOR32(A,B,Y1,Y2,R)=0.25*B/R*(ABS(A)-2.*A**2/R)*Y1/Y2
VDIV1(A1,A2,A3,R)=0.25*A2*(A3-A1)/R
VDIV2(A,B1,B2,B3,B4,R)=0.25*A*(B1+B2-B3-B4)/R
PP(Y)= AMAX1(0.,Y)
PN(Y)=-AMIN1(0.,Y)

N1 = N+1
N2 = M+1
N1M= N1-1
N2M= N2-1
allocate(V1(N1,N2M),V2(N1M,N2),F1(N1,N2M),F2(N1M,N2))
allocate(CP(N1M,N2M),CN(N1M,N2M),MX(N1M,N2M),MN(N1M,N2M))

! test dimensions:
if(n1m.ne.n.or.n2m.ne.m) then
  print*,' dimensions do not match in advection. stop'
  stop 'mpdata'
end if

IF(ISOR.EQ.3) IORD=MAX0(IORD,3)

V1 = U1
V2 = U2

IF(NONOS.EQ.1) THEN
  do J=2,N2-2
    do I=2,N1-2
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
      MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
    end do
  end do
END IF

do K=1,IORD
  do J=1,N2-1
    do I=2,N1-1
      F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
    end do
    F1( 1,J)=DONOR(X(N1-1,J),X(1,J),V1(1,J))
    !F1(N1,J)=DONOR(X(N1-1,J),X(1,J),V1(N1,J))
    F1(N1,J)=F1( 1,J)
    !if (k.eq.1) then
    !  F1( 1,J)=DONOR(X(N1-1,J),X(1,J),V1(1,J))
    !  F1(N1,J)=DONOR(X(N1-1,J),X(1,J),V1(N1,J))
    !else
    !  F1( 1,J)=0.
    !  F1(N1,J)=0.
    !end if
  end do

  do I=1,N1-1
    do J=2,N2-1
      F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
    end do
    if (k.eq.1.and.ifl.eq.4) then
      F2(I, 1)=DONOR(0.,X(I,1),V2(I,1))  ! fallout of rain
      F2(I,N2)=0.                        ! no fall-in of rain
    else
      F2(I, 1)=0.
      F2(I,N2)=0.
    end if
  end do

  do J=1,N2-1
    do I=1,N1-1
      X(I,J)=X(I,J)-(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
    end do
  end do
   
  IF (K.EQ.IORD) THEN
    RETURN
  END IF

  F1 = V1
  V1 = 0.
      
  F2 = V2
  V2 = 0
      
  DO J=2,N2-2
    DO I=2,N1-1
      V1(I,J)=VDYF(X(I-1,J),X(I,J),V1(I,J),.5*(H(I-1,J)+H(I,J)))      &
      +VCORR(V1(I,J), F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),        &
      ABS(X(I-1,J+1))+ABS(X(I,J+1))-ABS(X(I-1,J-1))-ABS(X(I,J-1)),    &
      ABS(X(I-1,J+1))+ABS(X(I,J+1))+ABS(X(I-1,J-1))+ABS(X(I,J-1))+EP, &
      .5*(H(I-1,J)+H(I,J)))
    END DO
  END DO

  IF(IDIV.EQ.1) THEN
    do J=2,N2-2
      do I=2,N1-1
        V1(I,J)=V1(I,J)                                          &
        -VDIV1(F1(I-1,J),F1(I,J),F1(I+1,J),.5*(H(I-1,J)+H(I,J))) &
        -VDIV2(F1(I,J),F2(I-1,J+1),F2(I,J+1),F2(I-1,J),F2(I,J),  &
        .5*(H(I-1,J)+H(I,J)))
      end do
    end do
  END IF

  do J=2,N2-1
    do I=2,N1-2
      V2(I,J)=VDYF(X(I,J-1),X(I,J),V2(I,J),.5*(H(I,J-1)+H(I,J)))      &
      +VCORR(V2(I,J), F1(I,J-1)+F1(I,J)+F1(I+1,J)+F1(I+1,J-1),        &
      ABS(X(I+1,J-1))+ABS(X(I+1,J))-ABS(X(I-1,J-1))-ABS(X(I-1,J)),    &
      ABS(X(I+1,J-1))+ABS(X(I+1,J))+ABS(X(I-1,J-1))+ABS(X(I-1,J))+EP, &
      .5*(H(I,J-1)+H(I,J)))
    end do
  end do

  IF(IDIV.EQ.1) THEN
    do J=2,N2-1
      do I=2,N1-2
        V2(I,J)=V2(I,J)                                           &
        -VDIV1(F2(I,J-1),F2(I,J),F2(I,J+1),.5*(H(I,J-1)+H(I,J)))  &
        -VDIV2(F2(I,J),F1(I+1,J),F1(I+1,J-1),F1(I,J-1),F1(I,J),   &
        .5*(H(I,J-1)+H(I,J)))
      end do
    end do
  END IF

  IF (ISOR.EQ.3) THEN
    do J=2,N2-2
      do I=3,N1-2
        V1(I,J) = V1(I,J) +  VCOR31(                                     &
        F1(I,J), X(I-2,J),X(I-1,J),X(I,J),X(I+1,J),.5*(H(I-1,J)+H(I,J)))
      end do
    end do

    do J=2,N2-2
      do I=3,N1-2
        V1(I,J) = V1(I,J)                                               &
        + VCOR32(F1(I,J),F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),       &
        ABS(X(I,J+1))-ABS(X(I,J-1))-ABS(X(I-1,J+1))+ABS(X(I-1,J-1)),    &
        ABS(X(I,J+1))+ABS(X(I,J-1))+ABS(X(I-1,J+1))+ABS(X(I-1,J-1))+EP, &
        .5*(H(I-1,J)+H(I,J)))
      end do
    end do

    do J=3,N2-2
      do I=2,N1-2
        V2(I,J) = V2(I,J) + VCOR31(                                     &
        F2(I,J),X(I,J-2),X(I,J-1),X(I,J),X(I,J+1),.5*(H(I,J-1)+H(I,J)))
      end do
    end do

     do J=3,N2-2
      do I=2,N1-2
        V2(I,J) = V2(I,J) + VCOR32(                                     &
        F2(I,J),F1(I,J-1)+F1(I+1,J-1)+F1(I+1,J)+F1(I,J),                &
        ABS(X(I+1,J))-ABS(X(I-1,J))-ABS(X(I+1,J-1))+ABS(X(I-1,J-1)),    &
        ABS(X(I+1,J))+ABS(X(I-1,J))+ABS(X(I+1,J-1))+ABS(X(I-1,J-1))+EP, &
        .5*(H(I,J-1)+H(I,J)))
      end do
      end do
  END IF

  IF (NONOS.NE.0) THEN
  ! NON-OSSCILATORY OPTION
    do J=2,N2-2
      do I=2,N1-2
        MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MX(I,J))
        MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MN(I,J))
      end do
    end do

    do J=2,N2-2
      do  I=2,N1-1
        F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
      end do
    end do

    do J=2,N2-1
      do I=2,N1-2
        F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
      end do
    end do

    do J=2,N2-2
      do I=2,N1-2
        CP(I,J) = (MX(I,J) - X(I,J))*H(I,J)/                       &
        (PN(F1(I+1,J))+PP(F1(I,J))+PN(F2(I,J+1))+PP(F2(I,J))+EP)
        CN(I,J)=(X(I,J) - MN(I,J))*H(I,J)/                         &
        (PP(F1(I+1,J))+PN(F1(I,J))+PP(F2(I,J+1))+PN(F2(I,J))+EP)
      end do
    end do
  
    do J=3,N2-2
      do I=3,N1-2
        V1(I,J) = PP(V1(I,J))*                                 &
        ( AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1., X(I-1,J)))   &
        + AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1.,-X(I-1,J))) ) &
        - PN(V1(I,J))*                                         &
        ( AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1., X(I ,J )))   &
        + AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1.,-X(I ,J ))) )

        V2(I,J) = PP(V2(I,J))*                                 &
        ( AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1., X(I,J-1)))   &
        + AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1.,-X(I,J-1))) ) &
        - PN(V2(I,J))*                                         &
        ( AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1., X(I ,J )))   &
        + AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1.,-X(I ,J ))) )
      end do
    end do
  END IF
end do
RETURN
END SUBROUTINE
