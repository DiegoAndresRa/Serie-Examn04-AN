PROGRAM metodo_Euler_sistema2x2
    IMPLICIT NONE
    REAL (KIND=8) t0,tf,y0,z0,dt,f1,f2,ti,yi,zi,yimas1,zimas1
    INTEGER n,i

OPEN(UNIT=54,FILE='Euler 2x2.txt',STATUS='replace',ERR=500)
    !*****Datos del problema
    t0=0.d0  !tiempo inicial
    y0=10.d0  !y(t0)
    z0=60.d0  !z(t0)
    tf=50.d0  !tiempo final
    !******Datos referentes al incremento
    n=1000               !Número de nodos (o divisiones) del dominio [y0,yf]]
    dt=(tf-t0)/Dble(n)  !incremento

33 FORMAT(I4,3x,F17.8,3x,F17.8,3x,F17.8)
    DO i=0,n,1
        ti=t0+DBLE(i)*dt
        yi=y0
        zi=z0
        WRITE(54,33) i,ti,yi,zi
        !Calculo de yimas1
        yimas1=yi+f1(yi,zi,ti)*dt
        zimas1=zi+f2(yi,zi,ti)*dt

        y0=yimas1
        z0=zimas1
    END DO
    WRITE(*,'(A,17x,A2,18x,A2,18x,A2)') 'i','ti','yi','zi'
    WRITE(*,33) i,ti,yi,zi

CLOSE(UNIT=54,STATUS='keep',ERR=500)

500 END PROGRAM

!*********f1(y,z,t)
FUNCTION f1(y,z,t)
    IMPLICIT NONE
    REAL (KIND=8) f1,y,z,t,L,C,R
    f1=(-0.16d0*y)+(0.012d0*y*z)
END FUNCTION
!*********f2(y,z,t)
FUNCTION f2(y,z,t)
    IMPLICIT NONE
    REAL (KIND=8) f2,y,z,t
    f2=(0.9d0*z)-(0.05d0*y*z)
END FUNCTION