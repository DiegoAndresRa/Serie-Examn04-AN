PROGRAM metodo_RK4 !Método de Runge-Kutta para una EDO con valor inicial de
    IMPLICIT NONE    !cuarto orden
    REAL (KIND=8) t0,tf,y0,dt,f,ti,yi,yimas1,k1,k2,k3,k4
    INTEGER n,i

OPEN(UNIT=54,FILE='RK4.txt',STATUS='replace',ERR=500)
    !*****Datos del problema
    t0=0.d0  !tiempo inicial
    y0=0.d0  !y(t0)
    tf=30.d0  !tiempo final
    !******Datos referentes al incremento
    n=60               !Número de nodos (o divisiones) del dominio [y0,yf]]
    dt=(tf-t0)/Dble(n)  !incremento

WRITE(54,'(A1,15x,A2,20x,A5)') 'i','ti','y(ti)'

33 FORMAT(I6,3x,F17.8,3x,F17.8)
    DO i=0,n,1
        ti=t0+DBLE(i)*dt
        yi=y0
       !WRITE(*,33) i,ti,yi
        WRITE(54,33) i,ti,yi
        !Términos de RK4
        k1=dt*f(yi,ti)
        k2=dt*f(yi+0.5d0*k1,ti+dt*0.5d0)
        k3=dt*f(yi+0.5d0*k2,ti+dt*0.5d0)
        k4=dt*f(yi+k3,ti+dt)
        !Cálculo del la función un punto adelante
        yimas1=yi+(1.d0/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
        y0=yimas1
    END DO
    WRITE(*,*) 'y(',tf,')=',yi


CLOSE(UNIT=54,STATUS='keep',ERR=500)

500 END PROGRAM

!*********f(y,t)
function f(v,t)
    implicit none
    real(kind=8) f,v,t
    f=(668.061d0-0.22d0*(v**2))/68.1
end function