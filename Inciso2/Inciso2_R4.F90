PROGRAM metodo_RK4 !Método de Runge-Kutta para una EDO con valor inicial de
    IMPLICIT NONE    !cuarto orden
    REAL (KIND=8) t0,tf,y0,dx,f,ti,yi,yimas1,k1,k2,k3,k4
    INTEGER n,i

OPEN(UNIT=54,FILE='RK4.txt',STATUS='replace',ERR=500)
    !*****Datos del problema
    t0=0.d0  !tiempo inicial
    y0=2.d0  !y(t0)
    tf=1.d0  !tiempo final
    !******Datos referentes al incremento
    n=100                !Número de nodos (o divisiones) del dominio [y0,yf]]
    dx=(tf-t0)/Dble(n)  !incremento

WRITE(54,'(A1,15x,A2,20x,A5)') 'i','ti','y(ti)'

33 FORMAT(I6,3x,F17.8,3x,F17.8)
    DO i=0,n,1
        ti=t0+DBLE(i)*dx
        yi=y0
       !WRITE(*,33) i,ti,yi
        WRITE(54,33) i,ti,yi
        !Términos de RK4
        k1=dx*f(yi,ti)
        k2=dx*f(yi+0.5d0*k1,ti+dx*0.5d0)
        k3=dx*f(yi+0.5d0*k2,ti+dx*0.5d0)
        k4=dx*f(yi+k3,ti+dx)
        !Cálculo del la función un punto adelante
        yimas1=yi+(1.d0/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
        y0=yimas1
    END DO
    WRITE(*,*) 'y(',tf,')=',yi


CLOSE(UNIT=54,STATUS='keep',ERR=500)

500 END PROGRAM

!*********f(y,t)
function f(y,x)
    implicit none
    real(kind=8) f,y,x
    f=-y+x+2.D0
end function