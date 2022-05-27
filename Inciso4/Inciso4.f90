PROGRAM metodo_Euler_sistema2x2
    IMPLICIT NONE
    REAL (KIND=8) t0,tf,I0,Q0,dt,f1,f2,ti,Ii,Qi,Iimas1,Qimas1
    INTEGER n,i

OPEN(UNIT=54,FILE='Euler 2x2.txt',STATUS='replace',ERR=500)
    !*****Datos del problema
    t0=0.d0  !tiempo inicial
    I0=0.d0  !I(t0)
    Q0=0.d0  !Q(t0)
    tf=5.d0  !tiempo final
    !******Datos referentes al incremento
    n=1000              !Número de nodos (o divisiones) del dominio [I0,If]]
    dt=(tf-t0)/Dble(n)  !incremento

33 FORMAT(I4,3x,F17.8,3x,F17.8,3x,F17.8)
    DO i=0,n,1
        ti=t0+DBLE(i)*dt
        Ii=I0
        Qi=Q0
        
        WRITE(54,33) i,ti,Ii,Qi
        !Calculo de Iimas1
        Iimas1=Ii+f1(Ii,Qi,ti)*dt
        Qimas1=Qi+f2(Ii,Qi,ti)*dt

        I0=Iimas1
        Q0=Qimas1
    END DO
    WRITE(*,'(A,17x,A2,18x,A2,18x,A2)') 'i','ti','Ii','Qi'
    WRITE(*,33) i,ti,Ii,Qi

CLOSE(UNIT=54,STATUS='keep',ERR=500)

500 END PROGRAM

!*********f1(I,Q,t)
FUNCTION f1(I,Q,t)
    IMPLICIT NONE
    REAL (KIND=8) f1,I,Q,t,L,C,R,E
    L=4.d0
    R=300.d0
    C=10.d0**(-3.d0)
    E=(127.d0*dsin(60.d0*t-(dacos(-1.d0)/2.d0)))/(1.d0+t**(2.d0))
    f1=-(1.d0/L)*((I*R)+(Q/C)-E)  
END FUNCTION
!*********f2(I,Q,t)
FUNCTION f2(I,Q,t)
    IMPLICIT NONE
    REAL (KIND=8) f2,I,Q,t
    f2=I
END FUNCTION