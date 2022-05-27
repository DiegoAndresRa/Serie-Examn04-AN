program euler ! para una ec con valor inicial
    implicit none
    real(kind=8) t0,tf,y0,dt,ti,yi,yimas1,f
    integer n,i

    open(unit=45,file='euler.txt',status='replace',err=500)
    !*+*** DATOS DEL PROBLEMA *****
    t0=0.d0 !tiempo inicial 
    y0=0.d0 !y(t0)

    tf=1.d0 !tiempo final


    !****** DATOS DEL INCREMENTO ******
    n=100 ! número de divisiones del dóminio 
    dt=(tf-t0)/dble(n)
    33 format(I4,3x,f17.8,3x,f17.8,3x)
    do i=0,n
        ti=t0+dble(i)*dt
        yi=y0
        write(*,33) i,ti,yi
        write(45,33) i,ti,yi
        !calculo del yimas1
        yimas1=yi+f(yi,ti)*dt
        y0=yimas1
    end do
    close(unit=45,status='keep',err=500)    
500 end program euler

function f(y,t)
    implicit none
    real(kind=8) f,y,t
    f=y**2+1.d0
end function