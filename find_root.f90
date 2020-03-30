program find_root
#this program using 4 methods to find roots of function. bisection, false position,  newton-raphson, and secant


implicit none
integer,parameter::dprec=kind(1.0d0)
real(dprec)::x,z,a,b,l
integer::input
integer::i

print*,"input function"
print*,"f(x)=cos(x)-x" !change this function with function that you want to find its root
print*,"select the methods"
write(*,*)"1 for bisection"
write(*,*)"2 for false position"
write(*,*)"3 for newton raphson"
write(*,*)"4 for secant"

read(*,*) input

select case(input)
        case(1)

            do
                write(*,*) 'submit 2 number for guess-root'
                read(*,*) a,b
                    if(f(a)*f(b)<=0.) then
                    exit
                    end if
            end do

        call bisection(z,a,b)
            
        case(2)

            do
                write(*,*) 'submit 2 number for guess-root'
                read(*,*) a,b
                    if(f(a)*f(b)<=0.) then
                    exit
                    end if
            end do
        
        call falseposition(z,a,b)

        case(3)

            write(*,*) 'submit 1 number for guess-root'
            read(*,*) a
            call newtonraphson(z,l,a)

        case(4)

            do
                write(*,*) 'submit 2 number for guess-root'
                read(*,*) a,b
                    if(f(a)*f(b)<=0.) then
                    exit
                    end if
            end do
        
        call secant(z,a,b)
        

end select

contains

function f(x) result (z)
!the function
implicit none
integer,parameter::dprec=kind(1.0d0)
real(dprec),intent(in)::x
real(dprec)::z

z= cos(x)-x !change this function with function that you want to find its root
return

end function

function fs(j) result (l)
!the derivative of the function

implicit none
integer,parameter::dprec=kind(1.0d0)
real(dprec),intent(in)::j
real(dprec)::l

l=-sin(j)-1 !change this with the derivative function that you want to find its root

return
end function


subroutine bisection(z,a,b) 
!kesrel is the relative error of the roots
implicit none
integer,parameter:: dprec=kind(1.0d0)
integer::i
real(dprec),parameter::minimum=1.0e-10_dprec
real(dprec)::kesrel,c,akar,h
real(dprec),intent(inout)::a,b,z

c=0.5_dprec*(a+b)

do i=1,1000
    if(f(a)*f(c)<=0.)then
        b=c 
    else
        a=c 
    end if

    akar=0.5_dprec*(a+b)
    kesrel=abs((c-akar)/(akar))

    if(kesrel-minimum<=0.)exit
    c=akar
end do

print*,"the root is" ,akar
h=f(akar)
print*, "the value of function on its root", h
print*,"the iteration",i

return
end subroutine



subroutine secant(z,a,b) 

implicit none
integer,parameter:: dprec=kind(1.0d0)
integer::i,j
real(dprec),parameter::minimum=1.0e-5_dprec
real(dprec)::kesrel,c,h
real(dprec),intent(inout)::a,b,z

do
write(*,*) 'masukkan input'
read(*,*) a,b
    if(f(a)*f(b)<=0.) then
    exit
end if
end do

do j=1,2
    c=((a*f(b)-b*f(a))/(f(b)-f(a)))
    a=c
end do

do i=1,1000
    c=a-((a-b)/(f(a)-f(b)))*f(a)

    kesrel=abs((a-c)/c)

    if(kesrel-minimum<=0.)then
    exit
    else
    a=c
    end if
end do

print*,"the root is" ,c
h=f(c)
print*, "the value of function on its root", h
print*,"the iteration",i

return
end subroutine


subroutine falseposition(a,b,z)

implicit none
integer,parameter:: dprec=kind(1.0d0)
integer::i
real(dprec),parameter::minimum=1.0e-5_dprec
real(dprec)::kesrel,c,h
real(dprec),intent(inout)::a,b,z

do
write(*,*) 'masukkan input'
read(*,*) a,b
    if(f(a)*f(b)<=0.) then
    exit
end if
end do

do i=1,1000
    c=((a*f(b)-b*f(a))/(f(b)-f(a)))
    
    kesrel=abs((a-c)/c)

    if(kesrel-minimum<=0.)then
    exit
    else
    a=c
    end if

end do

print*,"the root is" ,c
h=f(c)
print*, "the value of function on its root", h
print*,"the iteration",i

return
end subroutine

subroutine newtonraphson(z,l,a)

implicit none
integer,parameter:: dprec=kind(1.0d0)
integer::i
real(dprec),parameter::minimum=1.0e-5_dprec
real(dprec)::kesrel,c,h
real(dprec),intent(inout)::a,z,l
write(*,*) 'masukkan input'
read(*,*) a 

do i=1,1000
    c=a-(f(a)/(fs(a)))

    kesrel=abs((a-c)/c)

    if(kesrel-minimum<=0.)then
    exit
    else
    a=c
    end if

end do

print*,"the root is" ,c
h=f(c)
print*, "the value of function on its root", f(c)
print*,"the iteration",i

return
end subroutine

end program
