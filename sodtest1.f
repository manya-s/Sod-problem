      program sod
      implicit none
      integer :: i,n
      real :: gamma, ci 
      real :: x, dt, a
      real :: pl, rhol, ul, pr, rhor, ur
      real, dimension(1000,500) :: q0,q1,q4, p, rho, u, c, g,h
      real, dimension(1000,500) :: lambda1, lambda2, lambda3,delq0,delq1,delq4
      real, dimension(1000,500) :: f1,f2,f3,ua,ha,ca,qc1,qc2,qc3,fa1,fa2,fa3
      real, dimension(1000,500) :: e11,e12,e13,e21,e22,e23,e31,e32,e33

      ! q0 = rho*etot , q1 = rho*u q4= rho

      gamma = 7.0/5.0
      pl = 1
      rhol = 1
      ul = 0

      pr = 0.1
      rhor= 0.125
      ur = 0

      dt = 0.15

      do i = 1, 500
      	q0(i,1) = pl/(gamma - 1) + 0.5*rhol*(ul**2)
      	q1(i,1) = rhol*ul
      	q4(i,1) = rhol
      end do

      do i = 500,1000
      	q0(i,1) = pr/(gamma - 1) + 0.5*rhor*(ur**2)
      	q1(i,1) = rhor*ur
      	q4(i,1) = rhor
      end do


      
      open (1, file = "initial_sod.dat")
      do i = 1,1000
      p(i,1) = (gamma -1)*(q0(i,1) - 0.5*(q1(i,1)**2)/q4(i,1)) 
      rho(i,1) = q4(i,1)
      u(i,1) = q1(i,1)/q4(i,1)
      x = i
      write(1,*)  x, p(i,1), rho(i,1), u(i,1) 
      end do


      do n = 1,500
      	q0(1,n) = 2.5
      	q1(1,n) = 0
      	q4(1,n) = 1
      	h(1,n) = gamma*((q0(1,n)/q4(1,n)) + 0.5*(q1(1,n)**2)/(q4(1,n)**2)) + 0.5*(q1(1,n)**2)/q4(1,n)
      	p(1,n) = (gamma -1)*(q0(1,n) - 0.5*(q1(1,n)**2)/q4(1,n)) 
        rho(1,n) = q4(1,n)
        u(1,n) = q1(1,n)/q4(1,n)
        c(1,n) = (gamma*p(1,n)/rho(1,n))**0.5
      end do


      !time step begins


      do n = 1,499
      	do i = 2,1000

      		h(i,n) = gamma*((q0(i,n)/q4(i,n)) + 0.5*(q1(i,n)**2)/(q4(i,n)**2)) + 0.5*(q1(i,n)**2)/q4(i,n)
      		p(i,n) = (gamma -1)*(q0(i,n) - 0.5*(q1(i,n)**2)/q4(i,n)) 
            rho(i,n) = q4(i,n)
            u(i,n) = q1(i,n)/q4(i,n)
            c(i,n) = (gamma*p(i,n)/rho(i,n))**0.5

            f1(i,n) = rho(i,n)*u(i,n)*h(i,n)
            f2(i,n) = rho(i,n)*(u(i,n)**2) + p(i,n)
            f3(i,n) = rho(i,n)*u(i,n)

            ua(i,n) = ((q4(i,n)**0.5)*u(i,n) + (q4(i+1,n)**0.5)*u(i+1,n))/((q4(i,n)**0.5) + (q4(i+1,n)**0.5))
            ha(i,n) = ((q4(i,n)**0.5)*h(i,n) + (q4(i+1,n)**0.5)*h(i+1,n))/((q4(i,n)**0.5) + (q4(i+1,n)**0.5))
            ca(i,n) = (gamma-1)*(ha(i,n) - 0.5*(ua(i,n)**2))

            lambda1(i,n) = ua(i,n) - ca(i,n)
            lambda2(i,n) = ua(i,n) + ca(i,n)
            lambda3(i,n) = ua(i,n)

            e11(i,n) = ha(i,n) - ca(i,n)*ua(i,n)
            e12(i,n) = ua(i,n) - ca(i,n)
            e13(i,n) = 1

            e21(i,n) = ha(i,n) + ca(i,n)*ua(i,n)
            e22(i,n) = ua(i,n) + ca(i,n)
            e23(i,n) = 1

            e31(i,n) = 0.5*(ua(i,n)**2)
            e32(i,n) = ua(i,n)
            e33(i,n) = 1

            delq0(i,n) = q0(i+1,n) - q0(i,n)
            delq1(i,n) = q1(i+1,n) - q1(i,n)
            delq4(i,n) = q4(i+1,n) - q4(i,n)

            qc1(i,n) = (0.5*(gamma-1)/(ca(i,n)**2))*(0.5*(ua(i,n)**2)*delq4(i,n) - u(i,n)*delq1(i,n) + delq0(i,n)) - 0.5*(delq1(i,n) - ua(i,n)*delq4(i,n))/ca(i,n)
            qc2(i,n) = (0.5*(gamma-1)/(ca(i,n)**2))*(0.5*(ua(i,n)**2)*delq4(i,n) - u(i,n)*delq1(i,n) + delq0(i,n)) + 0.5*(delq1(i,n) - ua(i,n)*delq4(i,n))/ca(i,n)
            qc3(i,n) = (0.5*(gamma-1)/(ca(i,n)**2))*((ha(i,n) - ua(i,n)**2)*delq4(i,n) + u(i,n)*delq1(i,n) - delq0(i,n))

            fa1(i,n) =  0.5*(f1(i,n)+f1(i+1,n)) - 0.5*(abs(lambda1(i,n))*qc1(i,n)*e11(i,n) + abs(lambda2(i,n))*qc2(i,n)*e21(i,n) + abs(lambda3(i,n))*qc3(i,n)*e31(i,n))
            fa2(i,n) =  0.5*(f2(i,n)+f2(i+1,n)) - 0.5*(abs(lambda1(i,n))*qc1(i,n)*e12(i,n) + abs(lambda2(i,n))*qc2(i,n)*e22(i,n) + abs(lambda3(i,n))*qc3(i,n)*e32(i,n))
            fa3(i,n) =  0.5*(f3(i,n)+f3(i+1,n)) - 0.5*(abs(lambda1(i,n))*qc1(i,n)*e13(i,n) + abs(lambda2(i,n))*qc2(i,n)*e23(i,n) + abs(lambda3(i,n))*qc3(i,n)*e33(i,n))

            q0(i,n+1) = q0(i,n) - dt*(fa1(i,n)-fa1(i-1,n))
            q1(i,n+1) = q1(i,n) - dt*(fa2(i,n)-fa2(i-1,n))
            q4(i,n+1) = q4(i,n) - dt*(fa3(i,n)-fa3(i-1,n))
     
                   
           
            !if (lambdayyy(i,n) .ge. lambda1(i+1,n).and.lambda1(i,n).ge. a) then
            !	a = lambda1(i,n)
            !	else 
            	!	a = a
            !end if          
        end do
        !dt = 0.5/a
      end do
        !print*, a

        open(3, file = "hello.dat")
        do i = 1,1000
        	x = i
        	write(3,*)  q4(i,400) , p(i,400), u(i,400), x
        end do
        close(3)

      end program sod

 



