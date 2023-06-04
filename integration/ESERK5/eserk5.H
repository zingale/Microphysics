      subroutine ESERK(neqn,t,tend,dt,g,f,tol,spcrad,iwork,work,idid)
c
c     Stabilized Runge-Kutta (composite Euler steps)
c     with best ordering
c
c     ================
c
c
c no implicit variables
c
      implicit none
c
c
c     passed variables
c
      integer          neqn, idid, iwork(10), spcrad(2)
      double precision t, tend, dt, tol, work(1+5*neqn), g(neqn)
      external         f
c
c     local variables
c
      integer          i, j, k, l, m, i1, sec
      integer          stage, stage_intern, max_stages
      integer          total, accepted, rejected, rej_counter, feval
      integer          start   ! startindex in array b for SERK

      double precision eigmax, rho ! largest eigenvalue
      double precision b(15484)                                       !!! new
      double precision sum(neqn), err, fac
      double precision g_calc(neqn), g_help(neqn)
      double precision var_1(neqn), var_2(neqn), g_save(neqn), sc(neqn)
      double precision g_oth(neqn), g_rej(neqn), r, pas, x0n, y0n(neqn)
      double precision dtmax, ans, hmax, tfail, al1
      double precision, dimension(neqn,2001) :: g_work                !!! check!!! 2001 (maxstage+1)
      double precision, dimension(neqn,5) ::   y                      !!! new, old 4
     
      err = 0.0d0
      max_stages = 0

      m=0
      total = 0
      accepted = 0
      rejected = 0
      rej_counter = 11
      feval = 0

      tfail = t

      call f(neqn, t, g,  g_calc)
      feval = feval + 1
      total = total + 1
      
      do i=1, neqn
          sum(i) = 0.0d0
          g_help(i) = g(i)
          g_save(i) = g_calc(i)
          g_rej(i) = g(i)
          var_1(i) = g_calc(i)
      end do

c     number of accepted steps - needed in spectralradius computation using nonlinear power method
      iwork(7) = 0
c     number of rejected steps - needed in spectralradius computation using nonlinear power method
      iwork(8) = 0
c     number of evaluations of F used to estimate the spectral radius - needed in spectralradius computation using nonlinear power method
      iwork(9) = 0
c     needed in spectralradius computation using nonlinear power method
      hmax = abs(tend - t)
c     spcrad(1) = 0  -> subroutine rho is providing an upper bound on the spectralradius
c               = 1  -> compute an upper bound on the spectralradius using a nonlinear power method
c     spcrad(2) = 0  ->  The Jacobian may not be constant.
c               = 1  ->  The Jacobian is constant.
      if(spcrad(1)==0) then
          eigmax=rho(neqn,t,g)
          write(*,*) 'own rho =', eigmax
      else
          call SERKrho(neqn,t,f,g_help,g_calc,iwork,hmax,work,eigmax,
     *                 idid)
          write(*,*) 'internal rho =', eigmax
          write(*,*) 'number of func evals', iwork(9)
      end if
      dtmax = 0.98d0*2000.0d0*2000.0d0/eigmax
      
c     set stages and startindexes
      ans = sqrt( eigmax * dt/0.98d0)

      if(ans > 1800) then
            stage = 2000
            start = 13483
            stage_intern = 200
      else if(ans > 1600) then
            stage = 1800
            start = 11682
            stage_intern = 200
      else if(ans > 1400) then
            stage = 1600
            start = 10081
            stage_intern = 200
      else if(ans > 1200) then
            stage = 1400
            start = 8680
            stage_intern = 200
      else if(ans > 1000) then
            stage = 1200
            start = 7479
            stage_intern = 200
      else if(ans > 900) then
            stage = 1000
            start = 6478
            stage_intern = 100
      else if(ans > 800) then
            stage = 900
            start = 5577
            stage_intern = 100
      else if(ans > 700) then
            stage = 800
            start = 4776
            stage_intern = 100
      else if(ans > 600) then
            stage = 700
            start = 4075
            stage_intern = 100
      else if(ans > 500) then
            stage = 600
            start = 3474
            stage_intern = 100
      else if(ans > 450) then
            stage = 500
            start = 2973
            stage_intern = 50
      else if(ans > 400) then
            stage = 450
            start = 2522
            stage_intern = 50
      else if(ans > 350) then
            stage = 400
            start = 2121
            stage_intern = 50
      else if(ans > 300) then
            stage = 350
            start = 1770
            stage_intern = 50
      else if(ans > 250) then
            stage = 300
            start = 1469
            stage_intern = 50
      else if(ans > 200) then
            stage = 250
            start = 1218
            stage_intern = 50
      else if(ans > 150) then
            stage = 200
            start = 1017
            stage_intern = 50
      else if(ans > 100) then
            stage = 150
            start = 866
            stage_intern = 50
      else if(ans > 90) then
            stage = 100
            start = 765
            stage_intern = 10
      else if(ans > 80) then
            stage = 90
            start = 674
            stage_intern = 10
      else if(ans > 70) then
            stage = 80
            start = 593
            stage_intern = 10
      else if(ans > 60) then
            stage = 70
            start = 522
            stage_intern = 10
      else if(ans > 50) then
            stage = 60
            start = 461
            stage_intern = 10
      else if(ans > 45) then
            stage = 50
            start = 410
            stage_intern = 5
      else if(ans > 40) then
            stage = 45
            start = 364
            stage_intern = 5
      else if(ans > 35) then
            stage = 40
            start = 323
            stage_intern = 5
      else if(ans > 30) then
            stage = 35
            start = 287
            stage_intern = 5
      else if(ans > 25) then
            stage = 30
            start = 256
            stage_intern = 5
      else if(ans > 20) then
            stage = 25
            start = 230
            stage_intern = 5
      else if(ans > 19) then
            stage = 20
            start = 209
            stage_intern = 2
      else if(ans > 18) then
            stage = 19
            start = 189
            stage_intern = 2
      else if(ans > 17) then
            stage = 18
            start = 170
            stage_intern = 2
      else if(ans > 16) then
            stage = 17
            start = 152
            stage_intern = 2
      else if(ans > 15) then
            stage = 16
            start = 135
            stage_intern = 2
      else if(ans > 14) then
            stage = 15
            start = 119
            stage_intern = 2
      else if(ans > 13) then
            stage = 14
            start = 104
            stage_intern = 2
      else if(ans > 12) then
            stage = 13
            start = 90
            stage_intern = 2
      else if(ans > 11) then
            stage = 12
            start = 77
            stage_intern = 2
      else if(ans > 10) then
            stage = 11
            start = 65
            stage_intern = 2
      else if(ans > 9) then
            stage = 10
            start = 54
            stage_intern = 2
      else if(ans > 8) then
            stage = 9
            start = 44
            stage_intern = 2
      else if(ans > 7) then
            stage = 8
            start = 35
            stage_intern = 2
      else if(ans > 6) then
            stage = 7
            start = 27
            stage_intern = 2
      else if(ans > 5) then
            stage = 6
            start = 20
            stage_intern = 2
      else if(ans > 4) then
            stage = 5
            start = 14
            stage_intern = 2
      else if(ans > 3) then
            stage = 4
            start = 9
            stage_intern = 2
      else if(ans > 2) then
            stage = 3
            start = 5
            stage_intern = 2
      else if(ans > 1) then
            stage = 2
            start = 2
            stage_intern = 2
      else
            stage = 1
            start = 0
            stage_intern = 2
      end if

      !!!!!al1=1.0d0/(stage*stage/2.0d0)         !!! old
      al1=1.0d0/(stage*stage*49.0d0/100.0d0)     !!! new

      do while((t+dt .le. tend) .and.(m .le. 1).and.(total .le. 1.0d7))
c        Runge Kutta in 'stage' stagesd0, dividided in 10 of 'stage_inter' stages
          do i1=1,5                              !!! old 4
              pas=dt/i1
              x0n=t
              do l=1, neqn
                y0n(l) = g(l)
              end do
              if(t>0) then
                  call f(neqn, t, g,  g_calc)
                  feval = feval + 1
                  do l=1, neqn
                      g_save(l) = g_calc(l)
                  end do
              end if
              do j=1,i1
                  r=x0n
                  do l=1, neqn
                      g_work(l,1) = y0n(l)
                  end do
c                  g_work(:,1)=y0n;
                  do i=1,stage/stage_intern
                      if(j==1 .and. i==1) then
                          do l=1, neqn
                             g_calc(l) = g_save(l)
                          end do
                      else
                          do l=1, neqn
                             g_help(l) = g_work(l,1+(i-1)*stage_intern)
                          end do
c                          g_help = g_work(:,1+(i-1)*stage_intern);
                          call f(neqn, r, g_help,  g_calc)
                          feval = feval + 1
                      end if
                      do l=1, neqn
                         g_work(l,2+(i-1)*stage_intern)=
     *                            g_work(l,1+(i-1)*stage_intern)
     *                            +pas*al1*g_calc(l)
                      end do
                      r=x0n+(1+stage_intern*stage_intern*(i-1))*pas*al1
                      do k=2,stage_intern
                          do l=1, neqn
                               g_help(l) =
     *                               g_work(l,1+k+(i-1)*stage_intern-1);
                          end do
                          call f(neqn, r, g_help,  g_calc)
                          feval = feval + 1
                          do l=1, neqn
                              g_work(l,1+k+(i-1)*stage_intern)=
     *                              2.0d0*g_work(l,1+k+(i-1)*
     *                              stage_intern-1)-g_work(l,
     *                              1+k+(i-1)*stage_intern-2)
     *                              +2.0d0*pas*al1*g_calc(l)
                          end do
                          r=x0n+al1*(k*k+
     *                          stage_intern*stage_intern*(i-1))*pas
                      end do
                  end do
                  !!! new, if stage is odd
                  if (mod(stage,stage_intern)==1) then
                     do l=1, neqn
                        g_help(l) =
     *                       g_work(l,stage);
                     end do
                     call f(neqn, r, g_help,  g_calc)
                     feval = feval + 1
                     do l=1, neqn
                        g_work(l,stage+1)=
     *                    g_work(l,stage)+pas*al1*g_calc(l)
                     end do
                  end if
                  !!!
                  do k=1,stage+1
                     do l=1,neqn
                        sum(l) = sum(l)+b(start+k)*g_work(l,k)
                     end do
                  end do
                  do l=1,neqn
                     y(l,i1) = sum(l)
                     y0n(l) = sum(l)
                     sum(l) = 0.0d0
                  end do
                  x0n=x0n+pas
              end do
          end do
          err = 0.0d0
          do l=1, neqn
             g_oth(l) = g(l)
             !g(l) = 64.0d0/6.0d0*y(l,4)-81.0d0/6.0d0*y(l,3)+
             !*              24.0d0/6.0d0*y(l,2)-1.0d0/6.0d0*y(l,1)   !!! old
             g(l)=625.0d0/24.0d0*y(l,5)-1024.0d0/24.0d0*y(l,4)+       !!! new
     *      486.0d0/24.0d0*y(l,3)-64.0d0/24.0d0*y(l,2)+
     *      1.0d0/24.0d0*y(l,1)
             sc(l) = tol + tol*max(g(l),g_oth(l))
             !err=err+((-y(l,1)/6.0d0+2.0d0*y(l,2)-9.0d0/2.0d0*y(l,3) !!! old
             !*                +8.0d0/3.0d0*y(l,4))*2.0d0/sc(l))**2
             err=err+((y(l,1)/24.0d0-4.0d0*y(l,2)/3.0d0+              !!! new
     *                27.0d0/4.0d0*y(l,3)-32.0d0/3.0d0*y(l,4)+
     *                125.0d0/24.0d0*y(l,5))*2.0d0/sc(l))**2
          end do

c          do l=1, neqn
c             sc(l) = tol + tol*max(g(l),g_oth(l))
c          end do

c          do i=1, neqn
c             err=err+((-y(i,1)/6.0d0+2.0d0*y(i,2)-9.0d0/2.0d0*y(i,3)
c     *                 +8.0d0/3.0d0*y(i,4))*2.0d0/sc(i))**2
c          end do
          err = sqrt(err/neqn)
          
          t = t + dt
          call f(neqn, t, g, g_calc)
          feval = feval + 1
          do i=1, neqn
             var_2(i) = g_calc(i)
          end do

c     reject if error too big
          if(isnan(err) .or. 1.0d0/err .le. 1.0d0) then
              rejected = rejected + 1
              rej_counter = 0
c     reset time
              tfail = t
              t = t - dt
c     calculate new timestep
c              fac = 0.9d0*((tol/err)**0.35d0)*
c     &              ((err_old/tol)**0.2d0)
c              fac = min(10.0d0,fac)
c              fac = min(0.9d0*sqrt(tol/err),fac)
c              fac = max(fac,0.1d-3)
              fac = 0.8d0*(1.0d0/err)**(1.0d0/5.0d0)
              fac = min(10.0d0,fac)
              fac = max(fac,0.1d-3)
              dt = dt*fac
              dt = min(dtmax, dt)
c     reset g
              do i=1, neqn
                 g(i) = g_rej(i)
                 g_help(i) = g(i)
                 g_calc(i) = var_1(i)
              end do
c     accept step
         else
             accepted = accepted + 1
             if(total>1) then
c     calculate new timestep
c                fac = 0.9d0*((tol/err)**0.35d0)*
c     &              ((err_old/tol)**0.2d0)
                fac = 0.8d0*(1.0d0/err)**(1.0d0/5.0d0)
                if (t < tfail) then
                   sec = 1
                else
                   sec = 0
                end if
                
                if((rej_counter < 2) .and. sec == 1) then
                   rej_counter = rej_counter + 1
                   fac = min(1.0d0,fac)
                else if((rej_counter < 5) .and. (sec == 1)) then
                   rej_counter = rej_counter + 1
                   fac = min(2.5d0,fac)
                else
                   fac = min(10.0d0,fac)
                end if
                fac = max(fac,0.1d0)
                dt = dt*fac
                dt = min(dtmax, dt)
             end if
c     accept step just calculated
             do i=1, neqn
                g_rej(i) = g(i)
                g_help(i) = g(i)
                g_calc(i) = var_2(i)
                var_1(i) = g_calc(i)
             end do
             if(stage > max_stages) then
                 max_stages = stage
             end if
         end if
         total = total + 1
c     shorten last step if too big
         if(t+dt>tend) then
            dt=tend-t
            m=m+1
         end if
c     choose number of stages to use for next step
c     if h_new too bigd0, make it smaller
c     spcrad(1) = 0  -> subroutine rho is providing an upper bound on the spectralradius
c               = 1  -> compute an upper bound on the spectralradius using a nonlinear power method
c     spcrad(2) = 0  ->  The Jacobian may not be constant.
c               = 1  ->  The Jacobian is constant.
         if(spcrad(2)==0) then
            if(spcrad(1)==0) then
                eigmax=rho(neqn,t,g)
c                write(*,*) 'own rho =', eigmax
            else
                call SERKrho(neqn,t,f,g_help,g_calc,iwork,hmax,work,
     *                       eigmax,idid)
c                write(*,*) 'internal rho =', eigmax
c                write(*,*) 'number of func evals', iwork(9)
            end if
         end if
         dt = min(0.98d0*2000.0d0*2000.0d0/eigmax, dt)
         !write(*,*) dt
         ans = sqrt( eigmax * dt/0.98d0)

         if(ans > 1800) then
            stage = 2000
            start = 13483
            stage_intern = 200
         else if(ans > 1600) then
            stage = 1800
            start = 11682
            stage_intern = 200
         else if(ans > 1400) then
            stage = 1600
            start = 10081
            stage_intern = 200
         else if(ans > 1200) then
            stage = 1400
            start = 8680
            stage_intern = 200
         else if(ans > 1000) then
            stage = 1200
            start = 7479
            stage_intern = 200
         else if(ans > 900) then
            stage = 1000
            start = 6478
            stage_intern = 100
         else if(ans > 800) then
            stage = 900
            start = 5577
            stage_intern = 100
         else if(ans > 700) then
            stage = 800
            start = 4776
            stage_intern = 100
         else if(ans > 600) then
            stage = 700
            start = 4075
            stage_intern = 100
         else if(ans > 500) then
            stage = 600
            start = 3474
            stage_intern = 100
         else if(ans > 450) then
            stage = 500
            start = 2973
            stage_intern = 50
         else if(ans > 400) then
            stage = 450
            start = 2522
            stage_intern = 50
         else if(ans > 350) then
            stage = 400
            start = 2121
            stage_intern = 50
         else if(ans > 300) then
            stage = 350
            start = 1770
            stage_intern = 50
         else if(ans > 250) then
            stage = 300
            start = 1469
            stage_intern = 50
         else if(ans > 200) then
            stage = 250
            start = 1218
            stage_intern = 50
         else if(ans > 150) then
            stage = 200
            start = 1017
            stage_intern = 50
         else if(ans > 100) then
            stage = 150
            start = 866
            stage_intern = 50
         else if(ans > 90) then
            stage = 100
            start = 765
            stage_intern = 10
         else if(ans > 80) then
            stage = 90
            start = 674
            stage_intern = 10
         else if(ans > 70) then
            stage = 80
            start = 593
            stage_intern = 10
         else if(ans > 60) then
            stage = 70
            start = 522
            stage_intern = 10
         else if(ans > 50) then
            stage = 60
            start = 461
            stage_intern = 10
         else if(ans > 45) then
            stage = 50
            start = 410
            stage_intern = 5
         else if(ans > 40) then
            stage = 45
            start = 364
            stage_intern = 5
         else if(ans > 35) then
            stage = 40
            start = 323
            stage_intern = 5
         else if(ans > 30) then
            stage = 35
            start = 287
            stage_intern = 5
         else if(ans > 25) then
            stage = 30
            start = 256
            stage_intern = 5
         else if(ans > 20) then
            stage = 25
            start = 230
            stage_intern = 5
         else if(ans > 19) then
            stage = 20
            start = 209
            stage_intern = 2
         else if(ans > 18) then
            stage = 19
            start = 189
            stage_intern = 2
         else if(ans > 17) then
            stage = 18
            start = 170
            stage_intern = 2
         else if(ans > 16) then
            stage = 17
            start = 152
            stage_intern = 2
         else if(ans > 15) then
            stage = 16
            start = 135
            stage_intern = 2
         else if(ans > 14) then
            stage = 15
            start = 119
            stage_intern = 2
         else if(ans > 13) then
            stage = 14
            start = 104
            stage_intern = 2
         else if(ans > 12) then
            stage = 13
            start = 90
            stage_intern = 2
         else if(ans > 11) then
            stage = 12
            start = 77
            stage_intern = 2
         else if(ans > 10) then
            stage = 11
            start = 65
            stage_intern = 2
         else if(ans > 9) then
            stage = 10
            start = 54
            stage_intern = 2
         else if(ans > 8) then
            stage = 9
            start = 44
            stage_intern = 2
         else if(ans > 7) then
            stage = 8
            start = 35
            stage_intern = 2
         else if(ans > 6) then
            stage = 7
            start = 27
            stage_intern = 2
         else if(ans > 5) then
            stage = 6
            start = 20
            stage_intern = 2
         else if(ans > 4) then
            stage = 5
            start = 14
            stage_intern = 2
         else if(ans > 3) then
            stage = 4
            start = 9
            stage_intern = 2
         else if(ans > 2) then
            stage = 3
            start = 5
            stage_intern = 2
         else if(ans > 1) then
            stage = 2
            start = 2
            stage_intern = 2
         else
            stage = 1
            start = 0
            stage_intern = 2
         end if
        

         !al1=1.0d0/(stage*stage/2.0d0)    !!! old
         al1=1.0d0/(stage*stage*49.0d0/100.0d0)   !!! new
      end do

c     everthing fine
      idid = 1
c     number of function evaluations
      iwork(5) = feval
c     number of accepted steps
      iwork(7) = accepted
c     number of rejected steps
      iwork(8) = rejected
c     max number of stages used
      iwork(10) = max_stages

      return
      end

