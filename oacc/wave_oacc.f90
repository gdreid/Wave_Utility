      program wave_oacc
      implicit none
      integer, parameter       :: nx = 2050
      integer, parameter       :: ny = 2050
      integer, parameter       :: b_nx = 8
      integer, parameter       :: b_ny = 8
      integer, parameter       :: itterations = 10

      integer                  :: i, j, k, itt, max_itt
      double precision         :: dx, dy, dt
      double precision         :: xmin, xmax, ymin, ymax
      double precision         :: r, t_start, t_end
      double precision         :: res, jac
      real                     :: start_time, stop_time

      double precision, dimension(nx, ny) :: y0_n, y1_n, y0_np1
      double precision, dimension(nx,ny)  :: y1_np1, temp0, temp1
      double precision, dimension(nx)     :: x
      double precision, dimension(ny)     :: y


      xmin = -10.0
      xmax = 10.0
      ymin = -10.0
      ymax = 10.0
      max_itt = itterations * nx/b_nx

      dx = (xmax-xmin)/DBLE(nx-1)
      dy = (ymax-ymin)/DBLE(ny-1)
      dt = 0.1*(xmax-xmin)/DBLE(nx)
      
      do i = 1, nx
       x(i) = xmin + (i-1)*dx
      enddo

      do j = 1, ny
       y(j) = ymin + (j-1)*dy
      enddo

      do i = 1, nx
       do j = 1, ny
        y0_n(i,j)   = 0
        y1_n(i,j)   = 0
        y0_np1(i,j) = 0
        y1_np1(i,j) = 0 
       enddo
      enddo

      do i = 2, nx-1
       do j = 2, ny-1
        y0_n(i,j)   = dexp(-x(i)*x(i)-y(j)*y(j))
        y0_np1(i,j) = dexp(-x(i)*x(i)-y(j)*y(j))
       enddo
      enddo

      call cpu_time(start_time)

      !$acc data copy(y0_n,y1_n,y0_np1,y1_np1,temp0,temp1,x,y)      

      do itt = 0, max_itt - 1
       do k = 0, 10
        !$acc kernels
        do i = 2, nx-1
         do j = 2, ny-1
          res = (y0_np1(i,j) - y0_n(i,j))/dt & 
       - 0.5*(y1_np1(i,j)+y1_n(i,j))
          jac = 1d0/dt
          temp0(i,j) = y0_np1(i,j) - res/jac

          res = (y1_np1(i,j) - y1_n(i,j))/dt & 
       -0.5*(y0_np1(i-1,j)-2*y0_np1(i,j)+y0_np1(i+1,j))/(dx**2) & 
       -0.5*(y0_np1(i,j-1)-2*y0_np1(i,j)+y0_np1(i,j+1))/(dy**2) & 
       -0.5*(y0_n(i-1,j)-2*y0_n(i,j)+y0_n(i+1,j))/(dx**2) & 
       -0.5*(y0_n(i,j-1)-2*y0_n(i,j)+y0_n(i,j+1))/(dy**2)
          jac = 1d0/dt
          temp1(i,j) = y1_np1(i,j) - res/jac
         enddo
        enddo
        !$acc end kernels
        
        !$acc kernels
        do i = 2, nx-1
         do j = 2, ny-1
          y0_np1(i,j) = temp0(i,j)
          y1_np1(i,j) = temp1(i,j)
         enddo
        enddo
        !$acc end kernels
       enddo

       !$acc kernels
       do i = 1, nx
        do j = 1, ny
         res         = y0_n(i,j)
         y0_n(i,j)   = y0_np1(i,j)
         y0_np1(i,j) = res

         res         = y1_n(i,j)
         y1_n(i,j)   = y1_np1(i,j)
         y1_np1(i,j) = res
        enddo
       enddo
       !$acc end kernels
      enddo
      !$acc end data

      call cpu_time(stop_time)
      write(*,*) "total time: ", stop_time-start_time

      !call print2d(y0_n,nx,ny)
      !call print2d(y1_n,nx,ny)
      !call print2d(y0_np1,nx,ny)
      !call print2d(y1_np1,nx,ny)
      
      end

      subroutine print2d(y, nx, ny)
      implicit none
      integer                            ::  nx, ny, i, j
      double precision, dimension(nx,ny) ::  y

      do j = 1, ny
       do i = 1, nx
        if (y(i,j) > 0.05) then   
         write(*,"(i2)",advance="no") 1
        elseif (y(i,j) < -0.05) then
         write(*,"(i2)",advance="no") -1
        else
         write(*,"(i2)", advance="no") 0
        endif
       enddo
       write(*,*) ""
      enddo
      write(*,*) ""

      end
