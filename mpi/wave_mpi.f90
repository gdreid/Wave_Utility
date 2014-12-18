      program wave_single
      implicit none
      include 'mpif.h'

      integer, parameter       :: sizex = 2048
      integer, parameter       :: sizey = 2048
      integer, parameter       :: b_nx = 8
      integer, parameter       :: b_ny = 8
      integer, parameter       :: itterations = 10

      integer                  :: i, j, k, itt, max_itt
      integer                  :: nx, ny
      double precision         :: dx, dy, dt
      double precision         :: xmin, xmax, ymin, ymax
      double precision         :: r, t_start, t_end
      double precision         :: res, jac
      real                     :: start_time, stop_time

      double precision, dimension (:,:), allocatable :: &
      y0_n, y1_n, y0_np1, y1_np1, temp0, temp1
      double precision, dimension (:), allocatable   :: x, y

      integer                  :: num_procs, ierr, my_rank, tag
      integer                  :: status(MPI_STATUS_SIZE)
      integer                  :: starty, startx, f

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)

      f  = sizey/num_procs
      ny = f
      if (my_rank .lt. (sizey-f*num_procs)) then
         ny = ny + 1
      endif
      nx = sizex

      startx = 0
      starty = 1
      do i = 0, my_rank-1, 1
         starty = starty + f
         if (i .lt. (sizey-f*num_procs)) then
            starty = starty + 1
         endif
      enddo 
      starty = starty -1

      allocate(y0_n(0:nx+1,0:ny+1))
      allocate(y1_n(0:nx+1,0:ny+1))
      allocate(y0_np1(0:nx+1,0:ny+1))
      allocate(y1_np1(0:nx+1,0:ny+1))
      allocate(temp0(0:nx+1,0:ny+1))
      allocate(temp1(0:nx+1,0:ny+1))
      allocate(x(0:nx+1))
      allocate(y(0:ny+1))

      xmin = -3.0
      xmax = 3.0
      ymin = -3.0
      ymax = 3.0
      max_itt = itterations * sizex/b_nx

      dx = (xmax-xmin)/DBLE(sizex+1)
      dy = (ymax-ymin)/DBLE(sizey+1)
      dt = 0.1*(xmax-xmin)/DBLE(sizex)
      
      !write(*,*) "my rank is ", my_rank, " of ", num_procs

      xmin = xmin
      do i = 0, nx+1
       x(i) = xmin + (i)*dx
      enddo

      ymin = ymin + starty*dy
      do j = 0, ny+1
       y(j) = ymin + (j)*dy
      enddo

      do i = 0, nx+1
       do j = 0, ny+1
        y0_n(i,j)   = 0
        y1_n(i,j)   = 0
        y0_np1(i,j) = 0
        y1_np1(i,j) = 0 
       enddo
      enddo

      !write(*,*) "rank: ", my_rank, " 1"

      do i = 1, nx
       do j = 0, ny+1
        if ((j + starty .ne. 0) .and. (j + starty .ne. sizey + 1)) then
         y0_n(i,j)   = dexp(-x(i)*x(i)-y(j)*y(j))
         y0_np1(i,j) = dexp(-x(i)*x(i)-y(j)*y(j))
        endif
       enddo
      enddo      

      !write(*,*) "rank: ", my_rank, " 2"

      call cpu_time(start_time) 

      do itt = 0, max_itt - 1
       !write(*,*) "rank: ", my_rank, " 3 itt: ", itt
       do k = 0, 10
        do i = 1, nx
         do j = 1, ny
          if ((j + starty .ne. 0) .and. &
           (j + starty .ne. sizey + 1)) then 
           res = (y0_np1(i,j) - y0_n(i,j))/dt & 
            - 0.5*(y1_np1(i,j)+y1_n(i,j))
            jac = 1d0/dt
           temp0(i,j) = y0_np1(i,j) - res/jac
           !temp0(i,j) = y0_np1(i,j)

           res = (y1_np1(i,j) - y1_n(i,j))/dt & 
            -0.5*(y0_np1(i-1,j)-2*y0_np1(i,j)+y0_np1(i+1,j))/(dx**2) & 
            -0.5*(y0_np1(i,j-1)-2*y0_np1(i,j)+y0_np1(i,j+1))/(dy**2) & 
            -0.5*(y0_n(i-1,j)-2*y0_n(i,j)+y0_n(i+1,j))/(dx**2) & 
            -0.5*(y0_n(i,j-1)-2*y0_n(i,j)+y0_n(i,j+1))/(dy**2)
           jac = 1d0/dt
           temp1(i,j) = y1_np1(i,j) - res/jac
           !temp1(i,j) = y1_np1(i,j)
          endif
         enddo !j
        enddo !i

        do i = 1, nx
         do j = 0, ny+1
          if ((j + starty .ne. 0) .and. &
           (j + starty .ne. sizey + 1)) then
           y0_np1(i,j) = temp0(i,j)
           y1_np1(i,j) = temp1(i,j)
          endif
         enddo !j
        enddo !i

        !MPI Communication
        if (my_rank .gt. 0) then
        ! call mpi_sendrecv(temp0(0,1),nx+2,MPI_DOUBLE, my_rank-1, tag, &
        !  y0_np1(0,0),nx+2,MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &
        !  status, ierr)
        endif
        !write(*,*) "1: my rank: ", my_rank
        if (my_rank .lt. num_procs -1) then
        ! call mpi_sendrecv(temp0(0,ny),nx+2,MPI_DOUBLE, my_rank+1, &
        !  tag, y0_np1(0,ny+1),nx+2,MPI_DOUBLE, my_rank+1, tag, &
        !  MPI_COMM_WORLD, status, ierr)
        endif
        !write(*,*) "2: my rank: ", my_rank
        call mpi_barrier(MPI_COMM_WORLD, ierr)

        !MPI Communication
        if (my_rank .gt. 0) then
        ! call mpi_sendrecv(temp1(0,1),nx+2,MPI_DOUBLE, my_rank-1, tag, &
        !  y1_np1(0,0),nx+2,MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &
        !  status, ierr)
        endif
        !write(*,*) "1: my rank: ", my_rank
        if (my_rank .lt. num_procs -1) then
        ! call mpi_sendrecv(temp1(0,ny),nx+2,MPI_DOUBLE, my_rank+1, &
        !  tag, y1_np1(0,ny+1),nx+2,MPI_DOUBLE, my_rank+1, tag, &
        !  MPI_COMM_WORLD, status, ierr)
        endif
        !write(*,*) "2: my rank: ", my_rank
        call mpi_barrier(MPI_COMM_WORLD, ierr)

       enddo !k
        

       do i = 1, nx
        do j = 0, ny+1
         if ((j + starty .ne. 0) .and. &
          (j + starty .ne. sizey + 1)) then
          res         = y0_n(i,j)
          y0_n(i,j)   = y0_np1(i,j)
          y0_np1(i,j) = res

          res         = y1_n(i,j)
          y1_n(i,j)   = y1_np1(i,j)
          y1_np1(i,j) = res
         endif
        enddo
       enddo
       
      enddo

      call cpu_time(stop_time)
      if (my_rank .eq. 0) then 
       write(*,*) "num procs: ", num_procs, &
        " total time: ", stop_time-start_time
      endif
      !call print2d(y0_n, 0, nx+2, 0, ny+2)
      !call print2d(y1_n, 0, nx+2, 0, ny+2)
      !call print2d(y0_np1, 0, nx+2, 0, ny+2)
      !call print2d(y1_np1, 0, nx+2, 0, ny+2)

      call mpi_finalize(ierr)
      end

      subroutine print2d(y, sx, nx, sy, ny)
      implicit none
      integer                            ::  nx, ny, i, j, sx, sy
      double precision, dimension(0:nx-1,0:ny-1) ::  y

      do j = sy, sy+ny-1
       do i = sx, sx+nx-1
        !if (y(i,j) > 0.05) then   
        ! write(*,"(i2)",advance="no") 1
        !elseif (y(i,j) < -0.05) then
        ! write(*,"(i2)",advance="no") -1
        !else
        ! write(*,"(i2)", advance="no") 0
        !endif
        write(*,"(f7.3)",advance="no") y(i,j)
       enddo
       write(*,*) ""
      enddo
      write(*,*) ""

      end
