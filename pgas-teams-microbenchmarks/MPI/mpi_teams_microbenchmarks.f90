!
! Microbenchmarks for CAF
!
! (C) HPCTools Group, University of Houston
!
module mpi_teams_microbenchmarks
  implicit none

  integer, parameter :: SYNC_NITER = 10000
  integer, parameter :: BUFFER_SIZE = 1 !1024 !*1024
  integer, parameter :: RED_NITER = 1
  integer, parameter :: NITER = 1000
  integer, parameter :: N=1024
  integer, parameter :: NUM_TIME_G = 32

  integer            :: send_buffer(BUFFER_SIZE)
  integer            :: recv_buffer(BUFFER_SIZE) 

  integer :: rank, np
contains

!*************************barrier************************

  subroutine mpi_barrier_team(nteams)
    use mpi

    implicit none

    integer :: ierr
    integer :: team_id, new_comm
    integer :: i
    double precision :: t1, t2, time
    integer :: me, ne
    integer :: nteams
    if(nteams .eq. 1) then
       team_id = 1
    else
       team_id = mod(rank, 2)
    end if
    if(rank == 0)       print *, "-----------barrier------------"

    call mpi_comm_split(mpi_comm_world, team_id, rank, new_comm, ierr)

    call mpi_comm_rank(new_comm, me, ierr)
    call mpi_comm_size(new_comm, ne, ierr)

    t1 = mpi_wtime()
    do i = 1, SYNC_NITER
        call mpi_barrier(new_comm, ierr)
     end do
     t2 = mpi_wtime()
     t2 = 1000000*(t2-t1)/(SYNC_NITER)

     call mpi_reduce(t2, time, 1, mpi_double, mpi_sum, 0, new_comm, ierr)

     time = time / ne

     if (me == 0) then
        print *, "MPI Time for MPI_Barrier team:",  time, "us ", " for number of iterations:", SYNC_NITER
     end if
     if(rank == 0)       print *, "number of processes per team for MPI_barrier is: ", np/nteams
  end subroutine mpi_barrier_team

!*************************reduction************************

  subroutine mpi_co_sum_team(nteams)
    use mpi
    implicit none
    integer :: ierr
    integer :: team_id, new_comm
    integer :: i
    double precision :: t1, t2, time
    integer :: me, ne
    integer :: nrep, blksize
    integer :: num_time
    integer :: nteams
    if(nteams .eq. 1) then
       team_id = 1
    else
       team_id = mod(rank, 2)
    end if
    if(rank == 0)       print *, "-----------barrier------------"

    call mpi_comm_split(mpi_comm_world, team_id, rank, new_comm, ierr)

    call mpi_comm_rank(new_comm, me, ierr)
    call mpi_comm_size(new_comm, ne, ierr)

    blksize = 1
    num_time = 1
    do while (blksize <= BUFFER_SIZE)
       send_buffer(1:blksize) = 42
       nrep = RED_NITER

       t1 = mpi_wtime()
       do i = 1, nrep
          call mpi_allreduce(send_buffer, recv_buffer, &
               &blksize, MPI_INTEGER, MPI_SUM, new_comm, ierr)
       end do

       t2 = mpi_wtime()
       t2 = 1000000*(t2-t1)/(RED_NITER)

       call mpi_reduce(t2, time, 1, mpi_double, mpi_sum, 0, new_comm, ierr)
       time = time / np

       call mpi_barrier(new_comm, ierr)
       if (me == 0) then
          print *, "MPI Time for  MPI_ALLREDUCE SUM within team:", time, "us ", " for data size:", blksize, " and number of iterations:",  RED_NITER
       end if

       num_time = num_time + 1
       blksize = blksize + 2
    end do
     if(rank == 0)       print *, "number of processes per team for MPI_ALLREDUCE is: ", np/nteams

  end subroutine mpi_co_sum_team


end module mpi_teams_microbenchmarks

program main
  use mpi_teams_microbenchmarks
  use mpi
  implicit none
  integer :: ierr
  integer :: team_id, new_comm
  integer :: i
  double precision :: t1, t2, time
  integer :: me, ne

  team_id = mod(rank, 2)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  if(rank == 0)       print *, "------------begin MPI------------"

  call mpi_barrier_team(1)
  call mpi_co_sum_team(1)
  call mpi_barrier_team(2)
  call mpi_co_sum_team(2)

  if(rank == 0)       print *, "------------end MPI------------"

  call mpi_finalize(ierr)
end program main

