
module caf_teams
   type team_type
      private
      integer (kind=8) :: h = 0
   end type team_type
end module caf_teams

module caf_teams_microbenchmarks
  !use, intrinsic :: iso_fortran_env
  use caf_teams      
  implicit none

  include 'mpif.h'
  
  
  
  interface
     subroutine form_team(id, new_team, new_index)
       use caf_teams      
       !use, intrinsic :: iso_fortran_env
       implicit none
       integer, intent(in) :: id
       type(team_type), intent(out) :: new_team
       integer, optional, intent(in) :: new_index
     end subroutine form_team

     subroutine change_team(team)
       !use, intrinsic :: iso_fortran_env
       use caf_teams      
       implicit none
       type(team_type), intent(in) :: team
     end subroutine change_team
     
     subroutine end_change_team()
     end subroutine end_change_team

     function team_id()
       integer :: team_id
     end function team_id
  end interface

  
  integer :: tid
  integer, parameter :: SYNC_NITER = 1000
  integer, parameter :: BUFFER_SIZE = 1!1024 *1024
  integer, parameter :: RED_NITER = 10000
  integer, parameter :: NITER = 1000
  integer, parameter :: N=1024
  integer, parameter :: NUM_TIME = 32
  
  integer :: send_buffer(BUFFER_SIZE)[*]
  integer :: recv_buffer(BUFFER_SIZE)[*]
  integer :: A(N)[*]
  integer :: sync_correctness[*]
  double precision, allocatable :: time(:)[:]
  
  type(team_type) :: tm
  double precision :: t1, t2

  integer :: me, i, j
  
contains

!*********************************sync team*****************************!  
  subroutine sync_team(nteams)
    integer :: nteams
    me = this_image()
    if(nteams .eq. 1) then
       tid = 1
    else 
       tid = 1 + mod(this_image()-1, 2)
    end if
    if(this_image() == 1)       print *, "------------barrier------------"   
       
    call form_team (tid, tm)
    call change_team (tm)

    sync_correctness = 0;
    
    t1 = MPI_WTIME()
    do i = 1, SYNC_NITER
       sync all
       sync_correctness = sync_correctness + 1
    end do
    t2 = MPI_WTIME()

    sync all

    if(this_image()==1) then
          do j = 1, num_images()
             if (sync_correctness[j] /= SYNC_NITER) then
                write (*, '("sync_correctness wrong value is:",I20)'), sync_correctness[j]
                error stop "got wrong value"
             end if
          end do
    end if
    
    time(1) = 1000000*(t2-t1)/(SYNC_NITER)
    
    if (this_image() == 1) then
       do i = 2, num_images()
          time(1) = time(1) + time(1)[i]
       end do
      print *, "UHCAF Time for sync_team:",  time(1)/num_images(),"us ", " for number of iterations:", SYNC_NITER
   end if
   call end_change_team()
   if(this_image() == 1)       print *, "number of images per team for sync team is: ", num_images()/nteams

 end subroutine sync_team

!******************************** reduction *****************************
 subroutine co_sum_team(nteams)
   implicit none
   
   double precision :: t1, t2
   integer :: nteams
   integer ::  ni, nrep
   integer :: num_time
   integer :: blksize
   
   me = this_image()
   
   if(nteams .eq. 1) then
      tid = 1
   else 
      tid = 1 + mod(this_image()-1, 2)
   end if
   if(this_image() == 1)       print *, "------------reduction------------"   
   
   call form_team (tid, tm)
    
   call change_team (tm)
   
   num_time = 1
   blksize = 1
   do while (blksize <= BUFFER_SIZE)
      send_buffer(1:blksize) = 42
      nrep = RED_NITER
      
      t1 = MPI_WTIME()
      do i = 1, nrep
         call co_sum(send_buffer(1:blksize), recv_buffer(1:blksize));
         if (any(recv_buffer(1:blksize) /= 42*num_images())) then
            error stop "got wrong value"
            write (*, '("recv buffer wrong value is:",I20)'),recv_buffer(1)
         end if
      end do
      t2 = MPI_WTIME()

      time(num_time) = 1000000*(t2-t1)/(RED_NITER)
      
      sync all
      
      if (this_image() == 1) then
         do i = 2, num_images()
            time(num_time) = time(num_time) + &
                 time(num_time)[i]
         end do
         print *, "UHCAF Time for reduction sum within team:", time(num_time)/num_images(), "us ", " for data size:", blksize, " and number of iterations:", nrep
      end if

      num_time = num_time + 1
      blksize = blksize * 2
   end do
   call end_change_team()
   if(this_image() == 1)       print *, "number of images per team for reduce team is: ", num_images()/nteams


 end subroutine co_sum_team

!********************************broadcast *****************************
 subroutine broadcast_team(nteams)
   integer :: nteams
   me = this_image()
   
   if(nteams .eq. 1) then
      tid = 1
   else 
      tid = 1 + mod(this_image()-1, 2)
   end if

   if(this_image() == 1)       print *, "------------broadcast------------"   
   call form_team (tid, tm)
   call change_team (tm)
   
   sync all
   t1 = MPI_WTIME()
   do i=1,NITER
      A(:) = this_image()
      call co_broadcast(A,1)
      if (any(A(1:N) /= 1)) then
         error stop "got wrong value"
      end if
   end do
   t2 = MPI_WTIME()
   time(1) = 1000000*(t2-t1)/(SYNC_NITER)
   
   if (me == 1) then
      do i = 2, num_images()
         time(1) = time(1) + time(1)[i]
      end do
      print *, "UHCAF Time for broadcast team:",  time(1)/num_images(),"us", " for number of iterations:", NITER
   end if

   call end_change_team()
   if(this_image() == 1)       print *, "number of images per team for broadcast is: ", num_images()/nteams
 end subroutine broadcast_team

 !***********************form team ***********************************
 subroutine form_one_team_overhead
   me = this_image()
   
   tid = 1 
   if(this_image() == 1)       print *, "------------form team------------"   
   
   t1 = MPI_WTIME()
   do i = 1,NITER
      call form_team (tid, tm)
   end do
   t2 = MPI_WTIME()

   sync all
  
   time(1) = 1000000*(t2-t1)/(NITER)
   
   if (this_image() == 1) then
      do i = 2, num_images()
         time(1) = time(1) + time(1)[i]
      end do
      print *, "UHCAF Time for form team overhead:",  time(1)/num_images(),"us", " for number of iterations:", NITER
   end if
   call change_team (tm)
   sync all
   call end_change_team()
 end subroutine form_one_team_overhead
end module caf_teams_microbenchmarks

!***************************main program*********************************

program main
  use caf_teams_microbenchmarks

  implicit none
  allocate (time(NUM_TIME)[*] )
  if(this_image() == 1)       print *, "------------begin UHCAF------------"   
  call sync_team(1)
  call co_sum_team(1)
  call broadcast_team(1)
  call sync_team(2)
  call co_sum_team(2)
  call broadcast_team(2)
  call form_one_team_overhead()
  if(this_image() == 1)       print *, "------------end UHCAF------------"   
   
end program main

