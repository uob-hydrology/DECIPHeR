!% Toby Dunne
!% Mar 2016
! really behaves like a stack not a queue
! last in first out
module dta_stack
    use dta_utility
    implicit none
    integer:: block_size
    parameter(block_size=8192)

    private
    public :: stack_type, stack_test, new_stack
    integer, save :: next_block_id

    type stack_block
        integer :: block_id
        type(point_type), dimension(block_size) :: data
        type(stack_block), pointer :: next => null()
    end type

    type stack_type
        private

        type(stack_block), pointer :: data_head => null()

        type(stack_block), pointer :: tmp => null()

        integer :: head               ! head pointer

        integer, public :: item_count      ! current number of elements
    contains

        procedure, public :: to_array
        procedure, public :: info
        procedure, public :: push
        procedure, public :: pop
        procedure, public :: clear
        procedure, public :: free
    end type stack_type

contains

    subroutine stack_test()
        implicit none

        type(stack_type) :: qq
        type(point_type) :: qitem
        integer :: i

        qq = new_stack()
        print *, 'INIT STACK'
        call qq%info()

        print *, 'TEST ENQUEUE STACK'

        do i=1,10
            print *, 'enqueue : ', i
            qitem%x = i
            call qq%push(qitem)
            !call qq%info(qq)
            !print *, 'expect count : ', i ,' got : ', qq%item_count
        enddo

        call qq%info()

        print *, 'TEST LIFO STACK'

        do
            if(qq%item_count == 0) then
                print *, 'empty'
                exit
            endif
            !print *, 'dequeue : ', i
            qitem = qq%pop()
            !call qq%info(qq)
            !print *, 'expect count 2: ', qq%item_count
            print *, 'value: ', qitem%x
        enddo

        call qq%info()

        do i=1,10
            print *, 'push : ', i
            qitem%x = i
            call qq%push(qitem)
            !call qq%info(qq)
            !print *, 'expect count : ', i ,' got : ', qq%item_count
        enddo

        call qq%info()
        print *, 'TEST LIFO STACK'

        do
            if(qq%item_count == 0) then
                print *, 'empty '
                exit
            endif
            !print *, 'dequeue : ', i
            qitem = qq%pop()
            !call qq%info(qq)
            !print *, 'expect count 2: ', qq%item_count
            print *, 'value: ', qitem%x
        enddo

        call qq%info()

        print *, 'TEST ENQUEUE STACK'

        do i=1,10
            print *, 'enqueue : ', i
            qitem%x = i
            call qq%push(qitem)
            !call qq%info(qq)
            !print *, 'expect count : ', i ,' got : ', qq%item_count
        enddo

        call qq%info()


        print *, 'TEST LIFO STACK - one one'
        do i=1,20
            qitem%x = i
            call qq%push(qitem)
            qitem%x = -99

            qitem = qq%pop()
            print *, 'value: ', qitem%x

        enddo


    end subroutine stack_test

    function new_stack()
        implicit none
        type(stack_type) :: new_stack

        type(stack_block), pointer :: qblock

        allocate(qblock)
        qblock%block_id = next_block_id
        next_block_id = next_block_id + 1

        new_stack%data_head => qblock

        new_stack%head = 0 ! not available
        new_stack%item_count = 0

    end function

    subroutine to_array(q, array)
        implicit none
        class(stack_type) :: q
        type(point_type) :: array(q%item_count)
        integer :: i
        i = 1
        do
            if(q%item_count == 0) exit
            array(i) = q%pop()
            i = i+1
            !print *, 'i:',i, ' x:', array(i)%x, ' y:', array(i)%y
        enddo
    end subroutine to_array

    subroutine info(q)
        implicit none
        class(stack_type) :: q

        print *, 'count:', q%item_count &
            , 'head :', q%head &
            , 'data_head:', q%data_head%block_id


    end subroutine

    ! clear - ready for more data
    subroutine clear(q)
        implicit none
        class(stack_type) :: q

        type(stack_block), pointer :: next

        do while (associated(q%data_head%next))
            next => q%data_head%next
            deallocate(q%data_head)
            NULLIFY(q%data_head)
            q%data_head => next
        enddo

        q%head = 0
        q%item_count = 0
    end subroutine

!
    subroutine free(q)
        implicit none
        class(stack_type) :: q

        type(stack_block), pointer :: next

        do while (associated(q%data_head))
            next => q%data_head%next
            deallocate(q%data_head)
            !NULLIFY(q%data_head)
            q%data_head => next
        enddo

        if(associated(q%tmp)) then
            deallocate(q%tmp)
            NULLIFY(q%tmp)
        endif

        q%head = 0
        q%item_count = 0
    end subroutine


    !------------------------- Enqueue ----------------------------
    ! Syntax:       q%Enqueue(item)
    !
    ! Inputs:       item is a numeric value
    !
    ! Description:  Adds item to tail end of Q
    !--------------------------------------------------------------
    subroutine push(q, item)
        implicit none
        class(stack_type) :: q
        type(point_type) :: item

        if(q%head == block_size) then
            !print *, 'Add block:', next_block_id
            ! add another block on the head end
            if(associated(q%tmp)) then
                ! if the tmp block is available, reuse this
            else
                allocate(q%tmp)
            endif
            q%tmp%next => q%data_head
            q%data_head => q%tmp
            nullify( q%tmp )

            q%data_head%block_id = next_block_id
            next_block_id = next_block_id + 1

            ! set the current head to the newly allocated block

            q%head = 1
        else
            q%head = q%head + 1
        endif
        q%item_count = q%item_count + 1

        !print *, 'set index: ', q%data_head%block_id, ':' , q%head,   'value: ', item%value
        q%data_head%data(q%head) = item

    end subroutine

    function pop(q) result(item)
        implicit none
        class(stack_type) :: q
        type(point_type) :: item

        if(q%item_count > 0) then
            if(q%head < 1) then
                if(associated(q%data_head%next)) then
                    if (associated(q%tmp)) then
                        deallocate(q%tmp)
                    endif
                    q%tmp => q%data_head
                    q%data_head => q%data_head%next

                    ! set head to the last item on previous block
                    q%head = block_size
                endif
            endif

            item = q%data_head%data(q%head)

            q%head = q%head - 1
            q%item_count = q%item_count - 1


        else
            print *, 'Stack empty'
            stop
        endif

    end function


end module dta_stack
