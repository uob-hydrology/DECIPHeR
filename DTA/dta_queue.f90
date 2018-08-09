!% Toby Dunne
!% Mar 2016
module dta_queue
    use dta_utility
    implicit none
    integer:: block_size
    parameter(block_size=8192)

    private
    public :: queue_type, queue_test, new_queue
    integer, save :: next_block_id

    type queue_block
        integer :: block_id
        type(point_type), dimension(block_size) :: data
        type(queue_block), pointer :: next => null()
        type(queue_block), pointer :: prev => null()
    end type

    type queue_type
        private

        type(queue_block), pointer :: data_head => null()
        type(queue_block), pointer :: data_tail => null()

        type(queue_block), pointer :: tmp => null()

        integer :: tail               ! tail pointer
        integer :: head               ! head pointer

        integer, public :: item_count      ! current number of elements
    contains

        procedure, public :: to_array
        procedure, public :: info
        procedure, public :: enqueue
        procedure, public :: dequeue
        procedure, public :: pop
        procedure, public :: clear
        procedure, public :: free
    end type queue_type

contains

    subroutine queue_test()
        implicit none

        type(queue_type) :: qq
        type(point_type) :: qitem
        integer :: i

        qq = new_queue();
        print *, 'INIT QUEUE'
        call qq%info()

        print *, 'TEST ENQUEUE QUEUE'

        do i=1,10
            print *, 'enqueue : ', i
            qitem%x = i
            call qq%enqueue(qitem)
            !call qq%info(qq)
            !print *, 'expect count : ', i ,' got : ', qq%item_count
        enddo

        call qq%info()

        print *, 'TEST FIFO QUEUE'

        do
            if(qq%item_count == 0) then
                print *, 'empty'
                exit
            endif
            !print *, 'dequeue : ', i
            qitem = qq%dequeue()
            !call qq%info(qq)
            !print *, 'expect count 2: ', qq%item_count
            print *, 'value: ', qitem%x
        enddo

        call qq%info()

        do i=1,10
            print *, 'enqueue : ', i
            qitem%x = i
            call qq%enqueue(qitem)
            !call qq%info(qq)
            !print *, 'expect count : ', i ,' got : ', qq%item_count
        enddo

        call qq%info()
        print *, 'TEST FIFO QUEUE'

        do
            if(qq%item_count == 0) then
                print *, 'empty '
                exit
            endif
            !print *, 'dequeue : ', i
            qitem = qq%dequeue()
            !call qq%info(qq)
            !print *, 'expect count 2: ', qq%item_count
            print *, 'value: ', qitem%x
        enddo

        call qq%info()

        print *, 'TEST ENQUEUE QUEUE'

        do i=1,10
            print *, 'enqueue : ', i
            qitem%x = i
            call qq%enqueue(qitem)
            !call qq%info(qq)
            !print *, 'expect count : ', i ,' got : ', qq%item_count
        enddo

        call qq%info()


        print *, 'TEST LIFO QUEUE (STACK)'
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

        print *, 'TEST FIFO QUEUE - one one'
        do i=1,20
            qitem%x = i
            call qq%enqueue(qitem)
            qitem%x = -99

            qitem = qq%dequeue()
            print *, 'value: ', qitem%x

        enddo


    end subroutine queue_test

    function new_queue()
        implicit none
        type(queue_type) :: new_queue

        type(queue_block), pointer :: qblock

        allocate(qblock)
        qblock%block_id = next_block_id
        next_block_id = next_block_id + 1

        new_queue%data_head => qblock
        new_queue%data_tail => qblock

        new_queue%head = 0 ! not available
        new_queue%tail = 1
        new_queue%item_count = 0

    end function

    subroutine to_array(q, array)
        implicit none
        class(queue_type) :: q
        type(point_type) :: array(q%item_count)
        integer :: i
        i = 1
        do
            if(q%item_count == 0) exit
            array(i) = q%dequeue()
            i = i+1
        enddo
    end subroutine to_array

    subroutine info(q)
        implicit none
        class(queue_type) :: q

        print *, 'count:', q%item_count &
            , 'head :', q%head &
            , 'tail :', q%tail &
            , 'data_head:', q%data_head%block_id &
            , 'data_tail:', q%data_tail%block_id

        if(associated(q%data_tail%next)) then
            print *, 'tail.next:', q%data_tail%next%block_id
        else
            print *, 'tail.next: null'
        endif
        if(associated(q%data_tail%prev)) then
            print *, 'UNEXPECTED tail.prev:', q%data_tail%prev%block_id
        !else
            !print *, 'tail.prev: null'
        endif

        if(associated(q%data_head%next)) then
            print *, 'UNEXPECTED head.next:', q%data_head%next%block_id
        !else
            !print *, 'head.next: null'
        endif
        if(associated(q%data_head%prev)) then
            print *, 'head.prev:', q%data_head%prev%block_id
        else
            print *, 'head.prev: null'
        endif
    end subroutine

! clear queue - can be reused
    subroutine clear(q)
        implicit none
        class(queue_type) :: q

        type(queue_block), pointer :: prev

        do while (associated(q%data_head%prev))
            prev => q%data_head%prev
            deallocate(q%data_head)
            NULLIFY(q%data_head)
            q%data_head => prev
        enddo

        q%data_head => q%data_tail

        q%head = 0
        q%tail = 1
        q%item_count = 0
    end subroutine

    ! free data from queue
    subroutine free(q)
        implicit none
        class(queue_type) :: q

        type(queue_block), pointer :: prev

        do while (associated(q%data_head))
            prev => q%data_head%prev
            deallocate(q%data_head)
            NULLIFY(q%data_head)
            q%data_head => prev
        enddo

        NULLIFY(q%data_tail)

        if(associated(q%tmp)) then
            deallocate(q%tmp)
            NULLIFY(q%tmp)
        endif

        q%head = 0
        q%tail = 1
        q%item_count = 0
    end subroutine

    !------------------------- Enqueue ----------------------------
    ! Syntax:       q%Enqueue(item)
    !
    ! Inputs:       item is a numeric value
    !
    ! Description:  Adds item to tail end of Q
    !--------------------------------------------------------------
    subroutine enqueue(q, item)
        implicit none
        class(queue_type) :: q
        type(point_type) :: item


        if(q%head == block_size) then
            !print *, 'Add block:', next_block_id

            ! add another block on the head end
            if(associated(q%tmp)) then
                ! if the tmp block is available, reuse this
                q%data_head%next => q%tmp
                nullify( q%tmp )
            else
                allocate(q%data_head%next)
            endif

            q%data_head%next%block_id = next_block_id
            next_block_id = next_block_id + 1

            ! update the prev pointer to the current head
            q%data_head%next%prev => q%data_head
            ! set the current head to the newly allocated block
            q%data_head => q%data_head%next
            q%head = 1
        else
            q%head = q%head + 1
        endif
        q%item_count = q%item_count + 1

        !print *, 'set index: ', q%data_head%block_id, ':' , q%head,   'value: ', item%value
        q%data_head%data(q%head) = item

    end subroutine

    !------------------------- Dequeue ----------------------------
    ! Syntax:       xi = q%Dequeue()
    !
    ! Outputs:      xi is a numeric value
    !
    ! Description:  Removes tail item from Q FIFO
    !--------------------------------------------------------------
    function dequeue(q) result(item)
        implicit none
        class(queue_type) :: q
        type(point_type) :: item

        if(q%item_count > 0) then
            if(q%tail > block_size) then
                if(associated(q%data_tail%next)) then
                    q%data_tail => q%data_tail%next

                    ! set tail to the first item on next block
                    q%tail = 1

                    ! in order to prevent too many allocations
                    ! save one block value to tmp for reuse
                    if (associated(q%tmp)) then
                        deallocate(q%data_tail%prev)
                    else
                        q%tmp => q%data_tail%prev
                    endif

                    NULLIFY(q%data_tail%prev)

                endif
            endif

            item = q%data_tail%data(q%tail)
            !print *, 'get index: ',q%data_tail%block_id,':', q%tail, 'value: ', item%value

            q%tail = q%tail + 1
            q%item_count = q%item_count - 1
        else
            print *, 'Queue empty'
            stop
        endif

    end function

    !------------------------- Dequeue ----------------------------
    ! Syntax:       xi = q%Dequeue()
    !
    ! Outputs:      xi is a numeric value
    !
    ! Description:  Removes tail item from Q LIFO
    !--------------------------------------------------------------
    function pop(q) result(item)
        implicit none
        class(queue_type) :: q
        type(point_type) :: item

        if(q%item_count > 0) then
            if(q%head < 1) then
                if(associated(q%data_head%prev)) then
                    q%data_head => q%data_head%prev

                    ! set head to the last item on previous block
                    q%head = block_size

                    ! in order to prevent too many allocations
                    ! save one block value to tmp for reuse
                    if (associated(q%tmp)) then
                        deallocate(q%data_head%next)
                    else
                        q%tmp => q%data_head%next
                    endif
                    NULLIFY(q%data_head%next)
                endif
            endif

            item = q%data_head%data(q%head)

            q%head = q%head - 1
            q%item_count = q%item_count - 1


        else
            print *, 'Queue empty'
            stop
        endif

    end function


end module dta_queue
