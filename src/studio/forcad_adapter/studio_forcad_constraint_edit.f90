module studio_forcad_constraint_edit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk
    use forcad_utils, only: disjoint_set_union, disjoint_set_map, solve_least_squares_qr
    use studio_system_memory, only: working_memory_budget=>mesh_memory_budget
    implicit none(type, external)
    private

    public :: complete_constraint_edit

    character(*), parameter :: incompatible='The requested controls conflict with the existing continuity constraints.'
    character(*), parameter :: insufficient_memory='There is not enough available memory to solve the affected continuity constraints.'

contains

    subroutine complete_constraint_edit(rowptr,columns,values,pinned,delta,error)
        integer, intent(in) :: rowptr(:),columns(:)
        real(rk), intent(in) :: values(:)
        logical, intent(in) :: pinned(:)
        real(rk), intent(inout) :: delta(:,:)
        character(:), allocatable, intent(out) :: error
        integer, allocatable :: rows(:),indices(:),forest(:),component(:),row_column(:),row_next(:),control_next(:)
        integer, allocatable :: first_row(:),first_control(:)
        real(rk), allocatable :: coefficients(:),scales(:),remaining(:),balance(:)
        logical, allocatable :: active(:),driven(:)
        integer :: ncontrols,nrows,i,k,first,column,group,status

        error='ForCAD returned an invalid continuity system.'
        ncontrols=size(pinned)
        nrows=size(rowptr)-1
        if(nrows<0.or.ncontrols/=size(delta,1).or.size(delta,2)<1.or.size(columns)/=size(values))return
        if(rowptr(1)/=1.or.rowptr(nrows+1)<1)return
        if(rowptr(nrows+1)-1/=size(columns))return
        if(any(rowptr(2:)<rowptr(:nrows)).or.any(columns<1).or.any(columns>ncontrols))return
        if(.not.all(ieee_is_finite(values)).or..not.all(ieee_is_finite(delta)))return
        error=''
        if(nrows==0)return

        error=insufficient_memory
        allocate(rows(nrows+1),indices(size(columns)),coefficients(size(values)),forest(ncontrols),active(ncontrols), &
                 row_column(nrows),row_next(nrows),control_next(ncontrols),scales(size(delta,2)), &
                 remaining(size(delta,2)),balance(size(delta,2)),stat=status)
        if(status/=0)return
        if(.not.coalesce_constraints(rowptr,columns,values,ncontrols,rows,indices,coefficients))return
        forest=[(i,i=1,ncontrols)]
        active=.false.
        row_column=0
        do i=1,nrows
            first=0
            do k=rows(i),rows(i+1)-1
                column=indices(k)
                if(pinned(column))cycle
                active(column)=.true.
                if(first==0)first=column
                call disjoint_set_union(forest,first,column)
            end do
            row_column(i)=first
        end do
        call disjoint_set_map(forest,component,active=active)
        group=max(0,maxval(component))
        allocate(first_row(group),first_control(group),driven(group),stat=status)
        if(status/=0)return
        first_row=0
        first_control=0
        driven=.false.
        do i=ncontrols,1,-1
            group=component(i)
            if(group==0)cycle
            control_next(i)=first_control(group)
            first_control(group)=i
        end do
        do i=nrows,1,-1
            first=row_column(i)
            if(first==0)then
                scales=0._rk
                do k=rows(i),rows(i+1)-1
                    scales=max(scales,abs(delta(indices(k),:)))
                end do
                scales=merge(scales,1._rk,scales>0._rk)
                remaining=0._rk
                balance=0._rk
                do k=rows(i),rows(i+1)-1
                    remaining=remaining+coefficients(k)*(delta(indices(k),:)/scales)
                    balance=balance+abs(coefficients(k))*(abs(delta(indices(k),:))/scales)
                end do
                if(any(abs(remaining)>equation_tolerance(rows(i+1)-rows(i))*balance))then
                    error=incompatible
                    return
                end if
                cycle
            end if
            group=component(first)
            row_next(i)=first_row(group)
            first_row(group)=i
            if(driven(group))cycle
            do k=rows(i),rows(i+1)-1
                column=indices(k)
                if(.not.pinned(column))cycle
                driven(group)=any(delta(column,:)/=0._rk)
                if(driven(group))exit
            end do
        end do
        do group=1,size(driven)
            if(.not.driven(group))cycle
            call solve_component(rows,indices,coefficients,pinned,delta,first_row(group),row_next, &
                                 first_control(group),control_next,error)
            if(len(error)>0)return
        end do
        error=''
    end subroutine complete_constraint_edit

    logical function coalesce_constraints(rowptr,columns,values,ncontrols,rows,indices,coefficients) result(ready)
        integer, intent(in) :: rowptr(:),columns(:),ncontrols
        real(rk), intent(in) :: values(:)
        integer, intent(out) :: rows(:),indices(:)
        real(rk), intent(out) :: coefficients(:)
        integer, allocatable :: slots(:),terms(:)
        real(rk), allocatable :: sums(:),magnitudes(:)
        integer :: i,k,lo,hi,n,column,position,out,status
        real(rk) :: largest,value

        ready=.false.
        allocate(slots(ncontrols),terms(ncontrols),sums(ncontrols),magnitudes(ncontrols),stat=status)
        if(status/=0)return
        slots=0
        out=1
        rows(1)=out
        do i=1,size(rowptr)-1
            lo=rowptr(i)
            hi=rowptr(i+1)-1
            largest=maxval(abs(values(lo:hi)))
            n=0
            if(largest>0._rk)then
                do k=lo,hi
                    column=columns(k)
                    if(values(k)==0._rk)cycle
                    value=values(k)/largest
                    position=slots(column)
                    if(position==0)then
                        n=n+1
                        position=n
                        slots(column)=position
                        terms(position)=column
                        sums(position)=0._rk
                        magnitudes(position)=0._rk
                    end if
                    sums(position)=sums(position)+value
                    magnitudes(position)=magnitudes(position)+abs(value)
                end do
            end if
            largest=0._rk
            do k=1,n
                slots(terms(k))=0
                if(abs(sums(k))<=equation_tolerance(hi-lo+1)*magnitudes(k))cycle
                indices(out)=terms(k)
                coefficients(out)=sums(k)
                largest=max(largest,abs(sums(k)))
                out=out+1
            end do
            if(out>rows(i))coefficients(rows(i):out-1)=coefficients(rows(i):out-1)/largest
            rows(i+1)=out
        end do
        ready=.true.
    end function coalesce_constraints

    subroutine solve_component(rows,indices,coefficients,pinned,delta,row_head,row_next,control_head,control_next,error)
        integer, intent(in) :: rows(:),indices(:),row_head,row_next(:),control_head,control_next(:)
        real(rk), intent(in) :: coefficients(:)
        logical, intent(in) :: pinned(:)
        real(rk), intent(inout) :: delta(:,:)
        character(:), allocatable, intent(out) :: error
        integer, allocatable :: local(:),members(:)
        real(rk), allocatable :: matrix(:,:),rhs(:,:),scales(:),correction(:,:),residual(:,:),balance(:,:)
        real(rk) :: required,tolerance
        integer :: row,nrows,column,ncontrols,nrhs,j,k,status

        nrows=chain_length(row_head,row_next)
        ncontrols=chain_length(control_head,control_next)
        nrhs=size(delta,2)
        required=6._rk*real(nrows,rk)*real(ncontrols,rk)+4._rk*max(nrows,ncontrols)*real(nrhs,rk)
        error=insufficient_memory
        if(required>real(working_memory_budget(),rk)/(storage_size(0._rk)/8))return
        allocate(local(size(pinned)),members(ncontrols),matrix(nrows,ncontrols),rhs(nrows,nrhs),scales(nrhs),stat=status)
        if(status/=0)return
        column=control_head
        do j=1,ncontrols
            local(column)=j
            members(j)=column
            column=control_next(column)
        end do

        scales=0._rk
        row=row_head
        do while(row/=0)
            do k=rows(row),rows(row+1)-1
                column=indices(k)
                if(pinned(column))scales=max(scales,abs(delta(column,:)))
            end do
            row=row_next(row)
        end do
        scales=merge(scales,1._rk,scales>0._rk)
        matrix=0._rk
        rhs=0._rk
        row=row_head
        do j=1,nrows
            do k=rows(row),rows(row+1)-1
                column=indices(k)
                if(pinned(column))then
                    rhs(j,:)=rhs(j,:)-coefficients(k)*(delta(column,:)/scales)
                else
                    matrix(j,local(column))=coefficients(k)
                end if
            end do
            row=row_next(row)
        end do
        tolerance=equation_tolerance(max(nrows,ncontrols))
        error=incompatible
        if(.not.minimum_constraint_correction(matrix,rhs,tolerance,correction))return
        residual=matmul(matrix,correction)-rhs
        balance=matmul(abs(matrix),abs(correction))
        if(.not.all(ieee_is_finite(residual)).or..not.all(ieee_is_finite(balance)))return
        if(any(abs(residual)>tolerance*max(1._rk,abs(rhs),balance)))return
        correction=correction*spread(scales,1,ncontrols)
        if(.not.all(ieee_is_finite(correction)))return
        delta(members,:)=correction
        error=''
    end subroutine solve_component

    logical function minimum_constraint_correction(matrix,rhs,tolerance,correction) result(solved)
        real(rk), intent(in), contiguous :: matrix(:,:),rhs(:,:)
        real(rk), intent(in) :: tolerance
        real(rk), allocatable, intent(out) :: correction(:,:)
        real(rk), allocatable :: basis(:,:),row(:),reduced(:,:),solution(:,:)
        real(rk) :: norm,threshold
        integer :: i,j,pass,rank,nrows,ncontrols,status

        solved=.false.
        nrows=size(matrix,1)
        ncontrols=size(matrix,2)
        allocate(basis(ncontrols,min(nrows,ncontrols)),row(ncontrols),stat=status)
        if(status/=0)return
        threshold=tolerance*maxval(abs(matrix))
        rank=0
        do i=1,nrows
            row=matrix(i,:)
            do pass=1,2
                do j=1,rank
                    row=row-dot_product(row,basis(:,j))*basis(:,j)
                end do
            end do
            norm=norm2(row)
            if(norm<=threshold)cycle
            rank=rank+1
            basis(:,rank)=row/norm
            if(rank==ncontrols)exit
        end do
        if(rank==0)then
            allocate(correction(ncontrols,size(rhs,2)),source=0._rk,stat=status)
            solved=status==0
            return
        end if
        reduced=matmul(matrix,basis(:,:rank))
        solution=solve_least_squares_qr(reduced,rhs)
        if(size(solution,1)/=rank)return
        correction=matmul(basis(:,:rank),solution)
        solved=all(ieee_is_finite(correction))
    end function minimum_constraint_correction

    pure integer function chain_length(first,next) result(count)
        integer, intent(in) :: first,next(:)
        integer :: index

        count=0
        index=first
        do while(index/=0)
            count=count+1
            index=next(index)
        end do
    end function chain_length

    pure real(rk) function equation_tolerance(extent) result(tolerance)
        integer, intent(in) :: extent

        tolerance=min(sqrt(epsilon(1._rk)),256._rk*max(1,extent)*epsilon(1._rk))
    end function equation_tolerance
end module studio_forcad_constraint_edit
