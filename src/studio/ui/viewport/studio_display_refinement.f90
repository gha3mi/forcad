module studio_display_refinement
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QTimer,QThread_create,QThread_isInterruptionRequested
    use studio_state, only: studio_t,dp,scene_patch_drawn,viewport_fast_shaded,viewport_wireframe
    use studio_display_cache, only: display_cache_entry,display_sample_job,display_view
    use studio_edit_history, only: document_undo_stamp
    use studio_scene, only: studio_scene_count
    use studio_selection, only: selection_target_t
    use studio_display_sampling_policy, only: display_sampling_error,display_sample_counts
    use studio_single_patch, only: copy_single_patch_geometry,single_patch_revision
    use studio_forcad, only: rk,sample_display_grid
    use studio_viewport_camera, only: viewport_camera,prepare_camera,viewport_zoom_active,viewport_box_visible
    use studio_viewport_mesh, only: mesh_covered_patches
    implicit none(type,external)
    private

    public :: connect_display_refinement,request_display_refinement,finish_display_refinement
    public :: display_refinement_pending,display_grid_current

    integer,parameter :: settle_ms=90,worker_poll_ms=24,assessment_yield_ms=4
    real(dp),parameter :: minimum_adaptive_zoom=1.05_dp

contains

    subroutine connect_display_refinement(s)
        type(studio_t),intent(inout),target :: s

        if(.not.s%display_cache%timer%isNull())return
        s%display_cache%timer=QTimer(s%surface)
        call s%display_cache%timer%setSingleShot(.true.)
        call s%display_cache%timer%timeout%connect(refine_display,s)
    end subroutine connect_display_refinement

    subroutine request_display_refinement(s)
        type(studio_t),intent(inout) :: s
        type(viewport_camera) :: camera

        if(s%display_cache%timer%isNull())return
        if(refinement_paused(s))then
            if(refinement_available(s).and.viewport_zoom_active(s).and..not.s%display_cache%timer%isActive())then
                call s%display_cache%timer%start(settle_ms)
            end if
            return
        end if
        if(.not.prepare_refinement_view(s,camera))return
        if(.not.update_display_view(s,camera))return
        call s%display_cache%timer%start(settle_ms)
    end subroutine request_display_refinement

    logical function refinement_available(s) result(available)
        type(studio_t),intent(in) :: s

        available=.not.(s%home_active.or.s%mesh_only.or.s%zoom<=minimum_adaptive_zoom.or. &
            any(s%viewport_style==[viewport_fast_shaded,viewport_wireframe]))
    end function refinement_available

    logical function prepare_refinement_view(s,camera) result(ready)
        type(studio_t),intent(inout) :: s
        type(viewport_camera),intent(out) :: camera

        ready=refinement_available(s)
        if(ready)ready=prepare_camera(s,camera)
        if(.not.ready)s%display_cache%view=display_view()
    end function prepare_refinement_view

    logical function refinement_paused(s) result(paused)
        type(studio_t),intent(in) :: s

        paused=s%dragging.or.s%control_dragging.or.s%updating.or.viewport_zoom_active(s)
    end function refinement_paused

    logical function update_display_view(s,camera) result(changed)
        type(studio_t),intent(inout) :: s
        type(viewport_camera),intent(in) :: camera
        type(display_view) :: incoming
        integer :: index

        incoming%revisions=[s%document%document_revision,document_undo_stamp(s%document),s%document%mesh_generation]
        incoming%projection=[camera%center,camera%screen,camera%scale, &
            camera%yaw_sin,camera%yaw_cos,camera%pitch_sin,camera%pitch_cos]
        incoming%clip=[0._dp,0._dp,real(s%surface%width(),dp),real(s%surface%height(),dp)]
        incoming%mesh_options=[s%mesh_view,s%mesh_section_axis,merge(1,0,s%mesh_shaded)]
        incoming%visible_count=studio_scene_count(s%document%scene)

        associate(cache=>s%display_cache,view=>s%display_cache%view)
            changed=any(view%revisions/=incoming%revisions).or.any(view%projection/=incoming%projection).or. &
                any(view%clip/=incoming%clip).or.any(view%mesh_options/=incoming%mesh_options).or. &
                view%visible_count/=incoming%visible_count.or.view%isolation%count()/=s%isolation%count()
            if(.not.changed)then
                do index=1,s%isolation%count()
                    if(same_target(view%isolation%items(index),s%isolation%items(index)))cycle
                    changed=.true.
                    exit
                end do
            end if
            if(.not.changed)return
            incoming%isolation=s%isolation
            view=incoming
            cache%next_rank=1
            cache%next_variant=1
        end associate
    end function update_display_view

    logical function same_target(first,second) result(same)
        type(selection_target_t),intent(in) :: first,second

        same=first%category==second%category.and.first%rank==second%rank.and.first%id==second%id.and. &
            first%part==second%part.and.first%item==second%item.and.first%feature==second%feature
    end function same_target

    logical function display_refinement_pending(s) result(pending)
        type(studio_t),intent(in) :: s

        pending=associated(s%display_cache%job)
        if(.not.s%display_cache%timer%isNull())pending=pending.or.s%display_cache%timer%isActive()
    end function display_refinement_pending

    logical function display_grid_current(s,kind,variant) result(current)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: kind,variant

        current=.false.
        if(.not.refinement_available(s).or.s%display_cache%document/=s%document%document_revision)return
        associate(rank=>s%display_cache%ranks(kind))
            if(.not.allocated(rank%entries))return
            if(variant>size(rank%entries))return
            associate(entry=>rank%entries(variant))
                if(.not.allocated(entry%grids(1)%points))return
                current=entry%revision==single_patch_revision(kind,variant,catalog=s%document%patches)
                if(current)current=all(entry%base==s%document%scene%catalog(kind)%entries(variant)%patch%sample_counts)
            end associate
        end associate
    end function display_grid_current

    subroutine refine_display(data)
        class(*),intent(inout),optional :: data
        type(viewport_camera) :: camera
        integer :: kind,variant,counts(3)
        logical :: changed,deferred

        select type(s=>data)
        type is(studio_t)
            if(refinement_paused(s))then
                call s%display_cache%timer%start(worker_poll_ms)
                return
            end if
            if(associated(s%display_cache%job))then
                if(s%display_cache%worker%isRunning())then
                    call s%display_cache%timer%start(worker_poll_ms)
                    return
                end if
                call collect_samples(s,changed)
                if(changed)s%display_cache%dirty=.true.
            end if
            if(.not.prepare_refinement_view(s,camera))return
            if(update_display_view(s,camera))then
                call s%display_cache%timer%start(settle_ms)
                return
            end if
            call find_refinement(s,camera,kind,variant,counts,deferred)
            if(deferred)then
                call s%display_cache%timer%start(assessment_yield_ms)
            else if(kind>0)then
                call start_samples(s,kind,variant,counts)
            else if(s%display_cache%dirty)then
                s%display_cache%dirty=.false.
                call s%surface%update()
            end if
        end select
    end subroutine refine_display

    subroutine find_refinement(s,camera,kind,variant,counts,deferred)
        type(studio_t),intent(inout) :: s
        type(viewport_camera),intent(in) :: camera
        integer,intent(out) :: kind,variant,counts(3)
        logical,intent(out) :: deferred
        logical :: assessment_done
        integer :: rank,id,count,grid,first
        integer(int64) :: revision
        type(display_cache_entry),allocatable :: resized(:)
        logical,allocatable :: covered(:,:)

        kind=0; variant=0; counts=1
        assessment_done=.false.; deferred=.false.
        call mesh_covered_patches(s,covered)
        associate(cache=>s%display_cache)
            if(cache%document/=s%document%document_revision)then
                do rank=1,3
                    if(allocated(cache%ranks(rank)%entries))deallocate(cache%ranks(rank)%entries)
                end do
                cache%document=s%document%document_revision
            end if
            do rank=cache%next_rank,3
                if(.not.allocated(s%document%scene%catalog(rank)%entries))cycle
                count=size(s%document%scene%catalog(rank)%entries)
                first=1
                if(rank==cache%next_rank)first=cache%next_variant
                associate(list=>cache%ranks(rank))
                    if(.not.allocated(list%entries))allocate(list%entries(count))
                    if(size(list%entries)<count)then
                        allocate(resized(count))
                        do id=1,size(list%entries)
                            resized(id)%revision=list%entries(id)%revision
                            resized(id)%base=list%entries(id)%base
                            resized(id)%counts=list%entries(id)%counts
                            resized(id)%attempted=list%entries(id)%attempted
                            resized(id)%error=list%entries(id)%error
                            resized(id)%needs_pilot=list%entries(id)%needs_pilot
                            do grid=1,6
                                call move_alloc(list%entries(id)%grids(grid)%points,resized(id)%grids(grid)%points)
                            end do
                        end do
                        call move_alloc(resized,list%entries)
                    end if
                    do id=first,count
                        if(.not.scene_patch_drawn(s,rank,id))cycle
                        if(covered(id,rank))cycle
                        associate(patch=>s%document%scene%catalog(rank)%entries(id)%patch,entry=>list%entries(id))
                            if(rank==1.and.patch%degree(1)<2.and.patch%control_counts(1)<=2)cycle
                            if(any(patch%sample_counts(:rank)<2))cycle
                            if(.not.viewport_box_visible(camera,patch%controls,cache%view%clip))cycle
                            cache%next_rank=rank
                            cache%next_variant=id
                            revision=single_patch_revision(rank,id,catalog=s%document%patches)
                            if(entry%revision/=revision.or.any(entry%base/=patch%sample_counts))then
                                if(assessment_done)then
                                    deferred=.true.
                                    return
                                end if
                                assessment_done=.true.
                                entry=display_cache_entry()
                                entry%revision=revision
                                entry%base=patch%sample_counts
                                entry%counts=entry%base
                                entry%error=display_sampling_error(patch%points,entry%base,rank)
                                entry%needs_pilot=any(entry%base(:rank)==2)
                            end if
                            if(entry%needs_pilot)then
                                counts(:rank)=max(5,min(9,entry%base(:rank)))
                                kind=rank; variant=id
                                return
                            end if
                            counts=display_sample_counts(entry%base,rank,entry%error,camera%scale,entry%counts)
                            if(all(counts==entry%base))then
                                if(allocated(entry%grids(1)%points))then
                                    do grid=1,6
                                        if(allocated(entry%grids(grid)%points))deallocate(entry%grids(grid)%points)
                                    end do
                                    entry%counts=entry%base
                                    cache%dirty=.true.
                                end if
                                entry%attempted=0
                                cycle
                            end if
                            if(all(counts==entry%counts).or.all(counts==entry%attempted))cycle
                            kind=rank; variant=id
                            return
                        end associate
                    end do
                end associate
            end do
        end associate
    end subroutine find_refinement

    subroutine start_samples(s,kind,variant,counts)
        type(studio_t),intent(inout) :: s
        integer,intent(in) :: kind,variant,counts(3)

        associate(cache=>s%display_cache,entry=>s%display_cache%ranks(kind)%entries(variant))
            allocate(cache%job)
            cache%job%kind=kind
            cache%job%variant=variant
            cache%job%base=entry%base
            cache%job%counts=counts
            cache%job%document=cache%document
            cache%job%revision=entry%revision
            cache%job%bounds=s%document%scene%catalog(kind)%entries(variant)%patch%bounds
            cache%job%pilot=entry%needs_pilot
            if(entry%needs_pilot)then
                entry%needs_pilot=.false.
            else
                entry%attempted=counts
            end if
            call copy_single_patch_geometry(kind,variant,cache%job%geometry,catalog=s%document%patches)
            if(.not.allocated(cache%job%geometry))then
                deallocate(cache%job)
                call cache%timer%start(worker_poll_ms)
                return
            end if
            cache%worker=QThread_create(sample_display,cache%job)
            if(cache%worker%isNull())then
                deallocate(cache%job)
            else
                call cache%worker%start()
            end if
            call cache%timer%start(worker_poll_ms)
        end associate
    end subroutine start_samples

    subroutine sample_display(data)
        class(*),intent(inout),optional :: data
        real(rk),allocatable :: points(:,:)
        integer :: side

        select type(job=>data)
        type is(display_sample_job)
            do side=1,merge(6,1,job%kind==3)
                if(QThread_isInterruptionRequested())return
                if(.not.sample_display_grid(job%geometry,job%counts,real(job%bounds,rk),side,points))return
                job%grids(side)%points=real(points,dp)
                deallocate(points)
            end do
            if(job%pilot)call measure_pilot_grid(job)
            job%valid=.true.
        end select
    end subroutine sample_display

    subroutine measure_pilot_grid(job)
        type(display_sample_job),intent(inout) :: job
        real(dp) :: deviation(3),spacing_ratio(3)
        integer :: axes(2),counts(3),side,axis

        if(job%kind==3)then
            do side=1,6
                axis=(side+1)/2
                axes=[mod(axis,3)+1,mod(axis+1,3)+1]
                counts=[job%counts(axes),1]
                deviation=display_sampling_error(job%grids(side)%points,counts,2)
                job%error(axes)=max(job%error(axes),deviation(:2))
            end do
        else
            job%error=display_sampling_error(job%grids(1)%points,job%counts,job%kind)
        end if
        spacing_ratio=real(job%counts-1,dp)/real(max(1,job%base-1),dp)
        job%error=job%error*spacing_ratio**2
    end subroutine measure_pilot_grid

    logical function display_source_matches(s,job) result(matches)
        type(studio_t),intent(in) :: s
        type(display_sample_job),intent(in) :: job

        matches=.false.
        if(job%document/=s%document%document_revision.or.job%document/=s%display_cache%document)return
        associate(patches=>s%document%scene%catalog(job%kind))
            if(.not.allocated(patches%entries))return
            if(job%variant>size(patches%entries))return
            if(single_patch_revision(job%kind,job%variant,catalog=s%document%patches)/=job%revision)return
            matches=all(job%base==patches%entries(job%variant)%patch%sample_counts)
        end associate
    end function display_source_matches

    subroutine collect_samples(s,changed)
        type(studio_t),intent(inout) :: s
        logical,intent(out) :: changed
        integer :: rank,id,side,wanted(3)
        type(viewport_camera) :: camera

        changed=.false.
        associate(cache=>s%display_cache)
            call cache%worker%wait()
            call cache%worker%delete()
            associate(job=>cache%job)
                rank=job%kind; id=job%variant
                if(display_source_matches(s,job))then
                    associate(entry=>cache%ranks(rank)%entries(id))
                        if(job%pilot)then
                            if(job%valid)entry%error=max(entry%error,job%error)
                        else if(refinement_available(s).and.prepare_camera(s,camera))then
                            wanted=display_sample_counts(entry%base,rank,entry%error,camera%scale,entry%counts)
                            if(all(wanted==job%counts))then
                                if(job%valid)then
                                    do side=1,merge(6,1,rank==3)
                                        call move_alloc(job%grids(side)%points,entry%grids(side)%points)
                                    end do
                                    entry%counts=job%counts
                                    changed=.true.
                                end if
                            else
                                entry%attempted=0
                            end if
                        else
                            entry%attempted=0
                        end if
                    end associate
                end if
            end associate
            deallocate(cache%job)
        end associate
    end subroutine collect_samples

    subroutine finish_display_refinement(s)
        type(studio_t),intent(inout) :: s

        associate(cache=>s%display_cache)
            if(.not.cache%timer%isNull())then
                call cache%timer%stop()
                call cache%timer%delete()
            end if
            if(.not.associated(cache%job))return
            if(.not.cache%worker%isNull())then
                call cache%worker%requestInterruption()
                call cache%worker%wait()
                call cache%worker%delete()
            end if
            deallocate(cache%job)
        end associate
    end subroutine finish_display_refinement

end module studio_display_refinement
