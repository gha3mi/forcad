module studio_cli
    use studio_presentation_text, only: studio_product_name
    use studio_exchange, only: studio_exchange_suffix
    use studio_ui, only: studio_default_width, studio_default_height, studio_min_width, studio_min_height, studio_min_size_message
    use studio_preview, only: validate_preview
    implicit none(type, external)
    private
    public :: launch_options,read_launch_options

    type :: launch_options
        integer :: width=studio_default_width,height=studio_default_height
        logical :: fixed_size=.false.,home=.true.,enable_mcp=.false.,remember_workspace=.true.
        character(:),allocatable :: preview,screenshot,project_path
    end type launch_options
contains
    subroutine read_launch_options(options)
        type(launch_options),intent(out) :: options
        character(4096) :: argument,value
        integer :: index,count,status

        count=command_argument_count()
        options%remember_workspace=count==0
        index=1
        do while(index<=count)
            call get_command_argument(index,argument)
            select case(trim(argument))
            case('--screenshot')
                index=index+1
                if(index>count)error stop '--screenshot requires a filename'
                call get_command_argument(index,value)
                options%screenshot=trim(value)
            case('--project')
                index=index+1
                if(index>count)error stop '--project requires a filename'
                call get_command_argument(index,value)
                options%project_path=trim(value)
            case('--preview')
                index=index+1
                if(index>count)error stop '--preview requires a state name'
                call get_command_argument(index,value)
                call validate_preview(trim(value))
                options%preview=trim(value)
            case('--size')
                if(index+2>count)error stop '--size requires WIDTH HEIGHT'
                call get_command_argument(index+1,value)
                read(value,*,iostat=status)options%width
                if(status/=0)error stop 'Invalid width'
                call get_command_argument(index+2,value)
                read(value,*,iostat=status)options%height
                if(status/=0)error stop 'Invalid height'
                if(options%width<studio_min_width.or.options%height<studio_min_height)error stop studio_min_size_message
                options%fixed_size=.true.
                index=index+2
            case('--mcp')
                options%enable_mcp=.true.
            case('--fresh')
            case('--smoke-test')
                error stop 'Run the GUI smoke test with pixi run smoke'
            case('--help','-h')
                call print_usage()
                stop
            case default
                if(argument(1:1)/='-'.and.studio_exchange_suffix(trim(argument))=='.fcad')then
                    options%project_path=trim(argument)
                end if
            end select
            index=index+1
        end do

        if(allocated(options%project_path).and.allocated(options%preview))then
            error stop '--project cannot be combined with --preview'
        end if
        if(allocated(options%preview))then
            options%home=options%preview=='home'
        else
            options%preview='default'
            options%home=.not.allocated(options%screenshot)
        end if
        if(allocated(options%project_path))options%home=.true.
    end subroutine read_launch_options

    subroutine print_usage()
        print '(a)', studio_product_name//': [--screenshot FILE] [--size WIDTH HEIGHT]'
        print '(a)', 'Open a project: --project FILE or pass a .fcad filename directly'
        print '(a)', 'Home and new projects: --preview home|empty|units'
        print '(a)', 'Chess sample: --preview chess|chess-close|chess-parameters|canvas'
        print '(a)', 'Chess views: --preview top|wireframe|clean|inspector|assistant'
        print '(a)', 'Chess ribbon views: --preview ribbon-project|ribbon-create|ribbon-model'
        print '(a)', 'Chess ribbon views: --preview ribbon-parameters|ribbon-assembly|ribbon-view'
        print '(a)', 'Find Command dialog: --preview commands'
        print '(a)', 'Run the native GUI smoke test with pixi run smoke.'
    end subroutine print_usage
end module studio_cli
