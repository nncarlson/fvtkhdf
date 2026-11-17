! Small, deterministic fixtures for manual reader compatibility testing.
! Usage: generate output.vtkhdf mb|standalone type mode [type mode ...]
! Types: image, rectilinear, structured, unstructured. Modes: 0,1,2,3,
! or fixed-constant (mode 1 with only constant arrays, to check a workaround).
program generate
  use,intrinsic :: iso_fortran_env, only: int8, real64
  use vtkhdf_mb_file_type
  use vtkhdf_image_file_type, only: vtkhdf_image_file
  use vtkhdf_rectilinear_file_type, only: vtkhdf_rectilinear_file
  use vtkhdf_structured_file_type, only: vtkhdf_structured_file
  use vtkhdf_unstructured_file_type, only: vtkhdf_unstructured_file
  use vtkhdf_simple_dataset_class, only: vtkhdf_simple_dataset
  use vtkhdf_structured_base_class, only: vtkhdf_structured_base
  use vtkhdf_image_writer_type, only: vtkhdf_image_writer
  use vtkhdf_rectilinear_writer_type, only: vtkhdf_rectilinear_writer
  use vtkhdf_structured_writer_type, only: vtkhdf_structured_writer
  use vtkhdf_unstructured_writer_type, only: vtkhdf_unstructured_writer
  implicit none

  type :: block_info
    class(vtkhdf_simple_dataset), pointer :: dataset => null()
    integer :: mode
    logical :: constant_only = .false.
    type(vtkhdf_point_data_handle) :: point_scalar, point_vector
    type(vtkhdf_cell_data_handle) :: cell_scalar, cell_vector
    type(vtkhdf_field_data_handle) :: clock, table
  end type
  type(block_info), allocatable :: blocks(:)
  type(vtkhdf_mb_file) :: mb
  type(vtkhdf_image_file), target :: image_file
  type(vtkhdf_rectilinear_file), target :: rect_file
  type(vtkhdf_structured_file), target :: struct_file
  type(vtkhdf_unstructured_file), target :: unstructured_file
  type(vtkhdf_image_block), pointer :: image
  type(vtkhdf_rectilinear_block), pointer :: rect
  type(vtkhdf_structured_block), pointer :: grid
  type(vtkhdf_unstructured_block), pointer :: solid
  real(real64), parameter :: times(0:4) = [0.0_real64,0.2_real64,0.7_real64,1.5_real64,3.0_real64]
  integer, parameter :: origin(3) = [-2,3,1], pdims(3) = [3,4,2], cdims(3) = pdims-1
  integer :: stat, nblocks, b, step, mode
  character(1024) :: filename
  character(32) :: form, kind, argument, name
  character(:), allocatable :: errmsg
  logical :: composite, temporal

  if (command_argument_count() < 4 .or. mod(command_argument_count(),2) /= 0) &
      error stop 'usage: generate output mb|standalone type mode [type mode ...]'
  call get_command_argument(1,filename)
  call get_command_argument(2,form)
  if (form /= 'mb' .and. form /= 'standalone') error stop 'invalid output form'
  composite = form == 'mb'
  nblocks = (command_argument_count()-2)/2
  if (.not.composite .and. nblocks /= 1) error stop 'standalone requires one dataset'
  allocate(blocks(nblocks))
  if (composite) then
    call mb%create(trim(filename),stat,errmsg)
    if (stat /= 0) error stop errmsg
  end if
  do b = 1,nblocks
    call get_command_argument(2*b+1,kind)
    call get_command_argument(2*b+2,argument)
    blocks(b)%constant_only = argument == 'fixed-constant'
    if (blocks(b)%constant_only) then
      mode = VTKHDF_FIXED_MESH
      stat = 0
    else
      read(argument,*,iostat=stat) mode
    end if
    if (stat /= 0) error stop 'invalid mode'
    blocks(b)%mode = mode
    write(name,'(a,i0)') 'block_',b
    select case (kind)
    case ('image')
      if (composite) then
        call mb%add_block(trim(name),image,mode)
        blocks(b)%dataset => image
      else
        call image_file%create(trim(filename),stat,errmsg,mode)
        blocks(b)%dataset => image_file
      end if
    case ('rectilinear')
      if (composite) then
        call mb%add_block(trim(name),rect,mode)
        blocks(b)%dataset => rect
      else
        call rect_file%create(trim(filename),stat,errmsg,mode)
        blocks(b)%dataset => rect_file
      end if
    case ('structured')
      if (composite) then
        call mb%add_block(trim(name),grid,mode)
        blocks(b)%dataset => grid
      else
        call struct_file%create(trim(filename),stat,errmsg,mode)
        blocks(b)%dataset => struct_file
      end if
    case ('unstructured')
      if (composite) then
        call mb%add_block(trim(name),solid,mode)
        blocks(b)%dataset => solid
      else
        call unstructured_file%create(trim(filename),stat,errmsg,mode)
        blocks(b)%dataset => unstructured_file
      end if
    case default
      error stop 'invalid dataset type'
    end select
    if (.not.composite .and. stat /= 0) error stop errmsg
    call mesh(blocks(b),b,0,.true.)
    call setup_data(blocks(b),b)
  end do

  temporal = any(blocks%mode /= VTKHDF_STATIC)
  if (temporal) then
    do step = 0,4
      if (composite) then
        call mb%start_time_step(times(step))
      else
        call blocks(1)%dataset%start_time_step(times(step))
      end if
      do b = 1,nblocks
        if (blocks(b)%mode == VTKHDF_STATIC) cycle
        if (blocks(b)%mode >= VTKHDF_MOVING_MESH .and. (step == 2 .or. step == 4)) &
            call mesh(blocks(b),b,step,.false.)
        if (.not.blocks(b)%constant_only) call temporal_data(blocks(b),b,step)
      end do
      ! START finalizes step 0; CLOSE finalizes step 4.
      if (step == 0 .or. step == 4) cycle
      if (composite) then
        call mb%finalize_time_step()
      else
        call blocks(1)%dataset%finalize_time_step()
      end if
    end do
  end if
  if (composite) then
    call mb%close()
  else
    call blocks(1)%dataset%close()
  end if

contains

  subroutine mesh(block,id,step,initial)
    type(block_info), intent(inout) :: block
    integer, intent(in) :: id, step
    logical, intent(in) :: initial
    integer :: i, j, k, c, nc, epoch
    real(real64) :: xyz(3,3,4,2), tetra(3,4)
    real(real64), allocatable :: points(:,:)
    integer, allocatable :: connectivity(:), offsets(:)
    integer(int8), allocatable :: types(:)
    epoch = step/2
    select type (writer => block%dataset)
    class is (vtkhdf_image_writer)
      call writer%write_mesh(reshape([origin,origin+pdims-1],[3,2]), &
          [10.0_real64*id,-3.0_real64,2.0_real64],[0.5_real64,1.0_real64,2.0_real64])
    class is (vtkhdf_rectilinear_writer)
      if (initial) call writer%write_mesh(reshape([origin,origin+pdims-1],[3,2]))
      call writer%write_mesh_geometry([0.0_real64,0.6_real64,1.4_real64]+10*id, &
          [0.0_real64,1.1_real64,2.4_real64,3.9_real64], [0.0_real64,2.1_real64]+0.25_real64*epoch)
    class is (vtkhdf_structured_writer)
      if (initial) call writer%write_mesh(reshape([origin,origin+pdims-1],[3,2]))
      do k = 0,1
        do j = 0,3
          do i = 0,2
            xyz(:,i+1,j+1,k+1) = [10*id+0.5_real64*i+0.05_real64*j*k, &
                j+0.1_real64*i*k,2*k+0.25_real64*epoch]
          end do
        end do
      end do
      call writer%write_mesh_geometry(xyz,origin)
    class is (vtkhdf_unstructured_writer)
      nc = 1
      if (block%mode == VTKHDF_DYNAMIC_MESH .and. epoch == 1) nc = 2
      allocate(points(3,4*nc))
      tetra = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,4])
      do c = 1,nc
        points(:,4*c-3:4*c) = tetra
        points(1,4*c-3:4*c) = points(1,4*c-3:4*c)+10*id+2*(c-1)
      end do
      points(3,:) = points(3,:)+0.25_real64*epoch
      if (initial .or. block%mode == VTKHDF_DYNAMIC_MESH) then
        connectivity = [(i,i=1,4*nc)]
        offsets = [(4*i+1,i=0,nc)]
        types = spread(10_int8,1,nc)
        call writer%write_mesh(points,connectivity,offsets,types)
      else
        call writer%write_mesh_geometry(points)
      end if
    end select
  end subroutine

  subroutine setup_data(block,id)
    type(block_info), intent(inout) :: block
    integer, intent(in) :: id
    integer :: i, c
    real(real64) :: ps(3,4,2), pv(3,3,4,2), cs(2,3,1), cv(3,2,3,1)
    ps = reshape([(real(1000*id+i,real64),i=1,24)],pdims)
    cs = reshape([(real(1000*id+i,real64),i=1,6)],cdims)
    do c = 1,3
      pv(c,:,:,:) = ps+c
      cv(c,:,:,:) = cs+c
    end do
    select type (writer => block%dataset)
    class is (vtkhdf_structured_base)
      call writer%write_point_data('point_constant',ps,origin)
      call writer%write_point_data('point_constant_vector',pv,origin)
      call writer%write_cell_data('cell_constant',cs,origin)
      call writer%write_cell_data('cell_constant_vector',cv,origin)
      call writer%write_field_data('block_id',data=id)
      call writer%write_field_data('direction',data=[1.0_real64,2.0_real64,3.0_real64],as_vector=.true.)
      if (block%mode /= VTKHDF_STATIC .and. .not.block%constant_only) then
        block%point_scalar = writer%register_temporal_point_data('point_value',0.0_real64)
        block%point_vector = writer%register_temporal_point_data('point_vector',[0.0_real64,0.0_real64,0.0_real64])
        block%cell_scalar = writer%register_temporal_cell_data('cell_value',0.0_real64)
        block%cell_vector = writer%register_temporal_cell_data('cell_vector',[0.0_real64,0.0_real64,0.0_real64])
        block%clock = writer%register_temporal_field_data('clock',0.0_real64)
        block%table = writer%register_temporal_field_data('table',0.0_real64)
      end if
    class is (vtkhdf_unstructured_writer)
      if (block%mode /= VTKHDF_DYNAMIC_MESH) then
        block_scalar: block
          real(real64) :: s(4), v(3,4)
          s = [(real(1000*id+i,real64),i=1,4)]
          do c = 1,3
            v(c,:) = s+c
          end do
          call writer%write_point_data('point_constant',s)
          call writer%write_point_data('point_constant_vector',v)
          call writer%write_cell_data('cell_constant',[real(1000*id+1,real64)])
          call writer%write_cell_data('cell_constant_vector',v(:,:1))
        end block block_scalar
      end if
      call writer%write_field_data('block_id',data=id)
      call writer%write_field_data('direction',data=[1.0_real64,2.0_real64,3.0_real64],as_vector=.true.)
      if (block%mode /= VTKHDF_STATIC .and. .not.block%constant_only) then
        block%point_scalar = writer%register_temporal_point_data('point_value',0.0_real64)
        block%point_vector = writer%register_temporal_point_data('point_vector',[0.0_real64,0.0_real64,0.0_real64])
        block%cell_scalar = writer%register_temporal_cell_data('cell_value',0.0_real64)
        block%cell_vector = writer%register_temporal_cell_data('cell_vector',[0.0_real64,0.0_real64,0.0_real64])
        block%clock = writer%register_temporal_field_data('clock',0.0_real64)
        block%table = writer%register_temporal_field_data('table',0.0_real64)
      end if
    end select
  end subroutine

  subroutine temporal_data(block,id,step)
    type(block_info), intent(inout) :: block
    integer, intent(in) :: id, step
    integer :: i, c, np, nc, nt
    logical :: write_points, write_cells
    real(real64) :: ps(3,4,2), pv(3,3,4,2), cs(2,3,1), cv(3,2,3,1)
    real(real64), allocatable :: s(:), v(:,:), table(:,:)
    write_points = step == 0 .or. step == 1 .or. step == 3
    write_cells = step == 0 .or. step == 2 .or. step == 4
    if (block%mode == VTKHDF_DYNAMIC_MESH .and. (step == 2 .or. step == 4)) write_points = .true.
    ps = reshape([(real(1000*id+100*step+i,real64),i=1,24)],pdims)
    cs = reshape([(real(1000*id+100*step+i,real64),i=1,6)],cdims)
    do c = 1,3
      pv(c,:,:,:) = ps+c
      cv(c,:,:,:) = cs+c
    end do
    nt = 1+mod(step,3)
    table = reshape([(real(1000*id+100*step+i,real64),i=1,2*nt)],[2,nt])
    select type (writer => block%dataset)
    class is (vtkhdf_structured_base)
      if (write_points) then
        call writer%write_point_data(block%point_scalar,ps,origin)
        call writer%write_point_data(block%point_vector,pv,origin)
      end if
      if (write_cells) then
        call writer%write_cell_data(block%cell_scalar,cs,origin)
        call writer%write_cell_data(block%cell_vector,cv,origin)
      end if
      call writer%write_field_data(block%clock,data=times(step))
      call writer%write_field_data(block%table,data=table)
    class is (vtkhdf_unstructured_writer)
      nc = 1
      if (block%mode == VTKHDF_DYNAMIC_MESH .and. step/2 == 1) nc = 2
      np = 4*nc
      allocate(s(np),v(3,np))
      s = [(real(1000*id+100*step+i,real64),i=1,np)]
      do c = 1,3
        v(c,:) = s+c
      end do
      if (write_points) then
        call writer%write_point_data(block%point_scalar,s)
        call writer%write_point_data(block%point_vector,v)
      end if
      if (write_cells) then
        call writer%write_cell_data(block%cell_scalar,s(:nc))
        call writer%write_cell_data(block%cell_vector,v(:,:nc))
      end if
      call writer%write_field_data(block%clock,data=times(step))
      call writer%write_field_data(block%table,data=table)
    end select
  end subroutine
end program
