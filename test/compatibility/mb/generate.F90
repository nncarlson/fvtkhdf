! Small, deterministic fixtures for manual reader compatibility testing.
! Usage: generate output.vtkhdf mb|standalone ndim type mode [type mode ...]
! Types: image, rectilinear, structured, unstructured, hypertree. Modes: 0,1,2,3,
! or fixed-constant (mode 1 with only constant arrays, to check a workaround).
program generate
  use,intrinsic :: iso_fortran_env, only: int8, real64
  use vtkhdf_mb_file_type
  use vtkhdf_image_file_type, only: vtkhdf_image_file
  use vtkhdf_rectilinear_file_type, only: vtkhdf_rectilinear_file
  use vtkhdf_structured_file_type, only: vtkhdf_structured_file
  use vtkhdf_unstructured_file_type, only: vtkhdf_unstructured_file
  use vtkhdf_hypertree_file_type, only: vtkhdf_hypertree_file
  use vtkhdf_simple_dataset_class, only: vtkhdf_simple_dataset
  use vtkhdf_structured_base_class, only: vtkhdf_structured_base
  use vtkhdf_image_writer_type, only: vtkhdf_image_writer
  use vtkhdf_rectilinear_writer_type, only: vtkhdf_rectilinear_writer
  use vtkhdf_structured_writer_type, only: vtkhdf_structured_writer
  use vtkhdf_unstructured_writer_type, only: vtkhdf_unstructured_writer
  use vtkhdf_hypertree_writer_type, only: vtkhdf_hypertree_writer
#ifdef USE_MPI_F08
  use mpi_f08
#elif defined(USE_MPI)
  use mpi
#endif
  implicit none

  type :: block_info
    class(vtkhdf_simple_dataset), pointer :: dataset => null()
    integer :: mode
    logical :: constant_only = .false.
    integer :: pdims(3), cdims(3), local_points(3), local_cells(3)
    integer :: array_origin(3), selection(3,2)
    integer, allocatable :: tree_entities(:)
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
  type(vtkhdf_hypertree_file), target :: hypertree_file
  type(vtkhdf_image_block), pointer :: image
  type(vtkhdf_rectilinear_block), pointer :: rect
  type(vtkhdf_structured_block), pointer :: grid
  type(vtkhdf_unstructured_block), pointer :: solid
  type(vtkhdf_hypertree_block), pointer :: trees
  real(real64), parameter :: times(0:4) = [0.0_real64,0.2_real64,0.7_real64,1.5_real64,3.0_real64]
  integer, parameter :: origin(3) = [-2,3,1]
  integer :: stat, nblocks, b, step, mode, ndim, rank, nproc
  character(1024) :: filename
  character(32) :: form, kind, argument, name
  character(:), allocatable :: errmsg
  logical :: composite, temporal

  call get_command_argument(1,argument)
  if (argument == '--info') then
#ifdef USE_MPI
    print '(a)', '{"mpi":true}'
#else
    print '(a)', '{"mpi":false}'
#endif
    stop
  end if
  if (command_argument_count() < 5 .or. mod(command_argument_count()-3,2) /= 0) &
      error stop 'usage: generate output mb|standalone ndim type mode [type mode ...]'
#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,stat)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,stat)
#else
  rank = 0
  nproc = 1
#endif
  call get_command_argument(1,filename)
  call get_command_argument(2,form)
  if (form /= 'mb' .and. form /= 'standalone') error stop 'invalid output form'
  composite = form == 'mb'
  call get_command_argument(3,argument)
  read(argument,*,iostat=stat) ndim
  if (stat /= 0) error stop 'invalid dimensionality'
  if (ndim < 1 .or. ndim > 3) error stop 'dimensionality must be 1, 2, or 3'
  nblocks = (command_argument_count()-3)/2
  if (.not.composite .and. nblocks /= 1) error stop 'standalone requires one dataset'
  allocate(blocks(nblocks))
  if (composite) then
    call mb%create(trim(filename), &
#ifdef USE_MPI
        MPI_COMM_WORLD, &
#endif
        stat,errmsg)
    if (stat /= 0) error stop errmsg
  end if
  do b = 1,nblocks
    call get_command_argument(2*b+2,kind)
    call get_command_argument(2*b+3,argument)
    blocks(b)%constant_only = argument == 'fixed-constant'
    if (blocks(b)%constant_only) then
      mode = VTKHDF_FIXED_MESH
      stat = 0
    else
      read(argument,*,iostat=stat) mode
    end if
    if (stat /= 0) error stop 'invalid mode'
    blocks(b)%mode = mode
    call structured_partition(blocks(b),b)
    write(name,'(a,i0)') 'block_',b
    select case (kind)
    case ('image')
      if (composite) then
        call mb%add_block(trim(name),image,mode)
        blocks(b)%dataset => image
      else
        call image_file%create(trim(filename), &
#ifdef USE_MPI
            MPI_COMM_WORLD, &
#endif
            stat,errmsg,mode)
        blocks(b)%dataset => image_file
      end if
    case ('rectilinear')
      if (composite) then
        call mb%add_block(trim(name),rect,mode)
        blocks(b)%dataset => rect
      else
        call rect_file%create(trim(filename), &
#ifdef USE_MPI
            MPI_COMM_WORLD, &
#endif
            stat,errmsg,mode)
        blocks(b)%dataset => rect_file
      end if
    case ('structured')
      if (composite) then
        call mb%add_block(trim(name),grid,mode)
        blocks(b)%dataset => grid
      else
        call struct_file%create(trim(filename), &
#ifdef USE_MPI
            MPI_COMM_WORLD, &
#endif
            stat,errmsg,mode)
        blocks(b)%dataset => struct_file
      end if
    case ('unstructured')
      if (composite) then
        call mb%add_block(trim(name),solid,mode)
        blocks(b)%dataset => solid
      else
        call unstructured_file%create(trim(filename), &
#ifdef USE_MPI
            MPI_COMM_WORLD, &
#endif
            stat,errmsg,mode)
        blocks(b)%dataset => unstructured_file
      end if
    case ('hypertree')
      if (composite) then
        call mb%add_block(trim(name),trees,mode)
        blocks(b)%dataset => trees
      else
        call hypertree_file%create(trim(filename), &
#ifdef USE_MPI
            MPI_COMM_WORLD, &
#endif
            stat,errmsg,mode)
        blocks(b)%dataset => hypertree_file
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
#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

contains

  ! Each block has a different empty rank. Remaining ranks contribute unequal
  ! contiguous slabs/parts. Rank zero can be empty while still writing metadata.
  subroutine partition(id,phase,nlocal,start,total,owner)
    integer, intent(in) :: id, phase
    integer, intent(in), optional :: owner
    integer, intent(out) :: nlocal, start, total
    integer :: r, count, owner_
    owner_ = rank
    if (present(owner)) owner_ = owner
    nlocal = 0
    start = 0
    total = 0
    do r = 0,nproc-1
      count = r+1+phase
      if (nproc > 1 .and. r == mod(id-1,nproc)) count = 0
      if (r < owner_) start = start+count
      if (r == owner_) nlocal = count
      total = total+count
    end do
  end subroutine

  subroutine structured_partition(block,id)
    type(block_info), intent(inout) :: block
    integer, intent(in) :: id
    integer :: nlocal, start, total
    block%pdims = [3,4,2]
    if (ndim < 3) block%pdims(3) = 1
    if (ndim < 2) block%pdims(2) = 1
    call partition(id,0,nlocal,start,total)
    if (nproc == 1) then
      nlocal = 2
    else
      block%pdims(1) = total+1
    end if
    block%cdims = max(1,block%pdims-1)
    block%local_cells = block%cdims
    block%local_points = block%pdims
    block%local_cells(1) = nlocal
    block%local_points(1) = merge(nlocal+1,0,nlocal>0)
    block%array_origin = origin+[start,0,0]
    block%selection(:,1) = block%array_origin
    block%selection(:,2) = block%array_origin+block%local_points-1
    if (nlocal > 0 .and. start > 0) block%selection(1,1) = block%selection(1,1)+1
  end subroutine

  subroutine structured_values(block,id,step,ps,pv,cs,cv)
    type(block_info), intent(in) :: block
    integer, intent(in) :: id, step
    real(real64), allocatable, intent(out) :: ps(:,:,:), pv(:,:,:,:), cs(:,:,:), cv(:,:,:,:)
    integer :: i, j, k, c, start
    allocate(ps(block%local_points(1),block%local_points(2),block%local_points(3)))
    allocate(cs(block%local_cells(1),block%local_cells(2),block%local_cells(3)))
    allocate(pv(3,size(ps,1),size(ps,2),size(ps,3)),cv(3,size(cs,1),size(cs,2),size(cs,3)))
    start = block%array_origin(1)-origin(1)
    do k = 1,size(ps,3)
      do j = 1,size(ps,2)
        do i = 1,size(ps,1)
          ps(i,j,k) = 1000*id+100*step+start+i+block%pdims(1)*(j-1+block%pdims(2)*(k-1))
        end do
      end do
    end do
    do k = 1,size(cs,3)
      do j = 1,size(cs,2)
        do i = 1,size(cs,1)
          cs(i,j,k) = 1000*id+100*step+start+i+block%cdims(1)*(j-1+block%cdims(2)*(k-1))
        end do
      end do
    end do
    do c = 1,3
      pv(c,:,:,:) = ps+c
      cv(c,:,:,:) = cs+c
    end do
  end subroutine

  subroutine unstructured_values(block,id,step,ps,pv,cs,cv)
    type(block_info), intent(in) :: block
    integer, intent(in) :: id, step
    real(real64), allocatable, intent(out) :: ps(:), pv(:,:), cs(:), cv(:,:)
    integer :: i, c, nlocal, start, total, phase
    phase = 0
    if (block%mode == VTKHDF_DYNAMIC_MESH .and. step/2 == 1) phase = 1
    call partition(id,phase,nlocal,start,total)
    ps = [(real(1000*id+100*step+4*start+i,real64),i=1,4*nlocal)]
    cs = [(real(1000*id+100*step+start+i,real64),i=1,nlocal)]
    allocate(pv(3,size(ps)),cv(3,size(cs)))
    do c = 1,3
      pv(c,:) = ps+c
      cv(c,:) = cs+c
    end do
  end subroutine

  ! Root ownership follows x slabs and rotates between ranks in dynamic mode.
  ! Alternating roots have one refined level; dynamic phases swap that pattern.
  subroutine hypertree_topology(block,id,step,ids,xtree,refined,mask)
    type(block_info), intent(inout) :: block
    integer, intent(in) :: id, step
    integer, allocatable, intent(out) :: ids(:), xtree(:)
    logical, allocatable, intent(out) :: refined(:), mask(:)
    integer :: phase, nc, start, total, g, j, n, nchild
    phase = 0
    if (block%mode == VTKHDF_DYNAMIC_MESH) phase = step/2
    call partition(id,0,nc,start,total,owner=modulo(rank-phase,nproc))
    if (nproc == 1) nc = block%cdims(1)
    ids = [integer::]
    xtree = [1]
    refined = [logical::]
    mask = [logical::]
    block%tree_entities = [integer::]
    nchild = 2**ndim
    do g = 0,product(block%cdims)-1
      if (modulo(g,block%cdims(1)) < start .or. modulo(g,block%cdims(1)) >= start+nc) cycle
      ids = [ids,g+1]
      n = 1
      if (modulo(g+phase,2) == 0) n = 1+nchild
      xtree = [xtree,xtree(size(xtree))+n]
      refined = [refined,(j==0 .and. n>1,j=0,n-1)]
      mask = [mask,(j==n-1 .and. modulo(g,3)==0,j=0,n-1)]
      block%tree_entities = [block%tree_entities,((nchild+1)*g+j+1,j=0,n-1)]
    end do
  end subroutine

  subroutine hypertree_values(block,id,step,cs,cv)
    type(block_info), intent(in) :: block
    integer, intent(in) :: id, step
    real(real64), allocatable, intent(out) :: cs(:), cv(:,:)
    integer :: c
    cs = real(1000*id+100*step+block%tree_entities,real64)
    allocate(cv(3,size(cs)))
    do c = 1,3
      cv(c,:) = cs+c
    end do
  end subroutine

  subroutine mesh(block,id,step,initial)
    type(block_info), intent(inout) :: block
    integer, intent(in) :: id, step
    logical, intent(in) :: initial
    integer :: i, j, k, c, nc, epoch, start, total, phase
    real(real64) :: tetra(3,4)
    real(real64), allocatable :: points(:,:), xyz(:,:,:,:), x(:), y(:), z(:)
    integer, allocatable :: connectivity(:), offsets(:)
    integer, allocatable :: ids(:), xtree(:)
    logical, allocatable :: refined(:), mask(:)
    integer(int8), allocatable :: types(:)
    epoch = step/2
    select type (writer => block%dataset)
    class is (vtkhdf_image_writer)
      call writer%write_mesh(reshape([origin,origin+block%pdims-1],[3,2]), &
          [10.0_real64*id,-3.0_real64,2.0_real64],[0.5_real64,1.0_real64,2.0_real64])
    class is (vtkhdf_rectilinear_writer)
      if (initial) call writer%write_mesh(reshape([origin,origin+block%pdims-1],[3,2]))
      x = [(10*id+0.5_real64*i+0.1_real64*i*i,i=0,block%pdims(1)-1)]
      y = [(j+0.1_real64*j*j,j=0,block%pdims(2)-1)]
      z = [(2*k+0.1_real64*k*k+0.25_real64*epoch,k=0,block%pdims(3)-1)]
      call writer%write_mesh_geometry(x,y,z)
    class is (vtkhdf_structured_writer)
      if (initial) call writer%write_mesh(reshape([origin,origin+block%pdims-1],[3,2]))
      allocate(xyz(3,block%local_points(1),block%local_points(2),block%local_points(3)))
      start = block%array_origin(1)-origin(1)
      do k = 0,block%local_points(3)-1
        do j = 0,block%local_points(2)-1
          do c = 0,block%local_points(1)-1
            i = start+c
            xyz(:,c+1,j+1,k+1) = [10*id+0.5_real64*i+0.05_real64*j*k, &
                j+0.1_real64*i*k,2*k+0.25_real64*epoch]
          end do
        end do
      end do
      call writer%write_mesh_geometry(xyz,block%array_origin,block%selection)
    class is (vtkhdf_unstructured_writer)
      phase = 0
      if (block%mode == VTKHDF_DYNAMIC_MESH .and. epoch == 1) phase = 1
      call partition(id,phase,nc,start,total)
      allocate(points(3,4*nc))
      tetra = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,4])
      do c = 1,nc
        points(:,4*c-3:4*c) = tetra
        points(1,4*c-3:4*c) = points(1,4*c-3:4*c)+10*id+2*(start+c-1)
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
    class is (vtkhdf_hypertree_writer)
      if (initial) call writer%write_mesh(block%pdims)
      x = [(10*id+0.5_real64*i+0.1_real64*i*i,i=0,block%pdims(1)-1)]
      y = [(j+0.1_real64*j*j,j=0,block%pdims(2)-1)]
      z = [(2*k+0.1_real64*k*k+0.25_real64*epoch,k=0,block%pdims(3)-1)]
      call writer%write_mesh_geometry(x,y,z)
      if (initial .or. block%mode == VTKHDF_DYNAMIC_MESH) then
        call hypertree_topology(block,id,step,ids,xtree,refined,mask)
        call writer%write_mesh_topology(ids,xtree,refined,mask)
      end if
    end select
  end subroutine

  subroutine setup_data(block,id)
    type(block_info), intent(inout) :: block
    integer, intent(in) :: id
    real(real64), allocatable :: ps(:,:,:), pv(:,:,:,:), cs(:,:,:), cv(:,:,:,:)
    real(real64), allocatable :: ups(:), upv(:,:), ucs(:), ucv(:,:)
    select type (writer => block%dataset)
    class is (vtkhdf_structured_base)
      call structured_values(block,id,0,ps,pv,cs,cv)
      call writer%write_point_data('point_constant',ps,block%array_origin,block%selection)
      call writer%write_point_data('point_constant_vector',pv,block%array_origin,block%selection)
      call writer%write_cell_data('cell_constant',cs,block%array_origin)
      call writer%write_cell_data('cell_constant_vector',cv,block%array_origin)
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
    class is (vtkhdf_hypertree_writer)
      if (block%mode /= VTKHDF_DYNAMIC_MESH) then
        call hypertree_values(block,id,0,ucs,ucv)
        call writer%write_cell_data('cell_constant',ucs)
        call writer%write_cell_data('cell_constant_vector',ucv)
      end if
      call writer%write_field_data('block_id',id)
      call writer%write_field_data('direction',[1.0_real64,2.0_real64,3.0_real64],as_vector=.true.)
      if (block%mode /= VTKHDF_STATIC .and. .not.block%constant_only) then
        block%cell_scalar = writer%register_temporal_cell_data('cell_value',0.0_real64)
        block%cell_vector = writer%register_temporal_cell_data('cell_vector',[0.0_real64,0.0_real64,0.0_real64])
        block%clock = writer%register_temporal_field_data('clock',0.0_real64)
        block%table = writer%register_temporal_field_data('table',0.0_real64)
      end if
    class is (vtkhdf_unstructured_writer)
      if (block%mode /= VTKHDF_DYNAMIC_MESH) then
        call unstructured_values(block,id,0,ups,upv,ucs,ucv)
        call writer%write_point_data('point_constant',ups)
        call writer%write_point_data('point_constant_vector',upv)
        call writer%write_cell_data('cell_constant',ucs)
        call writer%write_cell_data('cell_constant_vector',ucv)
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
    integer :: i, nt
    logical :: write_points, write_cells
    real(real64), allocatable :: ps(:,:,:), pv(:,:,:,:), cs(:,:,:), cv(:,:,:,:)
    real(real64), allocatable :: ups(:), upv(:,:), ucs(:), ucv(:,:), table(:,:)
    write_points = step == 0 .or. step == 1 .or. step == 3
    write_cells = step == 0 .or. step == 2 .or. step == 4
    if (block%mode == VTKHDF_DYNAMIC_MESH .and. (step == 2 .or. step == 4)) write_points = .true.
    nt = 1+mod(step,3)
    table = reshape([(real(1000*id+100*step+i,real64),i=1,2*nt)],[2,nt])
    select type (writer => block%dataset)
    class is (vtkhdf_structured_base)
      call structured_values(block,id,step,ps,pv,cs,cv)
      if (write_points) then
        call writer%write_point_data(block%point_scalar,ps,block%array_origin,block%selection)
        call writer%write_point_data(block%point_vector,pv,block%array_origin,block%selection)
      end if
      if (write_cells) then
        call writer%write_cell_data(block%cell_scalar,cs,block%array_origin)
        call writer%write_cell_data(block%cell_vector,cv,block%array_origin)
      end if
      call writer%write_field_data(block%clock,data=times(step))
      call writer%write_field_data(block%table,data=table)
    class is (vtkhdf_hypertree_writer)
      call hypertree_values(block,id,step,ucs,ucv)
      if (write_cells) then
        call writer%write_cell_data(block%cell_scalar,ucs)
        call writer%write_cell_data(block%cell_vector,ucv)
      end if
      call writer%write_field_data(block%clock,times(step))
      call writer%write_field_data(block%table,table)
    class is (vtkhdf_unstructured_writer)
      call unstructured_values(block,id,step,ups,upv,ucs,ucv)
      if (write_points) then
        call writer%write_point_data(block%point_scalar,ups)
        call writer%write_point_data(block%point_vector,upv)
      end if
      if (write_cells) then
        call writer%write_cell_data(block%cell_scalar,ucs)
        call writer%write_cell_data(block%cell_vector,ucv)
      end if
      call writer%write_field_data(block%clock,data=times(step))
      call writer%write_field_data(block%table,data=table)
    end select
  end subroutine
end program
