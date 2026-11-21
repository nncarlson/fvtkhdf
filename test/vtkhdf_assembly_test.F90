! Shared MB/PDC construction, commitment, and lifetime fixtures.
program vtkhdf_assembly_test
  use,intrinsic :: iso_fortran_env, only: int8, real64
#ifdef PDC_COLLECTION
  use vtkhdf_pdc_file_type, composite_file => vtkhdf_pdc_file
#else
  use vtkhdf_mb_file_type, composite_file => vtkhdf_mb_file
#endif
  use vtkhdf_h5_c_binding, only: hid_t, H5Iget_file_id, H5Iis_valid, H5Lexists, H5Fclose
  use vtkhdf_test_env
#ifdef USE_MPI_F08
  use mpi_f08, only: MPI_COMM_WORLD
#elif defined(USE_MPI)
  use mpi, only: MPI_COMM_WORLD
#endif
  implicit none
#ifdef PDC_COLLECTION
  character(*), parameter :: prefix = 'pdc'
#else
  character(*), parameter :: prefix = 'mb'
#endif

  type(composite_file) :: file
  type(vtkhdf_unstructured_block), pointer :: fluid, other, tail
  type(vtkhdf_image_block), pointer :: solid
  type(vtkhdf_cell_data_handle) :: pressure, temperature
  integer :: rank, nproc, stat, case, regions, part, materials, node, j, ids(96), inner
  integer :: bounds(3,2), selected(3,2)
  integer(hid_t) :: root, file_id
  real :: values(1,1,1)
  character(:), allocatable :: errmsg
  character(24) :: label
  character(12), parameter :: cases(*) = [character(12) :: &
      'explicit', 'implicit', 'flat', 'empty', 'rootrefs', 'early', 'late', 'growth', 'rootgrowth']

  call test_env_init(rank,nproc)
  do case = 1, size(cases)
#ifdef USE_MPI
    call file%create(prefix//'_assembly_'//trim(cases(case))//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
    call file%create(prefix//'_assembly_'//trim(cases(case))//'.vtkhdf',stat,errmsg)
#endif
    if (stat /= 0) error stop errmsg
    select case (trim(cases(case)))
    case ('explicit','implicit')
      call file%add_block('Fluid.Name',fluid,mode=VTKHDF_FIXED_MESH)
      call write_triangle(fluid,1)
      pressure = fluid%register_temporal_cell_data('pressure',0.0)
      root = fluid%root()
      regions = file%add_node('Regions')
      node = file%add_node('Fluid Name',parent=regions)
      call file%add_block_ref(fluid,parent=regions)
      call file%add_block('Fluid Name',solid,mode=VTKHDF_FIXED_MESH)
      bounds = reshape([0,0,0,nproc,1,0],[3,2])
      selected = reshape([rank,0,0,rank,0,0],[3,2])
      call solid%write_mesh(bounds,[0.,2.,0.],[1.,1.,1.])
      call solid%write_field_data('id',2)
      temperature = solid%register_temporal_cell_data('temperature',0.0)
      call file%add_block_ref(solid,parent=regions)
      part = file%add_node('Part',parent=regions)
      node = file%add_node('Common',parent=part)
      call file%add_block_ref(fluid,parent=part)
      materials = file%add_node('Regions')
      call file%add_block_ref(fluid,parent=materials)
      node = file%add_node('Common',parent=materials)
      call file%add_block_ref(fluid,parent=0)
      node = file%add_node('Fluid Name')
      call file%add_block('Other',other,mode=VTKHDF_FIXED_MESH)
      call write_triangle(other,3)
      call file%add_block('Tail',tail,mode=VTKHDF_FIXED_MESH)
      call write_triangle(tail,4)
      ! Custom construction may continue after the timeline starts.
      call file%start_time_step(0.0_real64)
      node = file%add_node('Other')
      node = file%add_node('Other')
      call file%add_block_ref(fluid,parent=regions) ! duplicate retains its position
      call file%flush()
      file_id = H5Iget_file_id(root)
      if (H5Lexists(file_id,'/VTKHDF/Assembly/Regions') /= 0) error stop 'premature hierarchy'
      if (H5Lexists(file_id,'/VTKHDF/Assembly/Fluid_Name') <= 0) error stop 'flat link lost'
      if (H5Lexists(file_id,'/VTKHDF/Assembly/Fluid_Name_1') <= 0) error stop 'flat link lost'
      if (H5Lexists(file_id,'/VTKHDF/Assembly/Other') <= 0) error stop 'implicit link lost'
      if (trim(cases(case)) == 'explicit') then
        call file%set_assembly()
        if (H5Lexists(file_id,'/VTKHDF/Assembly/Regions/Fluid_Name') <= 0) error stop 'missing hierarchy'
      end if
      if (fluid%root() /= root) error stop 'writer target changed'
      ! These writes happen after explicit commitment, including in an open step.
      call fluid%write_field_data('after',42)
      call fluid%write_cell_data(pressure,[100.0+rank])
      values = 300.0+rank
      call solid%write_cell_data(temperature,values,[rank,0,0],selected)
      call file%start_time_step(1.0_real64)
      call fluid%write_cell_data(pressure,[101.0+rank])
      values = 301.0+rank
      call solid%write_cell_data(temperature,values,[rank,0,0],selected)
      stat = H5Fclose(file_id)
      if (stat < 0) error stop 'extra file reference'
      ! CLOSE commits the implicit case before finalizing the open block steps.
    case ('flat')
      call file%add_block('A',fluid)
      call write_triangle(fluid,1)
      call file%add_block('B',other)
      call write_triangle(other,2)
      call file%set_assembly()
      call fluid%write_field_data('after',42)
    case ('empty')
      regions = file%add_node('Empty')
      node = file%add_node('',parent=regions)
      node = file%add_node('Empty')
      ! No blocks, and no explicit commitment.
    case ('rootrefs')
      call file%add_block('A',fluid)
      call write_triangle(fluid,1)
      call file%add_block('B',other)
      call write_triangle(other,2)
      call file%add_block_ref(other)
      call file%add_block_ref(fluid)
      call file%add_block_ref(other) ! no ordinary nodes; still requires commitment
    case ('early','late')
      if (trim(cases(case)) == 'early') then
        regions = file%add_node('A')
        node = file%add_node('B')
      end if
      call file%add_block('A',fluid)
      call write_triangle(fluid,1)
      call file%add_block('B',other)
      call write_triangle(other,2)
      if (trim(cases(case)) == 'late') then
        regions = file%add_node('A')
        node = file%add_node('B')
      end if
      call file%add_block_ref(fluid,parent=regions)
      call file%set_assembly()
    case ('growth')
      call file%add_block('A',fluid)
      call write_triangle(fluid,1)
      do j = 1, size(ids)
        write(label,'(a,i3.3)') 'Node_',j
        ids(j) = file%add_node(trim(label))
        inner = file%add_node('Inside',parent=ids(j))
        if (j == 1) then
          regions = inner
          call file%add_block_ref(fluid,parent=regions)
        end if
      end do
      call file%add_block('B',other)
      call write_triangle(other,2)
      call file%add_block_ref(other,parent=inner)
      call file%add_block_ref(fluid,parent=regions) ! survives multiple container growths
      call file%add_block_ref(fluid)
      call file%set_assembly()
    case ('rootgrowth')
      call file%add_block('A',fluid)
      call write_triangle(fluid,1)
      do j = 2, 5
        label = achar(iachar('A')+j-1)
        node = file%add_node(trim(label))
        if (j == 2) regions = node
        call file%add_block(trim(label),other)
        call write_triangle(other,j)
      end do
      call file%add_block_ref(fluid,parent=regions)
      ! Implicit root references must grow the full root child array before
      ! resolving the colliding node names, while preserving the explicit order.
      call file%set_assembly()
    end select
    call file%close()
    call file%close() ! resets all construction state and remains idempotent
    if (trim(cases(case)) == 'explicit' .or. trim(cases(case)) == 'implicit') then
      if (H5Iis_valid(root) > 0) error stop 'block root leaked'
      if (H5Iis_valid(file_id) > 0) error stop 'file leaked'
    end if
    nullify(fluid,solid,other,tail)
  end do
  call test_env_finalize()

contains

  subroutine write_triangle(block, id)
    type(vtkhdf_unstructured_block), intent(inout) :: block
    integer, intent(in) :: id
    real :: points(3,3)
    points = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.],[3,3])
    points(1,:) = points(1,:) + rank
    call block%write_mesh(points,[1,2,3],[1,4],[5_int8])
    call block%write_field_data('id',id)
    call block%write_cell_data('constant',[real(id+rank)])
  end subroutine

end program
