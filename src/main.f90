Program generateInflow
    use, intrinsic :: ieee_arithmetic
    use mpi
    use utils
    implicit none

    ! MPI related 
    integer :: ierror, nprocs, myid
    ! Global variables
    integer, parameter :: dp = kind(1.0d0)
    integer ::  ierr, iunit
    real(dp) :: start_time
    character(len=512) :: inputfile
    logical :: filestatus
    ! Grid related (nslices == Nx)
    integer :: nslices, Ny, Nz
    integer(dp) :: i,j,k
    real(dp), allocatable, dimension(:) :: z, y
    real(dp), dimension(2) :: Lx, Ly, Lz
    ! Decomposition results
    integer :: startslice, endslice
    integer :: base_slices, remainder
    ! Synthetic Generator related
    real(dp), allocatable, dimension(:,:,:) :: random_sequence, spatially_correlated, previous_fluctuation, unscaled_fluctuation
    real(dp), allocatable, dimension(:,:,:) :: instantaneous_velocity
    real(dp) :: integralLengthScale, integralLengthScale_in, gridsizey, gridsizez
    real(dp), allocatable, dimension(:) :: b0j_y, b_y, b0j_z, b_z, dz, dy
    integer, dimension(:) :: seq_shape(3)
    integer(dp) :: n_y, n_z
    real(dp) :: pi, C_XC, deltaT, lagrangianTimeScale, alpha1, alpha2
    real(dp), allocatable :: inflowdata(:,:)
    real(dp), dimension(3,3) :: reynolds_stress_tensor, amplitude_tensor
    integer(dp) :: numcolumns, sizeofarray
    integer(dp) :: myslice, charwidth=20
    character(len=100) :: outfilename
    real(dp) :: ibulkvelocity, trueBulkVelocity
    ! Input parameters
    real(dp) :: Retau, Hin, viscIn, utau_in, bulk_input_velocity  
    ! Verbose I/O for decomposition
    logical :: verbose = .false.
!-------------------------------------------------------------------------------------------!
!                                   PROGRAM BEGINS                                          !
!-------------------------------------------------------------------------------------------!
    ! Start the MPI communication world
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)  
    ! Define pi
    pi = acos(-1.0d0)
    ! Snippet taken from CaNS - Credits to P. Costa
    iunit = 10
    open(newunit=iunit,file='parameters.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*,iostat=ierr) nslices, Ny, Nz
        read(iunit,*,iostat=ierr) Lx(2), Ly(2), Lz(2)
        read(iunit,*,iostat=ierr) integralLengthScale_in
        read(iunit,*,iostat=ierr) Retau, Hin, viscIn
        read(iunit,*,iostat=ierr) deltaT
        read(iunit,*,iostat=ierr) inputfile
      else
        error stop "parameters.in file encountered a problem!"
      end if
    close(iunit)
    ! Sync all MPI ranks before proceeding
    call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
    ! Compute the friction velocity based on input parameters
    utau_in = (Retau*viscIn)/Hin    
    ! Domain parameters  
    Lx(1) = 0.0
    Ly(1) = 0.0
    Lz(1) = 0.0
    ! Integral length scale in dimensional units
    integralLengthScale = integralLengthScale_in*(viscIn/utau_in)        ! Always 100 wall units
    ! Define the shape of the sequence
    seq_shape(1) = Ny
    seq_shape(2) = Nz
    seq_shape(3) = 3
    ! Allocate all arrays to be used
    allocate(y(Ny),z(Nz))
    allocate(random_sequence(Ny,Nz,3))
    allocate(spatially_correlated(Ny,Nz,3))
    allocate(previous_fluctuation(Ny,Nz,3))
    allocate(unscaled_fluctuation(Ny,Nz,3))
    allocate(instantaneous_velocity(Ny,Nz,3))
    ! Define the size of the array
    sizeofarray = 3
    ! Define the y and z arrays
    y(1) = Ly(1)
    do i=2,Ny
        y(i) = y(i-1) + (Ly(2)-Ly(1))/Ny
    end do
    gridsizey = y(10) - y(9)
    ! Read the file
    call read_file_skip_first(inputfile, inflowdata, numcolumns)  ! Covariance file
    ! Fix the mean velocity based on the input utau
    do k=1,Nz
        inflowdata(k,2) = inflowdata(k,2)*utau_in
    end do
    ! Define what z is
    z(:) = inflowdata(:,1)
    !Allocate the grid spacing arrays used to compute the flux
    allocate(dy(Ny),dz(Nz))
    ! Compute the grid spacing
    do k = 3,Nz
        dz(k) = z(k) - z(k-1)
    end do
    dz(2) = dz(3)
    dz(1) = dz(2)
    ! Compute the dz
    do j = 2,Ny
        dy(j) = y(j) - y(j-1)
    end do
    dy(1) = dy(2)
    ! Compute the bulk input velocity
    bulk_input_velocity = 0.0d0
    do k=1,Nz
        bulk_input_velocity = bulk_input_velocity + inflowdata(k,2)*dz(k)
    end do
    bulk_input_velocity = bulk_input_velocity/Lz(2)
    ! Compute the lagrangian time scale
    lagrangianTimeScale = integralLengthScale/(bulk_input_velocity)
    ! Fix deltaT based on user input
    if(deltaT .le. 0.0) then
        deltaT = 0.95*(Lx(2)/nslices)/bulk_input_velocity
    end if
    ! Print logo
    if(myid == 0) then
    	call printlogo()
    	print *, "*** Using ",nprocs, "MPI ranks ***"
	    print *, "-------------------------------------------------------"
    	print *, "Utau input:", utau_in
    	print *, "Integral length scale used", integralLengthScale
    	print *, "Input bulk velocity is", bulk_input_velocity
    	print *, "Lagrangian time scale: ", lagrangianTimeScale
    	print *, "Time step is ", deltaT
    	print *, "-------------------------------------------------------"
    	print *, "WARNING: The digital filter has isotropic kernel both in y and z!"
    endif
    if(myid == 0) then
        ! Check for the output file status
        inquire(file="slices", exist=filestatus)
        if (.not. filestatus) then
            ! Creating the folder for output
            print *, " - - - Output directory `slices` does not exist, creating it..."
            call system('mkdir slices')
        endif
    endif
    ! Find the value of n
    n_y = max(floor(integralLengthScale/gridsizey),1)
    n_z = max(floor(integralLengthScale/gridsizey),1)   ! assume isotropic turbulence
    ! Based on n_y allocate b0j_y and b_y
    allocate(b0j_y(-n_y:n_y),b_y(-n_y:n_y))
    allocate(b0j_z(-n_z:n_z),b_z(-n_z:n_z))
    ! Sync all MPI ranks before proceeding
    call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
    ! Calculate base number of slices per process and the remainder
    base_slices = nslices / nprocs
    remainder = mod(nslices, nprocs)
    ! Calculate the start and end slice for each process
    if (myid < remainder) then
        startslice = myid * (base_slices + 1) + 1
        endslice = startslice + base_slices
    else
        startslice = remainder * (base_slices + 1) + (myid - remainder) * base_slices + 1
        endslice = startslice + base_slices - 1
    end if
    if(verbose) print *, myid, startslice, endslice
    ! Sync all MPI ranks before proceeding
    call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
    ! Number of slices
    call cpu_time(start_time)
    do myslice=startslice,endslice
        ! Define what b model parameter is
        do i = -n_y, n_y
            b0j_y(i) = exp(-pi * abs(real(i, dp)) / (2.0*real(n_y, dp)))
        end do

        do i = -n_z, n_z
            b0j_z(i) = exp(-pi * abs(real(i, dp)) / (2.0*real(n_z, dp)))
        end do

        ! Now compute the bj model coefficient
        b_y = b0j_y/sqrt(sum(b0j_y**2))
        b_z = b0j_z/sqrt(sum(b0j_z**2))

        ! Generate the raw random sequences
        do k=1,3
            do i=1,Nz
                call generate_normal_random_numbers(0.0_dp,1.0_dp,Ny,random_sequence(:,i,k))
            end do
        end do

        ! Use the digital filter
        call digital_filter(random_sequence, b_y, b_z, spatially_correlated, seq_shape)

        ! Define the coefficients
        C_XC = pi/4.0
        alpha1 = exp(-C_XC*deltaT/lagrangianTimeScale)
        alpha2 = exp(-2.0d0*C_XC*deltaT/lagrangianTimeScale)
        if(myslice==1 .and. myid == 0) print *, "Value of alpha is:", alpha1, alpha2
        unscaled_fluctuation = alpha1 * previous_fluctuation + sqrt(1 - alpha2) * spatially_correlated
        ! Loop over z and assign the values
        do k=2,Nz
            ! Define the reynolds stress tensor
            reynolds_stress_tensor(1,1) = inflowdata(k,3)*utau_in*utau_in   ! Ruu
            reynolds_stress_tensor(1,2) = inflowdata(k,4)*utau_in*utau_in   ! Ruv
            reynolds_stress_tensor(1,3) = inflowdata(k,5)*utau_in*utau_in   ! Ruw
            reynolds_stress_tensor(2,1) = inflowdata(k,6)*utau_in*utau_in   ! Rvu
            reynolds_stress_tensor(2,2) = inflowdata(k,7)*utau_in*utau_in   ! Rvv
            reynolds_stress_tensor(2,3) = inflowdata(k,8)*utau_in*utau_in   ! Rvw
            reynolds_stress_tensor(3,1) = inflowdata(k,9)*utau_in*utau_in   ! Rwu 
            reynolds_stress_tensor(3,2) = inflowdata(k,10)*utau_in*utau_in  ! Rwv
            reynolds_stress_tensor(3,3) = inflowdata(k,11)*utau_in*utau_in  ! Rww
            ! Compute the amplitude tensor
            call cholesky(amplitude_tensor,reynolds_stress_tensor,sizeofarray)      
            do j=1,Ny
                call matrix_vector_dot_product(amplitude_tensor, unscaled_fluctuation(j,k,:), instantaneous_velocity(j,k,:))
                instantaneous_velocity(j,k,1) = instantaneous_velocity(j,k,1) + inflowdata(k,2)
            end do
        end do
        ! Correct the mass flux to make it divergence free U dot n_x / A_in
        ! NOTE: Mass flux is not corrected!
        ! Here U is the velocity vector, n_x is the x unit vector and A_in is the inlet area
        ibulkvelocity = 0.0d0
        do j=1,Ny
            do k=1,Nz
                ibulkvelocity = ibulkvelocity + instantaneous_velocity(j,k,1)*dz(k)*dy(j)
            end do 
        end do
        ibulkvelocity = ibulkvelocity/(sum(dz)*sum(dy))
        trueBulkVelocity = 0.0d0
        do k=1,Nz
            trueBulkVelocity = trueBulkVelocity + inflowdata(k,2)*dz(k)
        end do
        trueBulkVelocity = trueBulkVelocity/Lz(2)
        instantaneous_velocity(:,:,1) = instantaneous_velocity(:,:,1)*(trueBulkVelocity/ibulkvelocity)
        ! Set the previous fluctuation as the current one for time correlations
        previous_fluctuation = unscaled_fluctuation
        ! Write individual slices
        write(outfilename, '("slices/uslicedata_", I0, ".dat")') myslice
        call write_2d_array_to_file(outfilename,instantaneous_velocity(:,:,1))
        write(outfilename, '("slices/vslicedata_", I0, ".dat")') myslice
        call write_2d_array_to_file(outfilename,instantaneous_velocity(:,:,2))
        write(outfilename, '("slices/wslicedata_", I0, ".dat")') myslice
        call write_2d_array_to_file(outfilename,instantaneous_velocity(:,:,3))
        if(myid == 1) call show_progress(myslice-startslice+1,int(endslice-startslice+1,dp),charwidth,start_time)
    ! End my slice loop    
    end do
    ! Sync all MPI ranks before proceeding
    call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
    ! Finalise the MPI communication
    call MPI_FINALIZE(ierror)
End Program generateInflow
