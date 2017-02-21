workspace()

@everywhere include("simulate.jl")

function simulate_diffusion_and_binding()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Store start time.
	t_start::Int64 = convert(Int64, time_ns())
	
	# Experimental parameters. 
	delta_t::Float64 = 0.2650 # s
	number_of_post_bleach_images::Int64 = 10#4#40
	number_of_pixels::Int64 = 256 #  pixels
	r_bleach::Float64 = 32.0 #  pixels
	intensity_inside_bleach_region::Float64 = 0.6
	intensity_outside_bleach_region::Float64 = 0.9
	pixel_size::Float64 = 7.598e-7 # m

	# Particle parameters.
	D_SI::Float64 = 5.0e-10 # m^2 / s
	D::Float64 = D_SI / pixel_size^2 # pixels^2 / s
	k_on::Float64 = 0.5 # 1/s
	k_off::Float64 = 1.0#1.0 # 1/s
	mobile_fraction::Float64 = 0.9

	# Simulation parameters.
	number_of_pad_pixels::Int64 = 256 # pixels
	
	number_of_time_steps_fine_per_course::Int64 = 64
	
	number_of_particles_per_worker::Int64 = 10000000
	number_of_workers::Int64 = nworkers()

	# Generate particle trajectories and FRAP image data.
	image_data_post_bleach::Array{Int64, 3} = @parallel (+) for current_worker = 1:number_of_workers
		simulate(	D, 
					k_on, 
					k_off, 
					mobile_fraction,
					r_bleach,
					number_of_pixels,
					number_of_post_bleach_images,
					number_of_pad_pixels,
					delta_t,
					number_of_time_steps_fine_per_course,
					intensity_inside_bleach_region,
					intensity_outside_bleach_region,
					number_of_particles_per_worker)
	end

	# Save output.
	file_name_output::String = "simulated_frap_data.dat"
	file_stream_output::IOStream = open(file_name_output, "w")
	write(file_stream_output, D)
	write(file_stream_output, k_on)
	write(file_stream_output, k_off)
	write(file_stream_output, mobile_fraction)
	write(file_stream_output, delta_t)
	write(file_stream_output, r_bleach)
	write(file_stream_output, pixel_size)
	write(file_stream_output, number_of_pixels)
	write(file_stream_output, number_of_post_bleach_images)
	write(file_stream_output, number_of_pad_pixels)
	write(file_stream_output, intensity_inside_bleach_region)
	write(file_stream_output, intensity_outside_bleach_region)
	write(file_stream_output, number_of_particles_per_worker*number_of_workers)
	write(file_stream_output, image_data_post_bleach)
	close(file_stream_output)
	
	# Measure and print execution time.
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(join(("Execution time: ", round(t_exec/1e9 * 10.0) / 10.0, " seconds.")))

	nothing
end

simulate_diffusion_and_binding()
