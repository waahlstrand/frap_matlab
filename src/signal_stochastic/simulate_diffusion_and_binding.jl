workspace()

@everywhere include("simulate.jl")

function simulate_diffusion_and_binding(D_SI::Float64, k_on::Float64, k_off::Float64, mf::Float64)
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Store start time.
	t_start::Int64 = convert(Int64, time_ns())
	
	# Experimental parameters. 
	pixel_size::Float64 = 7.5e-7 # m
	delta_t::Float64 = 0.2 # s
	number_of_images::Int64 = 50
	number_of_pixels::Int64 = 256 #  pixels
	r_bleach::Float64 = 15e-6 / pixel_size # pixels corresponding to 15 µm radius (30 µm diameter)
	Ib::Float64 = 0.6
	Iu::Float64 = 1.0

	# Particle parameters.
	D::Float64 = D_SI / pixel_size^2 # pixels^2 / s

	# Simulation parameters.
	number_of_pad_pixels::Int64 = 128 # pixels
	number_of_time_steps_fine_per_course::Int64 = 32
	number_of_particles_per_worker::Int64 = 200000000
	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.

	# Generate particle trajectories and FRAP image data.
	data::Array{Int64, 3} = @parallel (+) for current_worker = 1:number_of_workers
		simulate(	D, 
					k_on, 
					k_off, 
					mf,
					Ib,
					Iu,
					r_bleach,
					number_of_pixels,
					number_of_images,
					number_of_pad_pixels,
					delta_t,
					number_of_time_steps_fine_per_course,
					number_of_particles_per_worker)
	end

	# Save output.
	file_name_output::String = join(("simulated_stochastic_data_", string(D_SI), "_", string(k_on), "_", string(k_off), ".dat"))
	file_stream_output::IOStream = open(file_name_output, "w")
	write(file_stream_output, D)
	write(file_stream_output, k_on)
	write(file_stream_output, k_off)
	write(file_stream_output, mf)
	write(file_stream_output, Ib)
	write(file_stream_output, Iu)
	write(file_stream_output, r_bleach)
	write(file_stream_output, number_of_pixels)
	write(file_stream_output, number_of_images)
	write(file_stream_output, number_of_pad_pixels)
	write(file_stream_output, delta_t)
	write(file_stream_output, pixel_size)
	write(file_stream_output, number_of_particles_per_worker * number_of_workers)
	write(file_stream_output, data)
	close(file_stream_output)
	
	# Measure and print execution time.
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(join(("Execution time: ", round(t_exec/1e9 * 10.0) / 10.0, " seconds.")))

	nothing
end

D_SI_VECTOR = [5e-12, 1e-11, 5e-11, 1e-10, 5e-10]
K_ON_VECTOR = [0.05, 0.1, 0.5, 1.0, 5.0]
K_OFF_VECTOR = [0.05, 0.1, 0.5, 1.0, 5.0]
for ind_D = 1:length(D_SI_VECTOR)
	for ind_k_on = 1:length(K_ON_VECTOR)
		for ind_k_off = 1:length(K_OFF_VECTOR)
			simulate_diffusion_and_binding(D_SI_VECTOR[ind_D], K_ON_VECTOR[ind_k_on], K_OFF_VECTOR[ind_k_off], 1.0)
		end
	end
end
