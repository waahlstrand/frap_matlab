using Random

@everywhere include("simulate.jl")

function run_simulate(	D_SI::Float64,
						k_on::Float64,
						k_off::Float64,
						mobile_fraction::Float64,
						alpha::Float64,
						beta::Float64,
						gamma::Float64,
						bleach_region_shape::Int64,
						r_bleach::Float64,
						lx_bleach::Float64,
						ly_bleach::Float64,
						pixel_size::Float64,
						number_of_prebleach_frames::Int64,
						number_of_bleach_frames::Int64,
						number_of_postbleach_frames::Int64,
						file_name_output::String)

	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	Random.seed!(random_seed)

	# Store start time.
	t_start::Int64 = convert(Int64, time_ns())

	# Experimental and simulation parameters.

	number_of_pixels::Int64 = 256 #  pixels

	delta_t::Float64 = 0.2 # s

	number_of_pad_pixels::Int64 = 128 # pixels
	number_of_time_steps_fine_per_course::Int64 = 32
	number_of_particles::Int64 = 20_000_000_000#40_000_000_000

	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.
	number_of_particles_per_worker::Array{Int64, 1} = convert(Array{Int64, 1}, floor(number_of_particles / number_of_workers) * ones(number_of_workers))
	number_of_particles_remaining::Int64 = number_of_particles - sum(number_of_particles_per_worker)
	number_of_particles_per_worker[1:number_of_particles_remaining] .+= 1

	D::Float64 = D_SI / pixel_size^2

	# Simulate data.
	data::Array{Int64, 3} = @distributed (+) for current_worker = 1:number_of_workers
		simulate(	D,
					k_on,
					k_off,
					mobile_fraction,
					alpha,
					beta,
					gamma,
					bleach_region_shape,
					r_bleach,
					lx_bleach,
					ly_bleach,
					number_of_pixels,
					number_of_pad_pixels,
					number_of_prebleach_frames,
					number_of_bleach_frames,
					number_of_postbleach_frames,
					delta_t,
					number_of_time_steps_fine_per_course,
					number_of_particles_per_worker[current_worker])
	end

	# Measure and print execution time.
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	t_exec_s::Float64 = convert(Float64, t_exec) / 1e9
	println(join(("Execution time: ", t_exec_s, " seconds.")))

	# Save output.
	file_stream_output::IOStream = open(file_name_output, "w")
	write(file_stream_output, D_SI)
	write(file_stream_output, k_on)
	write(file_stream_output, k_off)
	write(file_stream_output, mobile_fraction)
	write(file_stream_output, alpha)
	write(file_stream_output, beta)
	write(file_stream_output, gamma)
	write(file_stream_output, bleach_region_shape)
	write(file_stream_output, r_bleach)
	write(file_stream_output, lx_bleach)
	write(file_stream_output, ly_bleach)
	write(file_stream_output, number_of_pixels)
	write(file_stream_output, number_of_prebleach_frames)
	write(file_stream_output, number_of_bleach_frames)
	write(file_stream_output, number_of_postbleach_frames)
	write(file_stream_output, number_of_pad_pixels)
	write(file_stream_output, delta_t)
	write(file_stream_output, pixel_size)
	write(file_stream_output, number_of_particles)
	write(file_stream_output, t_exec_s)
	write(file_stream_output, data)

	close(file_stream_output)

	nothing
end

# Run simulations.
pixel_size = 7.5e-7 # m
number_of_prebleach_frames = 1
number_of_postbleach_frames = 1

# Small circle, no bleach = 4, D
D_SI = 1e-11
k_on = 0.0
k_off = 1.0
mobile_fraction = 1.0
alpha = 0.9
beta = 0.998
gamma = 0.0
bleach_region_shape = 0
r_bleach = 15e-6 / pixel_size
lx_bleach = 20e-6 / pixel_size
ly_bleach = 20e-6 / pixel_size
number_of_bleach_frames = 4
file_name_output = "sim_1.bin"
run_simulate(D_SI, k_on, k_off, mobile_fraction, alpha, beta, gamma, bleach_region_shape, r_bleach, lx_bleach, ly_bleach, pixel_size, number_of_prebleach_frames, number_of_bleach_frames, number_of_postbleach_frames, file_name_output)

# Large circle, no bleach = 2, DB
D_SI = 5e-10
k_on = 1.0
k_off = 0.5
mobile_fraction = 1.0
alpha = 0.8
beta = 1.0
gamma = 0.0
bleach_region_shape = 0
r_bleach = 25e-6 / pixel_size
lx_bleach = 20e-6 / pixel_size
ly_bleach = 20e-6 / pixel_size
number_of_bleach_frames = 4
file_name_output = "sim_2.bin"
run_simulate(D_SI, k_on, k_off, mobile_fraction, alpha, beta, gamma, bleach_region_shape, r_bleach, lx_bleach, ly_bleach, pixel_size, number_of_prebleach_frames, number_of_bleach_frames, number_of_postbleach_frames, file_name_output)

# Large circle, no bleach = 2, DB
D_SI = 1e-10
k_on = 0.0
k_off = 1.0
mobile_fraction = 0.8
alpha = 0.7
beta = 1.0
gamma = 2.0
bleach_region_shape = 1
r_bleach = 25e-6 / pixel_size
lx_bleach = 20e-6 / pixel_size
ly_bleach = 20e-6 / pixel_size
number_of_bleach_frames = 1
file_name_output = "sim_3.bin"
run_simulate(D_SI, k_on, k_off, mobile_fraction, alpha, beta, gamma, bleach_region_shape, r_bleach, lx_bleach, ly_bleach, pixel_size, number_of_prebleach_frames, number_of_bleach_frames, number_of_postbleach_frames, file_name_output)
