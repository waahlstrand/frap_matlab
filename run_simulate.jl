using Random

@everywhere include("simulate.jl")

function run_simulate()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	Random.seed!(random_seed)

	# Store start time.
	t_start::Int64 = convert(Int64, time_ns())

	# Experimental and simulation parameters.
	pixel_size::Float64 = 7.5e-7 # m
	number_of_pixels::Int64 = 256 #  pixels

	number_of_prebleach_frames::Int64 = 20
	number_of_bleach_frames::Int64 = 4
	number_of_postbleach_frames::Int64 = 50
	delta_t::Float64 = 0.2 # s

	number_of_pad_pixels::Int64 = 128 # pixels
	number_of_time_steps_fine_per_course::Int64 = 32
	number_of_particles_per_worker::Int64 = 100000#00
	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.

	r_bleach::Float64 = 15e-6 / pixel_size

	# System parameters.
	D_SI::Float64 = 5e-11
	D::Float64 = D_SI / pixel_size^2 # pixels^2 / s
	k_on::Float64 = 0.0
	k_off::Float64 = 1.0
	mobile_fraction::Float64 = 1.0
	alpha::Float64 = 0.6

	# Simulate data.
	data::Array{Int64, 3} = @distributed (+) for current_worker = 1:number_of_workers
		simulate(	D,
					k_on,
					k_off,
					mobile_fraction,
					alpha,
					r_bleach,
					number_of_pixels,
					number_of_pad_pixels,
					number_of_prebleach_frames,
					number_of_bleach_frames,
					number_of_postbleach_frames,
					delta_t,
					number_of_time_steps_fine_per_course,
					number_of_particles_per_worker)
	end
	C_prebleach::Array{Int64, 3} = data[:, :, 1:number_of_prebleach_frames]
	C_postbleach::Array{Int64, 3} = data[:, :, number_of_prebleach_frames+1:end]

	# Measure and print execution time.
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(join(("Execution time: ", round(t_exec/1e9 * 10.0) / 10.0, " seconds.")))

	# Save output.
	file_name_output::String = join(("simulated_stochastic_data_", string(D_SI), "_", string(k_on), "_", string(k_off), ".bin"))
	file_stream_output::IOStream = open(file_name_output, "w")
	write(file_stream_output, D)
	write(file_stream_output, k_on)
	write(file_stream_output, k_off)
	write(file_stream_output, mobile_fraction)
	write(file_stream_output, alpha)
	write(file_stream_output, r_bleach)
	write(file_stream_output, number_of_pixels)
	write(file_stream_output, number_of_prebleach_frames)
	write(file_stream_output, number_of_bleach_frames)
	write(file_stream_output, number_of_postbleach_frames)
	write(file_stream_output, number_of_pad_pixels)
	write(file_stream_output, delta_t)
	write(file_stream_output, pixel_size)
	write(file_stream_output, number_of_particles_per_worker * number_of_workers)
	write(file_stream_output, C_prebleach)
	write(file_stream_output, C_postbleach)
	write(file_stream_output, round(t_exec/1e9 * 10.0) / 10.0)
	close(file_stream_output)

	nothing
end

run_simulate()
