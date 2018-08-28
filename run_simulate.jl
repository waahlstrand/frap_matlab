using Random

@everywhere include("simulate.jl")

function run_simulate(	D::Float64,
						k_on::Float64,
						k_off::Float64,
						mobile_fraction::Float64,
						alpha::Float64,
						beta::Float64)

	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	Random.seed!(random_seed)

	# Store start time.
	t_start::Int64 = convert(Int64, time_ns())

	# Experimental and simulation parameters.
	pixel_size::Float64 = 7.5e-7 # m
	number_of_pixels::Int64 = 256 #  pixels

	number_of_prebleach_frames::Int64 = 10
	number_of_bleach_frames::Int64 = 2
	number_of_postbleach_frames::Int64 = 50
	delta_t::Float64 = 0.2 # s

	number_of_pad_pixels::Int64 = 128 # pixels
	number_of_time_steps_fine_per_course::Int64 = 32
	number_of_particles_per_worker::Int64 = 1000000
	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.

	r_bleach::Float64 = 15e-6 / pixel_size

	# Simulate data.
	data::Array{Int64, 3} = @distributed (+) for current_worker = 1:number_of_workers
		simulate(	D,
					k_on,
					k_off,
					mobile_fraction,
					alpha,
					beta,
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

	# Measure and print execution time.
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	t_exec_s::Float64 = convert(Float64, t_exec) / 1e9
	println(join(("Execution time: ", t_exec_s, " seconds.")))

	# Save output.
	file_name_output::String = join(("simulated_stochastic_data_", string(D), "_", string(k_on), "_", string(k_off), "_", string(mobile_fraction), "_", string(alpha), "_", string(beta), ".bin"))
	file_stream_output::IOStream = open(file_name_output, "w")
	write(file_stream_output, D)
	write(file_stream_output, k_on)
	write(file_stream_output, k_off)
	write(file_stream_output, mobile_fraction)
	write(file_stream_output, alpha)
	write(file_stream_output, beta)
	write(file_stream_output, r_bleach)
	write(file_stream_output, number_of_pixels)
	write(file_stream_output, number_of_prebleach_frames)
	write(file_stream_output, number_of_bleach_frames)
	write(file_stream_output, number_of_postbleach_frames)
	write(file_stream_output, number_of_pad_pixels)
	write(file_stream_output, delta_t)
	write(file_stream_output, pixel_size)
	write(file_stream_output, number_of_particles_per_worker * number_of_workers)
	write(file_stream_output, t_exec_s)
	write(file_stream_output, data)

	close(file_stream_output)

	nothing
end

#alpha = 0.6
#for D in (10.0, 100.0, 1000.0)
#	for k_on in (0.0, 1.0, 10.0)
#		for k_off in (1.0, 10.0)
#			for mobile_fraction in (0.5, 1.0)
#				for beta in (1.0, 0.999, 0.99)
#					println((D, k_on, k_off, mobile_fraction, alpha, beta))
#					run_simulate(D, k_on, k_off, mobile_fraction, alpha, beta)
#				end
#			end
#		end
#	end
#end
D = 1000.0
k_on = 1.0
k_off = 10.0
mobile_fraction = 1.0
alpha = 0.6
beta = 0.99
run_simulate(D, k_on, k_off, mobile_fraction, alpha, beta)
