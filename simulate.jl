function simulate(	D::Float64,
					k_on::Float64,
					k_off::Float64,
					mobile_fraction::Float64,
					alpha::Float64,
					beta::Float64,
					bleach_region_shape::Float64,
					r_bleach::Float64,
					lx_bleach::Float64,
					ly_bleach::Float64,
					number_of_pixels::Int64,
					number_of_pad_pixels::Int64,
					number_of_prebleach_frames::Int64,
					number_of_bleach_frames::Int64,
					number_of_postbleach_frames::Int64,
					delta_t::Float64,
					number_of_time_steps_fine_per_course::Int64,
					number_of_particles::Int64)

	# Simulation parameters.
	number_of_pixels_float::Float64 = convert(Float64, number_of_pixels)
	number_of_pad_pixels_float::Float64 = convert(Float64, number_of_pad_pixels)

	# Compute parameters of Markov chain from the on and off reaction rates.
	tau_u::Float64 = 1.0 / k_on
	tau_b::Float64 = 1.0 / k_off
	p_u::Float64 = k_off / ( k_on + k_off )
	p_b::Float64 = k_on / ( k_on + k_off )

	# If simulating pure diffusion:
	if k_on == 0.0
		number_of_time_steps_fine_per_course = 1
	end

	delta_t_fine::Float64 = delta_t / convert(Float64, number_of_time_steps_fine_per_course)
	sigma_fine::Float64 = 0.0

	p_bu::Float64 = delta_t_fine / tau_b
	p_ub::Float64 = delta_t_fine / tau_u

	# Initialization of variables.
	x::Float64 = 0.0
	y::Float64 = 0.0
	is_inside_bleach_region::Bool = false
	is_mobile::Bool = false
	is_unbound::Bool = false
	is_bleached::Bool = false
	ind_x::Int64 = 0
	ind_y::Int64 = 0

	C_prebleach::Array{Int64, 3} = zeros(number_of_pixels, number_of_pixels, number_of_prebleach_frames)
	C_postbleach::Array{Int64, 3} = zeros(number_of_pixels, number_of_pixels, number_of_postbleach_frames)

	# Simulate.
	for current_particle = 1:number_of_particles
		if mod(current_particle, 100_000) == 0
			println(current_particle)
		end

		# Random initial position.
		x = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()
		y = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()

		is_bleached = false

		# Is mobile?
		is_mobile = rand() <= mobile_fraction
		if is_mobile
			sigma_fine = sqrt( 2.0 * D * delta_t_fine )
		else
			sigma_fine = 0.0
		end

		# Initial state.
		if rand() <= p_u
			is_unbound = true
		else
			is_unbound = false
		end

		# Prebleach.
		for current_frame = 1:number_of_prebleach_frames
			# Propagate.
			for current_time_step_fine = 1:number_of_time_steps_fine_per_course
				if is_unbound
					x = x + sigma_fine * randn()
					y = y + sigma_fine * randn()
				end

				if is_unbound
					if rand() <= p_ub
						is_unbound = false
					end
				else
					if rand() <= p_bu
						is_unbound = true
					end
				end
			end
			x = mod(x, number_of_pixels_float + 2.0 * number_of_pad_pixels_float)
			y = mod(y, number_of_pixels_float + 2.0 * number_of_pad_pixels_float)

			# Imaging bleach.
			is_inside_bleach_region = false
			if (number_of_pad_pixels_float <= x <= number_of_pad_pixels_float + number_of_pixels_float) & (number_of_pad_pixels_float <= y <= number_of_pad_pixels_float + number_of_pixels_float)
				is_inside_bleach_region = true
			end

			if is_inside_bleach_region & (rand() <= (1.0 - beta))
				is_bleached = true
			end


			# Form image.
			if !is_bleached
				ind_x = convert(Int64, ceil(x - number_of_pad_pixels_float))
				if ind_x >= 1 && ind_x <= number_of_pixels
					ind_y = convert(Int64, ceil(y - number_of_pad_pixels_float))
					if ind_y >= 1 && ind_y <= number_of_pixels
						C_prebleach[ind_x, ind_y, current_frame] += 1
					end
				end
			end
		end

		# Bleach.
		for current_frame = 1:number_of_bleach_frames
			# Propagate.
			for current_time_step_fine = 1:number_of_time_steps_fine_per_course
				if is_unbound
					x = x + sigma_fine * randn()
					y = y + sigma_fine * randn()
				end

				if is_unbound
					if rand() <= p_ub
						is_unbound = false
					end
				else
					if rand() <= p_bu
						is_unbound = true
					end
				end
			end
			x = mod(x, number_of_pixels_float + 2.0 * number_of_pad_pixels_float)
			y = mod(y, number_of_pixels_float + 2.0 * number_of_pad_pixels_float)

			# Bleach.
			is_inside_bleach_region = false
			if bleach_region_shape == 0.0 # Circle
				if ( x - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 + ( y - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 <= r_bleach^2
					is_inside_bleach_region = true
				end
			else # Rectangle
				if (number_of_pad_pixels_float + 0.5 * number_of_pixels_float - 0.5 * lx_bleach <= x <= number_of_pad_pixels_float + 0.5 * number_of_pixels_float + 0.5 * lx_bleach) &
				   (number_of_pad_pixels_float + 0.5 * number_of_pixels_float - 0.5 * ly_bleach <= y <= number_of_pad_pixels_float + 0.5 * number_of_pixels_float + 0.5 * ly_bleach)
					is_inside_bleach_region = true
				end
			end

			if is_inside_bleach_region & (rand() <= (1.0 - alpha))
				is_bleached = true
			end

			# Imaging bleach.
			is_inside_bleach_region = false
			if (number_of_pad_pixels_float <= x <= number_of_pad_pixels_float + number_of_pixels_float) & (number_of_pad_pixels_float <= y <= number_of_pad_pixels_float + number_of_pixels_float)
				is_inside_bleach_region = true
			end

			if is_inside_bleach_region & (rand() <= (1.0 - beta))
				is_bleached = true
			end
		end

		# Postbleach.
		for current_frame = 1:number_of_postbleach_frames
			# Propagate.
			for current_time_step_fine = 1:number_of_time_steps_fine_per_course
				if is_unbound
					x = x + sigma_fine * randn()
					y = y + sigma_fine * randn()
				end

				if is_unbound
					if rand() <= p_ub
						is_unbound = false
					end
				else
					if rand() <= p_bu
						is_unbound = true
					end
				end
			end
			x = mod(x, number_of_pixels_float + 2.0 * number_of_pad_pixels_float)
			y = mod(y, number_of_pixels_float + 2.0 * number_of_pad_pixels_float)

			# Imaging bleach.
			is_inside_bleach_region = false
			if (number_of_pad_pixels_float <= x <= number_of_pad_pixels_float + number_of_pixels_float) & (number_of_pad_pixels_float <= y <= number_of_pad_pixels_float + number_of_pixels_float)
				is_inside_bleach_region = true
			end

			if is_inside_bleach_region & (rand() <= (1.0 - beta))
				is_bleached = true
			end

			# Form image.
			if !is_bleached
				ind_x = convert(Int64, ceil(x - number_of_pad_pixels_float))
				if ind_x >= 1 && ind_x <= number_of_pixels
					ind_y = convert(Int64, ceil(y - number_of_pad_pixels_float))
					if ind_y >= 1 && ind_y <= number_of_pixels
						C_postbleach[ind_x, ind_y, current_frame] += 1
					end
				end
			end
		end
	end

	return cat(dims = 3, C_prebleach, C_postbleach)
end
