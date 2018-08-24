function simulate(	D::Float64,
					alpha::Float64,
					r_bleach::Float64,
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

	delta_t_fine::Float64 = delta_t / convert(Float64, number_of_time_steps_fine_per_course)
	sigma_fine::Float64 = sqrt(2.0 * D * delta_t_fine )

	# Initialization of variables.
	x::Float64 = 0.0
	y::Float64 = 0.0
	is_inside_bleach_region::Bool = false
	is_bleached::Bool = false
	ind_x::Int64 = 0
	ind_y::Int64 = 0

	C_prebleach::Array{Int64, 3} = zeros(number_of_pixels, number_of_pixels, number_of_prebleach_frames)
	C_postbleach::Array{Int64, 3} = zeros(number_of_pixels, number_of_pixels, number_of_postbleach_frames)

	# Simulate.
	for current_particle = 1:number_of_particles
		if mod(current_particle, 1000000) == 0
			println(current_particle)
		end

		# Random initial position.
		x = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()
		y = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()
		# Prebleach.
		for current_frame = 1:number_of_prebleach_frames
			# Propagate.
			for current_time_step_fine = 1:number_of_time_steps_fine_per_course
					x = x + sigma_fine * randn()
					y = y + sigma_fine * randn()

			end

			# Form image.
			ind_x = convert(Int64, ceil(mod(x, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
			if ind_x >= 1 && ind_x <= number_of_pixels
				ind_y = convert(Int64, ceil(mod(y, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
				if ind_y >= 1 && ind_y <= number_of_pixels
					C_prebleach[ind_x, ind_y, current_frame] += 1
				end
			end
		end

		# Bleach.
		is_bleached = false
		for current_frame = 1:number_of_bleach_frames
			# Propagate.
			for current_time_step_fine = 1:number_of_time_steps_fine_per_course
					x = x + sigma_fine * randn()
					y = y + sigma_fine * randn()
			end

			# Bleach (with certain probability).
			is_inside_bleach_region = false
			if ( x - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 + ( y - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 <= r_bleach^2
				is_inside_bleach_region = true
			end

			if is_inside_bleach_region & !is_bleached
				is_bleached = rand() <= (1.0 - alpha)
			end
		end

		# Postbleach.
		if !is_bleached
			for current_frame = 1:number_of_postbleach_frames
				# Propagate.
				for current_time_step_fine = 1:number_of_time_steps_fine_per_course
						x = x + sigma_fine * randn()
						y = y + sigma_fine * randn()
				end

				# Form image.
				ind_x = convert(Int64, ceil(mod(x, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
				if ind_x >= 1 && ind_x <= number_of_pixels
					ind_y = convert(Int64, ceil(mod(y, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
					if ind_y >= 1 && ind_y <= number_of_pixels
						C_postbleach[ind_x, ind_y, current_frame] += 1
					end
				end
			end
		end
	end

	return cat(dims = 3, C_prebleach, C_postbleach)
end
