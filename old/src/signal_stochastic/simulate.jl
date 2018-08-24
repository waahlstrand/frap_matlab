function simulate(	D::Float64, 
					k_on::Float64, 
					k_off::Float64, 
					mf::Float64,
					Ib::Float64,
					Iu::Float64,
					r_bleach::Float64,
					number_of_pixels::Int64,
					number_of_images::Int64,
					number_of_pad_pixels::Int64,
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
	
	delta_t_fine::Float64 = delta_t / convert(Float64, number_of_time_steps_fine_per_course)
	sigma_fine::Float64 = sqrt( 2.0 * D * delta_t_fine )

	p_bu::Float64 = delta_t_fine / tau_b
	p_ub::Float64 = delta_t_fine / tau_u
	
	# Pre-work for picking random initial positions.
	area_inside_bleach_region::Float64 = pi * r_bleach^2
	integral_inside_bleach_region::Float64 = area_inside_bleach_region * Ib

	area_outside_bleach_region::Float64 = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float)^2 - area_inside_bleach_region
	integral_outside_bleach_region::Float64 = area_outside_bleach_region * Iu

	p_inside_bleach_region::Float64 = integral_inside_bleach_region / ( integral_inside_bleach_region + integral_outside_bleach_region )

	# Generate particle trajectories and FRAP image data.
	x::Float64 = 0.0
	y::Float64 = 0.0
	is_inside_bleach_region::Bool = false
	is_outside_bleach_region::Bool = false
	is_mobile::Bool = false
	is_unbound::Bool = false
	ind_x::Int64 = 0
	ind_y::Int64 = 0
	
	data::Array{Int64, 3} = zeros(number_of_pixels, number_of_pixels, number_of_images)

	for current_particle = 1:number_of_particles
		if mod(current_particle, 100000) == 0
			println(current_particle)
		end

		# Find random initial position.
		if rand() <= p_inside_bleach_region
			is_inside_bleach_region = false
			while !is_inside_bleach_region
				x = number_of_pad_pixels_float + 0.5 * number_of_pixels_float - r_bleach + 2.0 * r_bleach * rand()
				y = number_of_pad_pixels_float + 0.5 * number_of_pixels_float - r_bleach + 2.0 * r_bleach * rand()
				
				if (x - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 + (y - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 <= r_bleach^2
					is_inside_bleach_region = true
				end
			end
		else 
			is_outside_bleach_region = false
			while !is_outside_bleach_region
				x = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()
				y = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()
				
				if (x - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 + (y - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 > r_bleach^2
					is_outside_bleach_region = true
				end
			end
		end
		
		# Perform particle motion.
		is_mobile = rand() <= mf
		
		if is_mobile
			if rand() <= p_u
				is_unbound = true
			else
				is_unbound = false
			end
			
			for current_image = 1:number_of_images
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
				
				# Add 1 to the the pixel where the particle currently resides.
				ind_x = convert(Int64, ceil(mod(x, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
				if ind_x >= 1 && ind_x <= number_of_pixels
					ind_y = convert(Int64, ceil(mod(y, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
					if ind_y >= 1 && ind_y <= number_of_pixels
						data[ind_x, ind_y, current_image] += 1
					end
				end
			end
		else # Not mobile
			# Add 1 to the the pixel where the particle resides.
			ind_x = convert(Int64, ceil(x - number_of_pad_pixels_float))
			if ind_x >= 1 && ind_x <= number_of_pixels
				ind_y = convert(Int64, ceil(y - number_of_pad_pixels_float))
				if ind_y >= 1 && ind_y <= number_of_pixels
					for current_image = 1:number_of_images
						data[ind_x, ind_y, current_image] += 1
					end
				end
			end
		end        
	end

	return data
end
