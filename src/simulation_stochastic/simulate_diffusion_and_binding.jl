workspace()

function simulate_diffusion_and_binding()
	# Inititalization of random number generation device.
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	const t_start::Int64 = convert(Int64, time_ns())
	
	delta_t::Float64 = 0.25 # s
	number_of_post_bleach_images::Int64 = 10
	number_of_pixels::Int64 = 128 #  pixels
	number_of_pad_pixels::Int64 = 128 #  pixels
	r_bleach_region::Float64 = 32 #  pixels
	
	number_of_pixels_float::Float64 = convert(Float64, number_of_pixels)
	number_of_pad_pixels_float::Float64 = convert(Float64, number_of_pad_pixels)
	r_bleach_region_float::Float64 = convert(Float64, r_bleach_region)

	intensity_inside_bleach_region::Float64 = 0.6
	intensity_outside_bleach_region::Float64 = 0.9

	D::Float64 = 200.0 # pixels^2 / s
	k_on::Float64 = 0.05 # 1/s
	k_off::Float64 = 0.5 # 1/s
	mobile_fraction::Float64 = 0.5

	number_of_particles::Int64 = 1000000
	
	lambda_free::Float64 = 1.0 / k_on
	lambda_bound::Float64 = 1.0 / k_off
	p_free_marginal::Float64 = k_off / ( k_on + k_off )
	p_bound_marginal::Float64 = k_on / ( k_on + k_off )

	number_of_time_steps_fine_per_course::Int64 = 10
	delta_t_fine::Float64 = delta_t / convert(Float64, number_of_time_steps_fine_per_course)
	sigma_fine::Float64 = sqrt( 2.0 * D * delta_t_fine )

	p_bound_to_free::Float64 = delta_t_fine / lambda_free
	p_free_to_bound::Float64 = delta_t_fine / lambda_bound
	
	area_inside_bleach_region::Float64 = pi * r_bleach_region_float^2
	integral_inside_bleach_region::Float64 = area_inside_bleach_region * intensity_inside_bleach_region

	area_outside_bleach_region::Float64 = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float)^2 - area_inside_bleach_region
	integral_outside_bleach_region::Float64 = area_outside_bleach_region * intensity_outside_bleach_region

	p_inside_bleach_region::Float64 = integral_inside_bleach_region / ( integral_inside_bleach_region + integral_outside_bleach_region )

	x::Float64 = 0.0
	y::Float64 = 0.0
	is_inside_bleach_region::Bool = false
	is_outside_bleach_region::Bool = false
	is_mobile::Bool = false
	is_free::Bool = false
	ind_x::Int64 = 0
	ind_y::Int64 = 0
	
	image_data_post_bleach::Array{Int64, 3} = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images)

	for current_particle = 1:number_of_particles
		if mod(current_particle, 1000000) == 0
			println(current_particle)
		end
		
		# Find random initial position.
		if rand() <= p_inside_bleach_region
			is_inside_bleach_region = false
			while !is_inside_bleach_region
				x = number_of_pad_pixels_float + 0.5 * number_of_pixels_float - r_bleach_region_float + 2.0 * r_bleach_region_float * rand()
				y = number_of_pad_pixels_float + 0.5 * number_of_pixels_float - r_bleach_region_float + 2.0 * r_bleach_region_float * rand()
				
				if (x - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 + (y - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 <= r_bleach_region_float^2
					is_inside_bleach_region = true
				end
			end
		else 
			is_outside_bleach_region = false
			while !is_outside_bleach_region
				x = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()
				y = (number_of_pixels_float + 2.0 * number_of_pad_pixels_float) * rand()
				
				if (x - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 + (y - (number_of_pad_pixels_float + 0.5 * number_of_pixels_float) )^2 > r_bleach_region_float^2
					is_outside_bleach_region = true
				end
			end
		end
		
		# Perform particle motion.
		is_mobile = rand() <= mobile_fraction
		
		if is_mobile
			# Generate "hidden" Markov chain and observations.
			if rand() <= p_free_marginal
				is_free = true
			else
				is_free = false
			end
			
			for current_image_post_bleach = 1:number_of_post_bleach_images
				for current_time_step_fine = 1:number_of_time_steps_fine_per_course
					if is_free
						x = x + sigma_fine * randn()
						y = y + sigma_fine * randn()
					end
					
					if is_free
						if rand() <= p_free_to_bound
							is_free = false
						end
					else
						if rand() <= p_bound_to_free
							is_free = true
						end
					end
					
				end
				
				ind_x = convert(Int64, ceil(mod(x, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
				if ind_x >= 1 && ind_x <= number_of_pixels
					ind_y = convert(Int64, ceil(mod(y, number_of_pixels_float + 2.0 * number_of_pad_pixels_float) - number_of_pad_pixels_float))
					if ind_y >= 1 && ind_y <= number_of_pixels
						image_data_post_bleach[ind_x, ind_y, current_image_post_bleach] += 1
					end
				end
			end
		else # Not mobile
			ind_x = convert(Int64, ceil(x - number_of_pad_pixels_float))
			if ind_x >= 1 && ind_x <= number_of_pixels
				ind_y = convert(Int64, ceil(y - number_of_pad_pixels_float))
				if ind_y >= 1 && ind_y <= number_of_pixels
					for current_image_post_bleach = 1:number_of_post_bleach_images
						image_data_post_bleach[ind_x, ind_y, current_image_post_bleach] += 1
					end
				end
			end
		end        
	end

	file_name_output::String = join(("simulated_frap_data_", string(random_seed), ".dat"))
	file_stream_output::IOStream = open(file_name_output, "w")
	write(file_stream_output, image_data_post_bleach)
	close(file_stream_output)
	
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(t_exec/1e9)

	nothing
end

simulate_diffusion_and_binding()