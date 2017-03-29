function val = negloglik_db_px( D, ...
                                k_on, ...
                                k_off, ...
                                mf, ...
                                Ib, ...
                                Iu, ...
                                x_bleach, ...
                                y_bleach, ...
                                r_bleach, ...
                                delta_t, ...
                                number_of_pixels, ...
                                number_of_images, ...
                                number_of_pad_pixels, ...
                                data)

disp([D, k_on, k_off, mf, Ib, Iu])

model = signal_db(  D, ...
                    k_on, ...
                    k_off, ...
                    mf, ...
                    Ib, ...
                    Iu, ...
                    x_bleach, ...
                    y_bleach, ...
                    r_bleach, ...
                    delta_t, ...
                    number_of_pixels, ...
                    number_of_images, ...
                    number_of_pad_pixels);

val = sum( data(:) .* log(model(:)) - model(:) - gammaln(data(:) + 1) );

end

