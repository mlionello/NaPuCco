function tot_null_vx = get_nb_null_vx(eff_width, null_width)
tot_null_vx = 6*null_width*eff_width^2;
tot_null_vx = tot_null_vx + 8*null_width^3;
tot_null_vx = tot_null_vx + 12*eff_width*null_width^2;
end