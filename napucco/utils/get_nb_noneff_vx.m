function tot_noneff_vx = get_nb_noneff_vx(eff_width, noneff_width)
tot_noneff_vx = 6*noneff_width*eff_width^2;
tot_noneff_vx = tot_noneff_vx + 8*noneff_width^3;
tot_noneff_vx = tot_noneff_vx + 12*eff_width*noneff_width^2;
end