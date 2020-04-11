function values_pdf = lappdf(x,y,mx,my,var_lap);

b = sqrt(var_lap/2);

values_x = 1/2/b*exp(-abs(x-mx)/b);

values_y = 1/2/b*exp(-abs(y-my)/b);

values_x_total = prod(values_x(:));
values_y_total = prod(values_y(:));
values_pdf = values_x_total*values_y_total;