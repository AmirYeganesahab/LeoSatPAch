function display(b)

fprintf('Standard orbit interval:\n');
fprintf('Beginning:  %d %d %d %d %d %.1f\n', get(b.dTbeg,'year'), get(b.dTbeg,'month'), get(b.dTbeg,'day'), get(b.dTbeg,'hour'), get(b.dTbeg,'min'), get(b.dTbeg,'sec'));
fprintf('End:  %d %d %d %d %d %.1f\n', get(b.dTend,'year'), get(b.dTend,'month'), get(b.dTend,'day'), get(b.dTend,'hour'), get(b.dTend,'min'), get(b.dTend,'sec'));
fprintf('Polynomial order:  %d\n', b.iPolyOrder);
fprintf('Fit interval:  %.1f hours\n', b.dFitInt/3600);
fprintf('Validity interval:  %.1f hours\n', b.dInterval/3600);

