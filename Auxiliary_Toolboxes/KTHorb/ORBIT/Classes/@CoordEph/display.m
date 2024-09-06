function display(b)
for i = 1 : size(b.iPRN,2)
    fprintf('Epoch \n');
    b.dTime(i)
    for j = 1 : size(b.iPRN,1)
        if b.iPRN(j,i) == 0
            continue;
        end
        fprintf('Satellite number %d\n', b.iPRN(j,i));
        fprintf('X-coordinate [m]: %11.3f\n', b.dX(j,i));
        fprintf('Y-coordinate [m]: %11.3f\n', b.dY(j,i));
        fprintf('Z-coordinate [m]: %11.3f\n', b.dZ(j,i));
        fprintf('Clock correction: %g\n', b.dDts(j,i));
    end
end
