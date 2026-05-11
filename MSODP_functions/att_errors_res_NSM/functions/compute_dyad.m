function [dyad,dyad2,dyad_vec,dyad2_vec] = compute_dyad(vector)
    %COMPUTE_DYAD for a 3 x N vector and the squared values;
    Nt = length(vector(1, :));
    dyad      = nan(3, 3, Nt);
    dyad2     = nan(3, 3, Nt);
    dyad_vec  = nan(9, Nt);
    dyad2_vec = nan(9, Nt);
    for k = 1:Nt
        x1 = vector(1,k); x2 = vector(2,k); x3 = vector(3,k);
        dyad(:, :, k)   = [0,-x3,x2;x3,0,-x1;-x2,x1,0];
        dyad2(:, :, k)  = dyad(:, :, k) * dyad(:, : , k);
        dyad2_vec(:, k) = reshape(dyad2(:, :, k)', [9, 1]);
        dyad_vec(:, k)  = reshape(dyad(:, :, k)',  [9, 1]);
    end
end

