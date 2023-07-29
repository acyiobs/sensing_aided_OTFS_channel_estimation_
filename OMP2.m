function [x, omega] = OMP2(b, A, S, supp)
err= 1e-3;


[N,K] = size(A); % N:dim of signal, K:#atoms in dictionary
if (N ~= size(b))
    error('Dimension not matched');
end

x = zeros(K,1); 
if ~supp
    r = b;        
    omega = [];
    A_omega = [];    
else
    x(supp, :) = A(:, supp) \ b;
    r = b - A * x;
    omega = supp;
    A_omega = A(:,supp);
end

while size(omega, 1) < S 
    x_tmp = zeros(K,1);
    inds = setdiff([1:K],omega); 

    x_tmp(inds) = A(:, inds)' * r; 
    [~, ichosen] = max(abs(x_tmp)); 
    omega = [omega; ichosen];
    
    x = zeros(K, 1);
    x(omega, :) = A(:, omega) \ b; 
    r = b - A(:, omega) * x(omega, :); 

    if norm(r) <=err*norm(b)
        break
    end
    
end

end