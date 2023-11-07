function [ke, fe] = element1D(psi)
    d = length(psi)-1; %setting length and initializing
    ke.k = zeros(d+1);
    ke.b = zeros(d+1);
    fe = zeros([d+1,1]);
    for i = 1:d+1 %big loop
        for j = 1:i %little loop and use symmetry
            nk = ceil((2*d+1)/2); %number of quadrature points
            Qk = quadGaussLegendre(nk); %get the points
            ke.k(i,j) = dot(Qk.Weights,polyval(conv(psi(i).der, psi(j).der),Qk.Points)); %compute the integral using numerical methods
            ke.k(j,i) = ke.k(i,j); %symmetry!
            nb = ceil((2*d+1)/2); %repeat the above process
            Qb = quadGaussLegendre(nb);
            ke.b(i,j) = dot(Qb.Weights,polyval(conv(psi(i).fun, psi(j).fun),Qb.Points));
            ke.b(j,i) = ke.b(i,j);
        end
        n = ceil(2*d + 1/2); %same process as above but for f
        Q = quadGaussLegendre(n);
        fe(i) = dot(Q.Weights,polyval(psi(i).fun,Q.Points)); 
    end
end