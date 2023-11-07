function psi = polyLagrange(k)
%polyLagrange - builds Lagrange shape functions of degree k
%the end nodes are located at x=-1 and x=1
%struct psi with fields .fun and .der
%.fun holds the shape functions
%.der holds the derivatives of the shape functions

for i=1:k+1 %loop to find all the shape functions
    x = linspace(-1,1, k+1); % equally sets nodes between -1 and 1
    y = zeros(1,k+1); %initializes our values
    y(i) = 1; %value of 1 at the element corresponding to the shape function
    psi(i).fun = polyfit(x,y,k); %fitting polynomial to points
    psi(i).der = polyder(psi(i).fun); %taking derivative of found polynomial
end
end

