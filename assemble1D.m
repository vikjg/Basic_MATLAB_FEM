function [K,F] = assemble1D(K,F,ki,fi,nodelist)
%simply adds the i-th element matrix into the corresponding place (given by
%node list)
    K(nodelist(1):nodelist(2), nodelist(1):nodelist(2)) = K(nodelist(1):nodelist(2), nodelist(1):nodelist(2)) + ki;
    F(nodelist(1):nodelist(2)) = F(nodelist(1):nodelist(2)) + fi;
end