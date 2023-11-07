function [K,F] = enforceBC(K,F,boundaryNodes, boundaryValues)
    %purely Dirichlet case
    if(length(boundaryValues('Dirichlet'))==2)
        bv = boundaryValues('Dirichlet');
        K(1,:) = zeros([length(K),1]); %making necessary modifications to K
        K(end,:) = zeros([length(K),1]);
        K(:,1) = zeros([1,length(K)]);
        K(:,end) = zeros([1,length(K)]);
        K(1,1) = 1;
        K(end,end) = 1;
        F = F-(bv(1).*K(:,1))-(bv(2).*K(:,length(F))); %modifying F for the boundary values
        F(1) = bv(1);
        F(end) = bv(2);
    end
    % the one Neumann case
    if(length(boundaryValues('Dirichlet'))==1)
        bnodes = boundaryNodes('Dirichlet');
        if(bnodes(1)==1) %if first BC is Dirichlet
            bv = boundaryValues('Dirichlet');
            K(1,:) = zeros([length(K),1]); %modifying the first row and column of K
            K(:,1) = zeros([1,length(K)]);
            K(1,1) = 1;
            F = F-(bv(1).*K(:,length(K))); %modifying F
            F(1) = bv(1);
            F(end) = F(end) + bnodes(1) * boundaryValues('Neumann'); %adding the condition from the by parts integration
        else
            bv = boundaryValues('Dirichlet');
            K(end,:) = zeros([length(K),1]); %modifying the last row and column of K
            K(:,end) = zeros([1,length(K)]);
            K(end,end) = 1;
            F = F-(bv(1).*K(:,length(K))); %modifying F
            F(end) = bv(1);
            F(1) = F(1) - bnodes(1) * boundaryValues('Neumann'); %subtracting the condition from the by parts integration
        end
    end

end