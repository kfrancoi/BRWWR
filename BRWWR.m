function P = BRWWR(A, C, theta)
% - A is a square adjacency matrix.
%   The elements of A are positive and represent affinities between nodes.
%   If it is impossible to jump from node k to node l, A(k,l) = 0,
%   which corresponds to an infinite cost.
%   Each node j is supposed to be reachable from each node i.
%
% 
% - theta  = 1;       %% theta must lie between eps (= 0.00000001) and 20.0.
%                     %% If theta = 0, we obtain the expected cost between
%                     %% i and j.
%                     %% If theta = INFINITY, we obtain the shortest-path
%                     kernel.
%
% Returns R: the new transition probability matrix :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps   = 1000000 * realmin;
myMax = realmax;

[nr,nc] = size(A);

if (nr ~= nc)
    fprintf('Error: The cost matrix is not square !\n');
    return;
end;

if (theta < eps) || (theta > 20.0)
    fprintf('Error: The value of theta is out of the admissible range [%g,20.0] !\n',eps);
    return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %e = ones(nr,1);
    %I = eye(nr);
    %H = I - e*e'/nr;
    
    % Computation of the cost matrix Ca (inverse of affinities)
    Ca  = A;
    Ca(A >= eps) = 1./(A(A >= eps));
    Ca(Ca < eps)   = myMax;
    A(A < eps) = 0;

    % Computation of P, the reference transition probabilities matrix
    % representing the natural random walk on the graph
    %Pref = costToPtrans01(Ca,eps);
    %Pref = AtoP(A);
    Pref = A;
    % Computation of the W matrix
    W = exp(-theta * C) .* Pref;
    
    %Computation of the eigenvectors and eigenvalues of W
    % DO I NEED TO SCALE THIS ?
    %[Vl,Dl] = eig(W.');
    [Vr,Dr] = eig(W);
    
    % Test decomposition OK : W*Vr(:,1) == Dr(1,1) * Vr(:,1)
    
    for i=1:nr;
        for j=1:nr;
            P(i,j) = ( Vr(j,1) * W(i,j) ) / ( sum(Vr(:,1) .* W(i,1:nr)') );
        end
    end
    
    %P = ( Vr(1,1:nr) * W(1:nr,1:nr) ) / (sum(Vr(1,:) * W(1:nr,:)));
    
end
    
    