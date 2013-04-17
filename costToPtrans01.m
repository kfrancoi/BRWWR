function P0 = costToPtrans01(C0,eps);
% - C0 is a square cost matrix.
%   The elements of C0 are positive and must be C0(k,l) > eps (0.00000001).
%   If it is impossible to jump from node k to node l, C0(k,l) = realmax,
%   which corresponds to an infinite cost.
%   Node j is supposed to be reachable from node i.
%   Node j is made absorbing (once we have reached this node, we cannot
%   escape and no transition ever occur).
%
% - eps is the error tolerance (default: 1000000 * realmin).
%
% Returns P0: the transition matrix whose elements are proportional to 1/c0(k,l).	

	[nr,nc]   = size(C0);
	
	if (nr ~= nc)
	    fprintf('Error: The cost matrix is not square !\n');
	    return;
	end;
	
	myMin = min(min(C0));

    if (myMin <= eps)
        fprintf('Error: The cost matrix contains elements less than %g !\n',eps);
        return;
    end;
	
	% Computation of the affinity matrix A
	% eps = 1000000 * realmin;
	A   = zeros(nr,nc);
	A(C0 > eps)  = 1./(C0(C0 > eps));
	A(A <= eps) = 0;

	% Computation of P0, the reference transition probabilities matrix
	e  = ones(nr,1);
	P0 = A;
	Den = sum(P0,2)*e';
	P0(Den > eps) = P0(Den > eps)./Den(Den > eps);
	P0(P0 <= eps) = 0;
end
