function P = AtoP(A)
%Usage : [P] = AtoP(A)
%Compute the transistion probability matrix from adjacent matrix
%Author : Fran?oisse Kevin 08

% n =size(A,1);
% P = zeros(n,n);
% for i=1:n
%     for j=1:n
%         P(i,j) = A(i,j) / sum(A(i,:));
%     end
% end

s = sum(A,2);
n = size(A,1);
e = ones(n,1);
P = A ./ (s*e');
P = full(P);