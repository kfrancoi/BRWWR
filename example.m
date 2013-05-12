function example
%Launch BRWWR on a small example 

n = 30;
A = ones(n,n);

%Cinit = [1, 1, 1, 1, 3, 1, 1, 1, 1];

%Cinit = [1,2,3
%        4,5,6];

%1 Matrix 1 (Gaussian obstacle)
%Cinit = gauss2d(A, 3, [15,15]);
%Cinit = reshape(Cinit.*5', 1,n*n);

Cinit = ones(n);
Cinit = Cinit(1:20, 10);
Cinit - Cinit(10:30, 20);
Cinit = reshape(Cinit', 1,n*n);

%Cinit = [1,  1,  1,  1,  1,  1,  1,  1,  1 
%      1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  2,  2,  2,  1,  2,  2,  2,  1
%      1,  2,  3,  2,  1,  2,  3,  2,  1
%      1,  2,  2,  2,  1,  2,  2,  2,  1
%      1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  1,  1,  1,  1,  1,  1,  1,  1];
%Cinit = reshape(Cinit', 1,n*n);

%Cinit = [1,  1,  1,  3,  1,  1,  1,  1,  1
%      1,  1,  1,  3,  1,  1,  1,  1,  1
%      1,  1,  1,  3,  1,  1,  1,  1,  1
%      1,  1,  1,  3,  1,  1,  1,  1,  1
%      1,  1,  1,  3,  1,  1,  1,  1,  1
%      1,  1,  1,  3,  1,  1,  3,  1,  1
%      1,  1,  1,  1,  1,  1,  3,  1,  1
%      1,  1,  1,  1,  1,  1,  3,  1,  1
%      1,  1,  1,  1,  1,  1,  3,  1,  1];
%Cinit = reshape(Cinit', 1,n*n);

%Cinit = [1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  1,  5,  1,  1,  1,  1,  1,  1
%      1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  1,  1,  1,  1,  1,  1,  1,  1
%      1,  1,  1,  1,  1,  5,  5,  5,  5
%      1,  1,  1,  1,  1,  5,  5,  5,  5
%     1,  1,  1,  1,  1,  5,  5,  5,  5
%      1,  1,  1,  1,  1,  5,  5,  5,  5];
%Cinit = reshape(Cinit', 1,n*n);

C = ones(n*n, 1) * Cinit;


%Compute Markov chain transition matrix
M = Adj2Markov(A);
M = AtoP(M);
%M = costToPtrans01(ones(n*n,n*n),1000000 * realmin);

%Computation of Markov chain transition probability matrix
%nM = size(C,1)^2;
%MA = zeros(nM,nM);

%for i=1:nM;
%    for j=1:nM;
%        if (abs(mod(i,nM) - mod(j,nM)) == 1)
%            MA(i,j) = 1;
%        else 
%            MA(i,j) = 0;
%        end
%    end
%end

P = BRWWR(M, C, 0.1);

%Stationary disctribution of Pref
[Vl,Dl] = eig(M.');
PrefStat = reshape(-Vl(:,1), n,n)';


%Stationary distribution of BRWWR P
[Vl,Dl] = eig(P.');
PStat = reshape(Vl(:,1), n,n)';

h = figure;

%subplot(2,2,1); 
%imagesc(A);
%colorbar;

%subplot(2,2,2);
imagesc(reshape(Cinit, n, n)');
colorbar;
print(h, '-dpdf', 'A.pdf')

%subplot(2,2,3);
%imagesc(PrefStat);
%colorbar;

j = figure;
%subplot(2,2,4);
imagesc(abs(PStat));
colorbar;

print(j, '-dpdf', 'B.pdf');



