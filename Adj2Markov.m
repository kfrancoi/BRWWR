function M = Adj2Markov(A)
%Generate a markov chain transition matrix from an adjacency matrix
    n = size(A,1);
    nn = n*n;
    A2 = reshape([1:nn*nn], nn, nn)';
    A3 = reshape([1:n*n],n,n)';
    M = zeros(nn, nn);
    for i = 1:nn;
        for j = 1:nn;
            [r c] = find(A2==i*j);
            [fr fc] = find(A3 == i);
            [tr tc] = find(A3 == j);
            r1 = fr;
            c1 = fc;
             
            if (c1==1 && r1==1)
                if(~isempty(find([i+1 i+n]==j)))
                    M(i,j) = 1;
                end;
                
            elseif (c1==1 && r1==n)
                if(~isempty(find([i-n i+1]==j)))
                    M(i,j) = 1;
                end;
                
            elseif (c1 == n && r1==1)
                if(~isempty(find([i+n i-1]==j)))
                    M(i,j) = 1;
                end;
                
            elseif (c1 == n && r1==n)
                if(~isempty(find([i-1 i-n]==j)))
                    M(i,j) = 1;
                end;
                
            elseif (c1==1)
               if (~isempty(find([i-n i+1 i+n]==j)))
                   M(i,j) = 1;
               end;
               
            elseif (c1==n)
                if (~isempty(find([i-n i-1 i+n]==j)))
                    M(i,j) = 1;
                end;
                
            elseif (r1==1)
                if(~isempty(find([i-1 i+n i+1]==j)))
                    M(i,j) = 1;
                end;
                
            elseif (r1==n)
                if(~isempty(find([i+1 i-n i-1]==j)))
                    M(i,j) = 1;
                end;
                
            else
                if (~isempty(find([i-n i-1 i+1 i+n]==j)))
                	M(i,j) = 1;
                end;
            end;
        end;
    end;