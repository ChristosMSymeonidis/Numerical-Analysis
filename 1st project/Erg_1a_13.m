function Erg_1a_13
M= 8; %gia tyxaia M
B = randn(M);
det(B) - rec_det(B) %prepei na einai 0!
k=1;
for N=2:10
    trC=0;
    tC=0;
    tG=0;
    for i=1:10
       A = randn(N);
       b = randn(N,1);
       tic;
       xrC = rec_Cramer(A,b);
       trC = toc + trC;
       tic;
       xC = Cramer(A,b);
       tC = toc + tC;
       tic;
       xG = gaussianElimination(A,b);
       tG = toc + tG;
   end
   tarC(k) = trC/10;
   taC(k) = tC/10;
   taG(k) = tG/10;
   k=k+1;
end
plot(2:10,tarC,2:10,taC,2:10,taG);grid on; 
legend('Rec Crammer','Crammer','Gauss');
 xlabel('matrix dimension NxN') 
 ylabel('tC,tG,trC') 

end

function value = rec_det(A) 
    if size(A,2) == 2  %N==2  
    value = A(1,1)*A(2,2)-A(1,2)*A(2,1);
	else
    value=0;
		for i = 1:size(A,2)%1:N
        M = A(1,i);%stoxeio ij
        sa = A(2:size(A,2),1:size(A,2));% petaw tin prwtti grammi 
        sa(: , i)=[];% petaw kai tin ekastote stili 
        value = value + (-1)^(1+i)*M*rec_det(sa);
        end
	end 
end

function x = rec_Cramer(A, b)


    n = length(b);
	d = rec_det(A); 
    
	x = zeros(n, 1);
	for j = 1:n
		x(j) = rec_det([A(:,1:j-1) b A(:,j+1:end)]) / d;
        
	end
end

function x = gaussianElimination(A, b)
	[m, n] = size(A);
	if m ~= n
		error('Matrix A must be square!');
	end
	n1 = length(b);
	if n1 ~= n
		error('Vector b should be equal to the number of rows and columns of A!');
	end
	Aug = [A b]; % build the augmented matrix
	C = zeros(1, n + 1);
	
	% elimination phase
	for k = 1:n - 1
		% ensure that the pivoting point is the largest in its column
		[pivot, j] = max(abs(Aug(k:n, k)));
		C = Aug(k, :);
		Aug(k, :) = Aug(j + k - 1, :);
		Aug(j + k - 1, :) = C;
		if Aug(k, k) == 0
			error('Matrix A is singular');
		end
		for i = k + 1:n
			r = Aug(i, k) / Aug(k, k);
			Aug(i, k:n + 1) = Aug(i, k:n + 1) - r * Aug(k, k: n + 1);
		end
	end
	
	% back substitution phase
	x = zeros(n, 1);
	x(n) = Aug(n, n + 1) / Aug(n, n);
	for k = n - 1:-1:1
		x(k) = (Aug(k, n + 1) - Aug(k, k + 1:n) * x(k + 1:n)) / Aug(k, k);
	end
end

function x = Cramer(A, b)
    n = length(b);
	d = det(A); 
    % d = rec_det(A);
	x = zeros(n, 1);
	for j = 1:n
		x(j) = det([A(:,1:j-1) b A(:,j+1:end)]) / d;
        % x(j) = rec_det([A(:,1:j-1) b A(:,j+1:end)]) / d;
	end
end

