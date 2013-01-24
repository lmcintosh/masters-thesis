function x=simple_process(T,M,q)
% q is probability of staying, 1 - q is probability of switching


%% make a markov process with x \in {0,1} and p(x_t+1 = x_t|x_t) = q

x=zeros(T,M);
x(1,:) = round(rand(1,M));
for k=1:M
for j=2:T
	a=rand(1,1);
	if (a<=q)
		x(j,k) = x(j-1,k);
	else
		x(j,k) = 1-x(j-1,k);
	end;
end;
end;