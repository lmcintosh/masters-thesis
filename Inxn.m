function [h, Mut_Info, countsame]=Inxn(a,b,bins,min,max) %returning HA
% function to compute the MI between past and future. countsame is now the
% probability of the bin staying the same in the next time instance.

%% pA, pB, and pAB
pA=zeros(bins, 1);
pB=zeros(bins, 1);    
pAB = zeros(bins,bins);

%% Mins, maxs, and binwidths 
			
binwidth=(max-min)/(bins);    


%% Counting

countsame = 0;
for i = 1:length(a)
    index_a = floor((a(i) - min)/binwidth) + 1;
	if (index_a < 1)
		index_a = 1;        % HACK for when current starts at zero
	end
	if (a(i) == max)
		index_a = bins;
	end;
    index_b = floor((b(i) - min)/binwidth) + 1;
    %if (b(i) == 0)
    %    index_b = 1;        % HACK for when current starts at zero
    %end
    if (b(i) == max)
        index_b = bins;
    end;
	%indexOFa(i)=index_a;

	if(index_a == index_b)
		countsame = countsame + 1;
		%index_a
		%index_b
	end;
    if(index_a > bins || index_b > bins)
		index_a
	    index_b
	end;
    %count
    pA(index_a) = pA(index_a) + 1;
    pB(index_b) = pB(index_b) + 1;
    pAB(index_a,index_b) = pAB(index_a,index_b) + 1;
end

%countsame = countsame / length(a);

sumA = sum(pA);
sumB = sum(pB);
sumAB = sum(sum(pAB));

pA = pA/sumA;
pB = pB/sumB;
pAB = pAB/sumAB;

%% Entropies
HA = 0;
HB = 0;
HAB = 0;

for i = 1:length(pB)
    if pB(i)~= 0
        HB = HB - pB(i)*log2(pB(i));
    end
end

for i = 1:length(pA)
    if pA(i)~= 0
        HA = HA - pA(i)*log2(pA(i));
        for j = 1:length(pB)
            if pAB(i,j)~=0
                HAB = HAB - pAB(i,j)*log2(pAB(i,j));
            end
        end
    end
end

%p = [pA pB pAB];
h = [HA HB HAB];
Mut_Info = HA + HB - HAB;
countsame = countsame/length(a);
