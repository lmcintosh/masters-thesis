function [h, Mut_Info]=muti2(a,b,binsB,minB,maxB) %returning HA
% function to compute the joint entropy H
% and mutual information Mut_Info of two
% variables a and b

%% pA, pB, and pAB
pA=zeros(2, 1);
pB=zeros(binsB, 1);    
pAB = zeros(2,binsB);

%% Mins, maxs, and binwidths  
a=a';
minA=min(a);   
spikeA = 20;        % A is neuron voltage!
              
binwidthB=(maxB-minB)/(binsB);      


%% Counting


for i = 1:length(a)
    index_a = 1;%floor((a(i) - minA)/binwidthA) + 1
    if (a(i) == spikeA)
        index_a = 2;
    end;
    index_b = floor((b(i) - minB)/binwidthB) + 1;
    if (b(i) == 0)
        index_b = 1;        % HACK for when current starts at zero
    end
    if (b(i) == maxB)
        index_b = binsB;
    end;
    pA(index_a) = pA(index_a) + 1;
    pB(index_b) = pB(index_b) + 1;
    pAB(index_a,index_b) = pAB(index_a,index_b) + 1;
end

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

