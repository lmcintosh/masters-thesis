function [p, h, Mut_Info]=muti(a,b,binsA,binsB,minB,maxB) %returning HA
% function to compute the joint entropy H
% and mutual information Mut_Info of two
% variables a and b

%% pA, pB, and pAB
pA=zeros(binsA, 1);
pB=zeros(binsB, 1);    
pAB = zeros(binsA,binsB);

%% Mins, maxs, and binwidths  
a=a';
minA=min(a);   
spikeA = 20; %A is neuron voltage!
k=find(a<spikeA);
maxA=max(a(k));  
binsA = binsA - 1;
binwidthA=(maxA-minA)/(binsA);    % it's binsA-1 otherwise would have binsA + 1 bins

%minB=min(min(b));               
%maxB=max(max(b));               
binwidthB=(maxB-minB)/(binsB);    % it's binsB-1 otherwise would have binsB + 1 bins  


%% Counting


for i = 1:length(a)
    index_a = floor((a(i) - minA)/binwidthA) + 1;
    if (a(i) == maxA)
        index_a = binsA;
    end;
    if (a(i) == spikeA)
        index_a = binsA +1;
    end;
    index_b = floor((b(i) - minB)/binwidthB) + 1;
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

p = [pA pB pAB];
h = [HA HB HAB];
Mut_Info = HA + HB - HAB;

