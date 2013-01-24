function [] = AdEx_masterscript(M,T);
%% parameters
a = linspace(0,40,20);
vBins = [2,3,4,10,50,100];
iBins = vBins;
itype = [0,1];
statetype = [0,1];
noisetype = [1,3]; % 0 is step, 1 is ou_318, 2 is o_u, 3 is shot noise, 4 is white noise, 5 is Markov
current = [200, 300, 400, 500, 600, 700, 800];
%M = 1000;

%% outputs
%i_mem = cell(length(a),length(vBins),length(iBins),length(itype),length(statetype),length(noisetype),length(current));
rootname = ['functa_',date];
extension = '.mat';


for j = 1:length(vBins)
    for k = 1:length(iBins)
        for l = 1:length(itype)
            for m = 1:length(statetype)
                for n = 1:length(noisetype)
                    for o = 1:length(current)
                        i_mem = cell(length(a),1);
                        i_pred = i_mem;
                        non_pred = i_mem;
                        structureInI = i_mem;
                        for i = 1:length(a)
                            [i_mem{i},i_pred{i},non_pred{i},structureInI{i}] = AdEx_perm(a(i),M,T,vBins(j),iBins(k),itype(l),statetype(m),noisetype(n),current(o));
                        end
                        data = int2str([j,k,l,m,n,o]);
                        data = data(data~=' ');
                        filename = [rootname, data, extension];
                        everything = struct('i_mem',i_mem,'i_pred',i_pred,'non_pred',non_pred,'structureInI',structureInI);
                        eval(['save ', filename,' everything']);
                        clear i_mem i_pred non_pred structureInI
                    end
                end
            end
        end
    end
end
