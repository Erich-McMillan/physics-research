function [ enslcc, ensavglcc, ensavglccpres ] = loadLCC(files, prescfiles, ntwsize, enssize)
% loadLCC takes a list of files and the size of the networks and loads all 
%   the LCC values and corresponding degrees also will find the avg values 
%   per degree and std over the ensemble
% Places the lcc values and their corresponding degrees in a 2d array then
%   concatenates all sequences together into enslcc. At the same time the
%   avg values per degree per sequence are calculated and placed into
%   ensavglcc so that the effects over the ensemble may be seen.
%   pres stands for prescribed sequence... 

    enselcc = []; %enslcc = {[seqnum];[degree];[val]}
    ensavglcc = zeros(ntwsize,enssize); % ensavglcc = {[degree];[seqnum]}
    enstotlcc = zeros(ntwsize,enssize);
    ensavglccpres = zeros(ntwsize,enssize); % ensavglccpres = {[degree];[seqnum]}
    enstotlccpres = zeros(ntwsize,enssize); 
    
    for i = 1:1:size(files)
        name = files(i).name;
        disp(name);
        data = dlmread(name);
        name = prescfiles(i).name;
        disp(name);
        sequ = dlmread(name);
        enslcc(i,:,:) = data;
        tot = 1;
         for j = 1:1:enssize
            for k = 1:1:ntwsize
                if(data(tot,1) > 0)
                    ensavglcc(data(tot,1),i) = ensavglcc(data(tot,1),i) + data(tot,2)/size(combnk([1:data(1,1)],data(1,1)),2); 
                    enstotlcc(data(tot,1),i) = enstotlcc(data(tot,1),i) + 1;
                    ensavglccpres(sequ(k)) = ensavglcc(data(tot,1),i) + data(tot,2)/size(combnk([1:data(1,1)],data(1,1)),2); 
                    enstotlccpres(sequ(k)) = enstotlcc(data(tot,1),i) + 1;
                    tot = tot + 1;
                end
            end
        end
        
    end
    
    ensavglcc = ensavglcc ./ enstotlcc;
    
end

