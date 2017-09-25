% will determine the degree sequences of each adj matrix in each file
% all files must be adj files
clear size temp_adj degree;
files = dir('AdjStrg_316*');
size = 316;

for f = 1:1:length(files)
    file = files(f);
    fn = file.name;
    adjmat = dlmread(fn);
    
    j = 1;
    for e = 1:1:2500
        temp_adj = adjmat((e-1)*size+1:e*size,1:size);
        
        for k = 1:1:size
           tempd(k,1) = sum(temp_adj(k,:));
          
        end
        tempd = flip(sort(tempd));
        degree((j-1)*size+1:j*size,1) = tempd;
        j = j+1;
    end
    tempn = strcat('degseq_',fn);
    dlmwrite(tempn, degree);
    
end

clear fn adjmat e file j k temp_adj size tempn ans f files degree tempd