function [ tridist, sqrdist ] = loadTriSqr(files)
% loadTriSqr takes a list of files and loads the number of squares and
%   triangles from each file.
%   each file's data is loaded then appended into a single distribution over
%   the ensemble into tridist and sqrdist respectively
    tridist = [];
    sqrdist = [];

    for i = 1:1:size(files)
       name = files(i).name;
       data = dlmread(name);
       tridist = [tridist; data(:,1)'];
       sqrdist = [sqrdist; data(:,2)'];
    end    

end

