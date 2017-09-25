% will read the global clustering coefficent values and place into matlab
% variable
clear gcc;
files = dir('*ccg_N32r3.5*');
for i = 1:1:size(files) 
    filename = files(i).name;
    temp = dlmread(filename, '\n');
    gcc(i,1) = mean(temp);
    gcc(i,2) = std(temp);
end

clear files filename i n temp;