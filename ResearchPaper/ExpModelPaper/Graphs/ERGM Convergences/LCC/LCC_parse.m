% takes the data, enss, and size variables and uses them to parse for data
count = 1;
output = zeros(316,1);
for i = 1:1:enss
    
    for j = 1:1:size
        if isnan(data(count))
            curr = 0;
        else 
            curr = data(count);
        end
        output(j) = output(j) + curr;
        count = count + 1;
    end
    
    
end

output = output./enss;
clear count i j curr