% will find the avg lcc as a function of degree
degseqs = dir('data_degreeseq/degseq_AdjStrg_32**-2.000*');
lccs = dir('data_lcc/LCC_N32_G-2.00*');
size = 32;
lcc_freq = zeros(10,size);
lcc_val = zeros(10,size);
j = 1;

for i = [2,3,4,5,6,7,8,9,10,1]
   cd data_degreeseq;
   tseqn = degseqs(i).name;
   tseq = dlmread(tseqn);
   cd ..;
   
   cd data_lcc;
   tlccn = lccs(i).name;
   tlcc = dlmread(tlccn);
   tlcc(isnan(tlcc))=0;
   cd ..;
   
   for k = 1:1:length(tseq)
        j = tseq(k)+1;
        lcc_freq(i,j) = lcc_freq(i,j) + 1;
        lcc_val(i,j) = lcc_val(i,j) + tlcc(k);
   end 
end

avglcc = lcc_val./lcc_freq;

clear lcc_val lcc_freq i j k lccs size tlcc tlccn tseq tseqn degseqs
disp('still need standard deviation')

