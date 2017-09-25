% will find the histogram of the differences between the convergences of
% the lagrange values found by Erich and Weibin.
er = '/Users/erich/Documents/Programming/Research/Bassler/Summer 2015-Summer 2016/ExponentialModel/LagrangianSolver/results_2017_07_03/';
er20s = 'lagrMultiN316r2.0Sq**.txt';
er25s = 'lagrMultiN316r2.5Sq**.txt';
er30s = 'lagrMultiN316r3.0Sq**.txt';
er35s = 'lagrMultiN316r3.5Sq**.txt';

wb = '/Users/erich/Documents/Programming/Research/Bassler/Summer 2015-Summer 2016/WeibinResults/seq316_dbroot_NoSelfLoop/';
wb20s = 'dbroot_N316r2.0Sq*';
wb25s = 'dbroot_N316r2.5Sq*';
wb30s = 'dbroot_N316r3.0Sq*';
wb35s = 'dbroot_N316r3.5Sq*';

final_20_lagr = [];
final_25_lagr = [];
final_30_lagr = [];
final_35_lagr = [];

final_20_erFV = [];
final_25_erFV = [];
final_30_erFV = [];
final_35_erFV = [];

final_20_wbFV = [];
final_25_wbFV = [];
final_30_wbFV = [];
final_35_wbFV = [];

cd(er);
er20f = dir(er20s);
er25f = dir(er25s);
er30f = dir(er30s);
er35f = dir(er35s);

cd(wb);
wb20f = dir(wb20s);
wb25f = dir(wb25s);
wb30f = dir(wb30s);
wb35f = dir(wb35s);

for i = 1:1:1000
    % er
    cd(er);
    er20_name = er20f(i).name;
    disp(er20_name);
    er20_data = dlmread(er20_name);
    er20_lagr = er20_data(:,2);
    
    er25_name = er25f(i).name;
    disp(er25_name);
    er25_data = dlmread(er25_name);
    er25_lagr = er25_data(:,2);
    
    er30_name = er30f(i).name;
    disp(er30_name);
    er30_data = dlmread(er30_name);
    er30_lagr = er30_data(:,2);
    
    er35_name = er35f(i).name;
    disp(er35_name);
    er35_data = dlmread(er35_name);
    er35_lagr = er35_data(:,2);
    
    % wb
    cd(wb);
    wb20_name = wb20f(i).name;
    disp(wb20_name);
    wb20_data = dlmread(wb20_name);
    wb20_lagr = wb20_data(:,2);
    
    wb25_name = wb25f(i).name;
    disp(wb25_name);
    wb25_data = dlmread(wb25_name);
    wb25_lagr = wb25_data(:,2);
    
    wb30_name = wb30f(i).name;
    disp(wb30_name);
    wb30_data = dlmread(wb30_name);
    wb30_lagr = wb30_data(:,2);
   
    wb35_name = wb35f(i).name;
    disp(wb35_name);
    wb35_data = dlmread(wb35_name);
    wb35_lagr = wb35_data(:,2);

    % append to final
    final_20_lagr = horzcat(final_20_lagr, (er20_lagr-wb20_lagr)');
    final_25_lagr = horzcat(final_25_lagr, (er25_lagr-wb25_lagr)');
    final_30_lagr = horzcat(final_30_lagr, (er30_lagr-wb30_lagr)');
    final_35_lagr = horzcat(final_35_lagr, (er35_lagr-wb35_lagr)');
    
    final_20_erFV = horzcat(final_20_erFV, er20_data(:,4)');
    final_25_erFV = horzcat(final_25_erFV, er25_data(:,4)');
    final_30_erFV = horzcat(final_30_erFV, er30_data(:,4)');
    final_35_erFV = horzcat(final_35_erFV, er35_data(:,4)');
    
    final_20_wbFV = horzcat(final_20_wbFV, wb20_data(:,4)');
    final_25_wbFV = horzcat(final_25_wbFV, wb25_data(:,4)');
    final_30_wbFV = horzcat(final_30_wbFV, wb30_data(:,4)');
    final_35_wbFV = horzcat(final_35_wbFV, wb35_data(:,4)');

end
 
clear *name *data *f *s i er20_lagr er25_lagr er30_lagr er35_lagr wb20_lagr wb25_lagr wb30_lagr wb35_lagr wb er
 
 
 
 