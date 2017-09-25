% will plot the lagrangian multipliers 
% found by Erich vs the predicted values
er = '/Users/erich/Documents/Programming/Research/Bassler/Summer 2015-Summer 2016/ExponentialModel/LagrangianSolver/results_2017_08_04/';
er20s = 'lagrMultiN316r2.0Sq*';
er25s = 'lagrMultiN316r2.5Sq*';
er30s = 'lagrMultiN316r3.0Sq*';
er35s = 'lagrMultiN316r3.5Sq*';
%lagrMultiN316r2.0Sq90.txt

cd(er);
er20f = dir(er20s);
er25f = dir(er25s);
er30f = dir(er30s);
er35f = dir(er35s);

figure(1);
hold on;
% figure(2);
% hold on;
% figure(3);
% hold on;
% figure(4);
% hold on;
for i = 1:1:1000
    cd(er);
    er20_name = er20f(i).name;
    disp(er20_name);
    er20_data = dlmread(er20_name);
%     figure(1);
    subplot(2,2,1);
    hold on;
    plot(er20_data(:,1),er20_data(:,2));
    
    er25_name = er25f(i).name;
    disp(er25_name);
    er25_data = dlmread(er25_name);
%     figure(2);    
    subplot(2,2,2);
    hold on;
    plot(er25_data(:,1),er25_data(:,2)); 
    
    er30_name = er30f(i).name;
    disp(er30_name);
    er30_data = dlmread(er30_name);
%     figure(3);
    subplot(2,2,3);
    hold on;
    plot(er30_data(:,1),er30_data(:,2));
    
    er35_name = er35f(i).name;
    disp(er35_name);
    er35_data = dlmread(er35_name);
%     figure(4);
    subplot(2,2,4);
    hold on;
    plot(er35_data(:,1),er35_data(:,2));

end

for i = 1:1:1000
    cd(er);
    er20_name = er20f(i).name;
    disp(er20_name);
    er20_data = dlmread(er20_name);
%     figure(1);
    subplot(2,2,1);
    hold on;
    plot(er20_data(:,1),er20_data(:,5),'k');
    
    er25_name = er25f(i).name;
    disp(er25_name);
    er25_data = dlmread(er25_name);
%     figure(2);
    subplot(2,2,2);
    hold on;
    plot(er25_data(:,1),er25_data(:,5),'k');
    
    er30_name = er30f(i).name;
    disp(er30_name);
    er30_data = dlmread(er30_name);
%     figure(3);
    subplot(2,2,3);
    hold on;
    plot(er30_data(:,1),er30_data(:,5),'k');

    
    er35_name = er35f(i).name;
    disp(er35_name);
    er35_data = dlmread(er35_name);
%     figure(4);
    subplot(2,2,4);
    hold on;
    plot(er35_data(:,1),er35_data(:,5),'k');

    
end