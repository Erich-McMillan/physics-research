% first order transition
clear all; clc; %#ok<CLSCR>

% load data
data1 = dlmread('/Users/erich/Dropbox/Research/Bassler/ChungLuDouble/ChungLuDoubleGammaTrans.txt');
data2 = dlmread('/Users/erich/Dropbox/Research/Bassler/ChungLuDouble/Output Files/firstordertrans.txt');

% plot data
figure(1);
plot(data1(:,1), data1(:,2));
figure(2);
plot(data2(:,1), data2(:,2));

