% first order transition
clear all; clc; %#ok<CLSCR>

% load data
data1 = dlmread('/Users/erich/Documents/CollegeWork/Research/C Programming/ChungLuDouble/ChungLuDoubleGammaTrans.txt');

% plot data
figure(1);
plot(data1(:,1), data1(:,2));

