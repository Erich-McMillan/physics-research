% ChungLuDouble
% Clear all previous work
clear all; clc; close all; %#ok<CLSCR>

% load data
ensemble = dlmread( '/Users/erich/Documents/CollegeWork/Research/C Programming/ChungLuDouble/ChungLuDoubleEnsemble.txt' );
expected = dlmread( '/Users/erich/Documents/CollegeWork/Research/C Programming/ChungLuDouble/ChungLuDoubleProb.txt' );
removed = dlmread( '/Users/erich/Documents/CollegeWork/Research/C Programming/ChungLuDouble/ChungLuDoubleRemoved.txt' );
data = dlmread( '/Users/erich/Documents/CollegeWork/Research/C Programming/ChungLuDouble/ChungLuDoubleData.txt' );

percentRemoved = data(1);
C = data(2);
gamma = data(3);
n = data(4);
y = @(x) C * x.^gamma;

sumEN = sum(ensemble);
sumEX = sum(expected);
sumRE = sum(removed);

figure(1);
hold on;
%plot( ensemble./sumEN, 'ro' );
plot( expected./sumEX, 'b*' );
plot( removed./sumRE, 'mo' );
fplot( @(x) y(x), [.5,n-.5], '--b' );

fun = y(1:1:n);
