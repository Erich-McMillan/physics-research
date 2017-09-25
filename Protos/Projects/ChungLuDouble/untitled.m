% ChungLuDouble
% Clear all previous work
clear all; clc; close all; %#ok<CLSCR>

% load data
ensemble = dlmread( '/Users/erich/Documents/CollegeWork/Research/C Programming/ChungLuDouble/ChungLuDoubleEnsemble.txt' );
expected = dlmread( '/Users/erich/Documents/CollegeWork/Research/C Programming/ChungLuDouble/ChungLuDoubleProb.txt' );

figure(1);
hold on;
plot( ensemble(1:100),'ro-');
plot( expected(1:100), 'bo-' );