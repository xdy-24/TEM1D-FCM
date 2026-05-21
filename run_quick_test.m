% run_quick_test.m
% Quick test for TEM1D-FCM
%
% This script runs the example inversion using the provided realmodel.txt
% data file and the main program inv_main.m.

clear;
clc;
close all;

addpath(genpath(pwd));

fprintf('Running quick test for TEM1D-FCM...\n');
fprintf('This test uses the example data file realmodel.txt.\n');

if ~isfile('realmodel.txt')
    error('The file realmodel.txt was not found in the current folder.');
end

if ~isfile('inv_main.m')
    error('The file inv_main.m was not found in the current folder.');
end

inv_main;

fprintf('Quick test completed successfully.\n');
