%___________________________________________________________________%
%  The Single-parameter Bees Algorithm (BA1) source codes           %
%  version 1.0                                                      %
%                                                                   %
%  Developed in MATLAB R2022a                                       %
%                                                                   %
%  Programmer: Hamid Furkan Suluova                                 %
%  Authors: H.F. Suluova and D.T. Pham                              %
%         e-Mail: hfsuluova@gmail.com                               %
%                 hfs158@bham.ac.uk                                 %
%___________________________________________________________________%

clc
clear

IndependentRuns = 1; % Number of Runs
Function = 'F1';     % Benchmark Function
% Choose from {"F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10",
% "F11", "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19", "F20",
% "F21", "F22", "F23"}

PopulationSize = 100;
MaximumNFE = 500000;

for tour=1:IndependentRuns
[Result_BA1(tour).it , Result_BA1(tour).OptCost, Result_BA1(tour).NFE] = Function_BA1(Function, PopulationSize, MaximumNFE);

Result_BA1(tour).BestCost = Result_BA1(tour).OptCost(end);
end

