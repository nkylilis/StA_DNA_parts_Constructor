%% Information

% Author: 
% - Nicolas Kylilis, PhD in Synthetic Biology & Molecular Biology

% Date:
% - 20th Aug 2020

% Description:
% - Computer program that transfroms a part sequence to a StA_Level0 DNA 
% assembly format sequence. 
%   -- Inputs: 1) the part sequence; 2) the StA format prefix&suffix
%   sequences
%   -- Outputs: 2) xls files that are compatible with IDT's online ordering
%   service for oligo synthesis

% Development Environoment & Dependencies
% - Matlab 2020a
% - Bioinformatics Toolbox

clc; clear all; close all;

%% Read input data from file and transform to appropriate data type

filename1 = 'RBS_Library_Calculator_Designs';
T_data = readtable(filename1);
T_data = T_data(:, 1:6);


filename2 = 'StA_prefix_sufffix';
T_StA = readtable(filename2);
% T_StA.Promoter{}
% T_StA.RBS{}
% T_StA.CDS{}
% T_StA.Terminator{}

%% Produce output file with sequences

plate_rows = ['A'; 'B' ;'C'; 'D'; 'E'; 'F'; 'G'; 'H';];

w = 0;
for z = 1:size(plate_rows)
    for x = 1:12
        w = w+1;
        well{w,1} = [plate_rows(z,1)  num2str(x)];
        
    end
    
end


output = {};
for i = 1:size(T_data,1)
    output(i,1) = {well{i,1}};
    output(i,2) = {convertCharsToStrings(T_data.Name(i))};
    output(i,3) = {convertCharsToStrings(T_StA.RBS{1}) + rna2dna(convertCharsToStrings(T_data.RBS_sequence{i})) + convertCharsToStrings(T_StA.RBS{2})};
    output(i,4) = {seqrcomplement(output{i,3})};
end



% **** uncomment as needed *********

% Plate doublex
T_output = array2table(output,...
    'VariableNames',{'Well Position', 'Name','Sequence','Sequence 2'});
writetable(T_output,'StA_parts Plate_dublex.xlsx');

 
% Plate primer pair
% output2 = [output(:,1:3) ; output(:,[1:2,4])];
% T_output = array2table(output2,...
%     'VariableNames',{'Well Position', 'Name','Sequence'});
% writetable(T_output,'StA_parts Plate_primer_pair.xlsx')

