%% Clear the workspace

clear all
close all
clc

fclose('all');

%% Read in the subject ID data from the text file "BRUN.txt"

fid = fopen('BRUN.txt', 'rt');
C = textscan(fid, '%s %s');
fclose(fid);

PrismaNumber = C{1};
BRUN         = C{2};

N = size(PrismaNumber, 1);

%% Write out the data to an XLSX file with a suitable header, as in early 2019

Head = { '3T Number', 'BRUN', 'Genotype', 'Genotype' };

Type = repmat({ 'Negative' }, [N, 2]);

Data = horzcat(PrismaNumber, BRUN, Type);

Full = vertcat(Head, Data);

xlswrite('GenScan II Subjects Directory with Genotype - Extended.xlsx', Full, 'GenScan II Subjects');
