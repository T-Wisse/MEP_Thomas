function [fitness] = readingFitnessData(standardGeneID)
% readingFitnessData reads in fitness data and aligns it to the feature
% table
%   Detailed explanation goes here
% path = D:\Users\Thomas\Studie\MEP\MEP_Thomas\data
    filename = 'CostanzoFitnessReduced.xlsx'; % Name of the excel file with fitness data - Obtained from Costanzo 2022 - Environmental robustness of the global yeast genetic interaction network
    fullFitnessTable = readtable(filename); % read in excel file as table
    fitnessTable = removevars(fullFitnessTable,{'GeneName','x5DayIncubation'}); % Drop unneccesary columns
    
    fitness = zeros(length(standardGeneID),1);

    for ii = 1:length(fitnessTable.SystematicName)
        if sum(strcmp(fitnessTable.SystematicName{ii}, standardGeneID)) == 1
            temp = find(strcmp(fitnessTable.SystematicName{ii}, standardGeneID));
            fitness(temp) = fitnessTable.x3DayIncubation(ii); %NOTE THAT THIS METHOD OVERWRITES THE FITNESS VALUE BY THE LATEST VALUE IF 2 OR MORE ENTRIES MAP TO THE SAME LOCATION WHICH DOES OCCUR
        end
    end
end