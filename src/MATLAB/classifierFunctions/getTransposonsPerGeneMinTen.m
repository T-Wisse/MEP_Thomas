function [tnPerGeneMinTen,readPerGeneMinTen,tnDensityMinTen] = getTransposonsPerGeneMinTen(T,geneLength,geneCount)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

tnPerGeneMinTen = zeros(geneCount,1);
readPerGeneMinTen = zeros(geneCount,1);

for ii=1:geneCount %Loop through each gene. May have an issue with genes on seperate chromosomes here mapping to the same location...
    insertionLocations = str2num(T.InsertionLocations{ii});
    readsPerInsertionLocations = str2num(T.ReadsPerInsertionLocation{ii});
    %Find all reads (tncoordinates_copy) that start at a basepair number between the start and end of a gene (start_coor = gene.coordinates(ii,1) and end_coor = gene.coordinates(ii,2))
    xx=insertionLocations>=T.StartLocation(ii)+geneLength(ii)*0.1 & insertionLocations<=T.EndLocation(ii)-geneLength(ii)*0.1;
    %determine how many reads there are per gene (index of each read is stored in xx). 
    tnPerGeneMinTen(ii) = sum(xx);
    readPerGeneMinTen(ii) = sum(readsPerInsertionLocations(xx));
end

tnDensityMinTen = tnPerGeneMinTen./geneLength;
end