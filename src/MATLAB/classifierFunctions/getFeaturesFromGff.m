function [genes,features,essential] = getFeaturesFromGff(gff)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%CDSs (stored as struct elements)
features.genes=find(strcmp(gff(:,3),'gene')); %get all genes from list
genes.chr=gff(features.genes,1); %get the chromosome where the genes belong to
genes.coordinates=cell2mat(gff(features.genes,[4 5])); %get the start and end basepair numbers of the genes
genes.annotation=gff(find(strcmp(gff(:,3),'gene')),9); %% add annotation to genes. struct

%%%essential CDSs
features.essential=find(strcmp(gff(:,2),'YeastMine'));
essential.chr=gff(features.essential,1);
essential.coordinates=cell2mat(gff(features.essential,[4 5]));
essential.annotation=cell2mat(gff(features.essential,9));
end