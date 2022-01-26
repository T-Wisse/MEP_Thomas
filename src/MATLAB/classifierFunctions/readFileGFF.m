function [coordinates,chromosomePerGeneNumber,stopChrCoordinates,geneID,standardGeneID,strandPerGene] = readFileGFF()
%readFileGFF reads the supplied GFF file and returns it content in useful arrays/matrices
%   coordinates: Begin (:,1) and end (:,2) position of each gene on it's
%   chromosome
%   chromosomePerGeneNumber: Lists the chromosome each gene is located on.
%   Chromosome 17 corresponds with mitochondrial DNA
%   stopChrCoordinates: The length of each chromosome
%   geneID: name of each gene
%   standardGeneID: standardized name of each gene
%   strandPerGene: strand of each gene, forward (1) or reverse (0)


    gffObj = GFFAnnotation('Saccharomyces_cerevisiae.gff3.txt');
    %reference: Chr
    %Start: bp
    %Stop: bp
    %Feature: chr, gene, mRNA, exon
    %strand: + (forward) or - (reverse)
    % frame: no idea
    %attributes: name and more
    
    attributes = gffObj.Attributes;
    startCoordinates = gffObj.Start;  
    stopCoordinates = gffObj.Stop;        
    features = gffObj.Feature;
    chromosomePerGene = gffObj.Reference;
    strand = gffObj.Strand;
       
    geneRows = zeros(length(features),1);
    chrRows = zeros(length(features),1);

    for ii = 1:length(features)
        geneRows(ii) = strcmp(features{ii},'gene');
        chrRows(ii) = strcmp(features{ii}, 'chromosome');
    end
    
    geneRows = logical(geneRows);
    chrRows = logical(chrRows);
    
    stopChrCoordinates = stopCoordinates(chrRows);

    startGeneCoordinates = startCoordinates(geneRows);
    stopGeneCoordinates = stopCoordinates(geneRows);
    chromosomePerGene = chromosomePerGene(geneRows);
    chromosomePerGeneNumber = roman2num(chromosomePerGene); %returns nan for mitochodrial chr.
    chromosomePerGeneNumber(isnan(chromosomePerGeneNumber)) = 17;%set mito chr to 17
    strandPerGene = strand(geneRows);
    strandPerGene = (cell2mat(strandPerGene) == '+');

    coordinates = [startGeneCoordinates, stopGeneCoordinates];
    
    %CAN GET STRAND + or -, could be useful?
    geneAttributes = attributes(geneRows);
    geneID = cell(length(geneAttributes),1);
    standardGeneID = cell(length(geneAttributes),1);

    for ii = 1:length(startGeneCoordinates)
        partsA = strsplit(geneAttributes{ii}, {'gene:',';biotype'});
        temp = partsA{2};
        partsB = strsplit(temp, {'Name='});
        geneID{ii} = partsB{end};
        partsC = strsplit(geneAttributes{ii}, {'gene:',';'});
        standardGeneID{ii} = partsC{2};
    end


end