function [coordinates,chromosomePerGene,stopChrCoordinates,geneID] = readFileGFF()
%readFileGFF reads the supplied GFF file and returns it content in useful arrays/matrices
%   Detailed explanation goes here

    gffObj = GFFAnnotation('Saccharomyces_cerevisiae.gff3.txt');
    %reference: Chr
    %Start: bp
    %Stop: bp
    %Feature: chr, gene, mRNA, exon
    %strand: + or -
    % frame: no idea
    %attributes: name and more
    
    attributes = gffObj.Attributes;
    startCoordinates = gffObj.Start;  
    stopCoordinates = gffObj.Stop;        
    features = gffObj.Feature;
    chromosomePerGene = gffObj.Reference;
    
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

    coordinates = [startGeneCoordinates, stopGeneCoordinates];
    
    %CAN GET STRAND + or -, could be useful?
    geneAttributes = attributes(geneRows);
    geneID = cell(length(geneAttributes),1);

    for ii = 1:length(startGeneCoordinates)
        partsA = strsplit(geneAttributes{ii}, {'gene:',';biotype'});
        tempGeneID = partsA{2};
        partsB = strsplit(tempGeneID, {'Name='});
        geneID{ii} = partsB{end};
    end


end