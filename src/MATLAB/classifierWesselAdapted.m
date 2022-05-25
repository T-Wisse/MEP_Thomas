%https://sites.google.com/site/satayusers/

clear all
close all

%ORDER OF FILES LOADED: BEDFILE, perInsertionsFile, perGeneFile, essentialGeneFile
% File locations: C:\Users\Thomas\Desktop\Uni\SATAYDATA_211210\SATAY_Data\SATAY_Analysis

%%renamed the WT .wig file to fit this naming scheme. should check if it is
%%from the exact same run. Also renamed yLIC bed file

filebase = 'WT_merged-techrep-a_techrep-b_trimmed.sorted.bam'; % Wild-type for training classifier
% filebase = 'all_cleaned_fw_reads_trimmed.sorted.bam'; %This is yLIC137_7
% filebase = 'yLIC137_8_merged_cleaned_forward_reads_trimmed.sorted.bam';
% filebase = 'yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam';
% filebase = 'yTW001_6_merged_cleaned_forward_reads_trimmed.sorted.bam';
% filebase = 'yWT03a_16_merged_cleaned_forward_reads_trimmed.sorted.bam'; %bem1;
% filebase = 'yWT03a_21_merged_cleaned_forward_reads_trimmed.sorted.bam'; %bem1
% filebase = 'yWT04a_14_merged_cleaned_forward_reads_trimmed.sorted.bam'; %bem1bem3
% filebase = 'yWT04a_23_merged_cleaned_forward_reads_trimmed.sorted.bam'    ; %bem1bem3


%% Load satay data from transposonmapper & annotation file

% Should be able to do this with path and regexp...
perInsertionsFile = append(filebase,'_pergene_insertions.txt');
perGeneFile = append(filebase,'_pergene.txt');
essentialGeneFile = append(filebase, '_peressential.txt');
wigFile = append(filebase, '.wig');
bedFile = append(filebase,'.bed.txt');
[coordinates, chromosomePerGene, chromosomeEndPos, geneID, standardGeneID, geneSenseBinary] = readFileGFF(); %Read in the gff annotation file

% [bedFile, path]=uigetfile('Pick a file'); %<-point to the adjusted .bed
% file ALTERNATIVE METHOD


%pergene_insertions.txt: Gene name	Chromosome	Start location	End location  Insertion locations	Reads per insertion location

Table1 = readtable(perInsertionsFile); %load .insertions file
Table2 = readtable(perGeneFile); %load .pergene file
Table3 = readtable(essentialGeneFile); %load peressential file
data = load(bedFile, '-ascii'); %.bed file: Chromosome name start location end location . score. Note: file must be changed manually to be numerical only
readData = loadWigfile(wigFile); %Takes a long time, can we shorten this somehow?
essentialGenes = readtable('Cerevisiae_AllEssentialGenes_List.txt'); %Read in the essential genes file % THIS IS BASED ON A COMBINED LIST TAKEN FROM: https://rdrr.io/bioc/SLGI/man/essglist.html and http://www.essentialgene.org/

geneCount = height(Table1);
essentialGeneList = zeros(length(standardGeneID),1);
for ii = 1:height(standardGeneID)
    essentialGeneList(ii) = sum(strcmp(standardGeneID{ii},essentialGenes.AllEssentialGenesFoundInLists___Cerevisiae_EssentialGenes_List_)); %For each gene check if it is listed as essential
end
data = data(:,[1,2]); %We are only interested in the first 2 entries for the .bed file.
data(:,1) = transformChromosomeIdentifiers(data(:,1)); % transformation of used chromosome identifiers to 1:17. 

%% put all features and coordinates on one single huge concatenated chromosome
% Why do we have non-unique entries to tn locations? probably they go in reverse directions?

lengthPreviousChr = [0; cumsum(chromosomeEndPos)]; %get the length in bp of each chromosome. The first chromosome starts at 0
lengthPreviousChr = double(lengthPreviousChr(1:end-1)); % Remove the last entry as it does not act as a starting point to another chromosome
%force type double to allow for addition with other double type data. May
%be useful to force double immediately upon reading in, also for other data
tnCoordinates = data;

[tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat] = concatCoordinatesIntoSingleChr(Table1,tnCoordinates,geneCount,lengthPreviousChr);


%% reading in wig file

chrEntries = find(isnan(readData(:,1)));
chrEntries = [chrEntries; length(readData)];

for ii = 1:17 %ARE THEY IN THE SAME ORDER? OTHERWISE THIS DOES NOT MAKE SENSE. Readdata counts up 1133 - 1148. Mito last
    readData(chrEntries(ii):chrEntries(ii+1),1) = readData(chrEntries(ii):chrEntries(ii+1),1) + lengthPreviousChr(ii); 
end
readData(isnan(readData(:,1)),:) = []; %Remove rows which contain a nan value


%% Get number of reads & transposon per gene
tnPerGene = Table2.NumberOfTransposonsPerGene;
readPerGene = Table2.NumberOfReadsPerGene;

geneLength = Table1.EndLocation - Table1.StartLocation; %CAN USE GFF DATA FOR THIS
tnDensity = tnPerGene./geneLength; %define transposon density per gene

% Count number of transposon per gene minus end and beginning
[tnPerGeneMinTen,readPerGeneMinTen,tnDensityMinTen] = getTransposonsPerGeneMinTen(Table1,geneLength,geneCount); %the read count from table1 does not correspond with read count from table2. 

% Find the longest transposon free interval within each gene
tnFreeInterval = getTnFreeIntervalPerGene(tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat,geneCount);

%% intergenic regions of 20kb 
% RUNS, BUT DOES IT actuaLLY WORK AS INTENDED?

[intergenicLength20kbPerGene,intergenicTn20kbPerGene,intergenicTnDensity20kb] = getIntergenicTnDensity(geneLength,chromosomeEndPos,tnCoordinates,coordinates,chromosomePerGene);


%% creating a list with the number of transposons and reads in the 100bp 5' of each gene (promotor region)

geneCoordinatesConcat(:,1) = geneStartCoordinatesConcat;
geneCoordinatesConcat(:,2) = geneEndCoordinatesConcat;

% creating a concatinated verson of all transposons and reads.
readPosConcat = zeros(max(max(geneCoordinatesConcat)),1);
readPosConcat(readData(:,1)) = readData(:,2);
tnPosConcat = zeros(max(max(geneCoordinatesConcat)),1);
tnPosConcat(tnCoordinatesConcat) = 1;

[promotorRead] = promotorData(readPosConcat,geneCoordinatesConcat,geneCount,geneSenseBinary);
[promotorTn] = promotorData(tnPosConcat,geneCoordinatesConcat,geneCount,geneSenseBinary);

promotorReadPerTn=promotorRead./promotorTn;

%% Combine training/testing data into a single table
readnumb = 1; %TEMPORARY SO I CAN CONTINUE. VERY VERY WRONG! WHAT IS READNUMB? %READNUMB IS PROBABLY TOTAL READ COUNT
% 

readPerGenePerTn = readPerGene./(tnPerGene.*sum(readnumb)); %unused?
readPerGenePerTnMinTen = readPerGeneMinTen./(tnPerGeneMinTen.*sum(readnumb)); %unused?
tnFreeIntervalPerGeneLength = (double(tnFreeInterval)./geneLength);
intergenicTn20kbUp = intergenicTn20kbPerGene(:,1);
intergenicTn20kbDown = intergenicTn20kbPerGene(:,2);
NI20kb = tnDensity./((1+intergenicTnDensity20kb(:,1)+intergenicTnDensity20kb(:,2))/2); %IntergenicTnDensity20kb can be 0... --> inf values --> added 1
%NI20kb is compares the density of transposons in the local region to the3
%transposons in the gene
readsPerTnAverage = readPerGenePerTn;

dataTableOfFeatures = table(tnPerGeneMinTen, readPerGeneMinTen, geneLength, tnDensity, tnFreeIntervalPerGeneLength, promotorTn, intergenicTn20kbUp, intergenicTn20kbDown, NI20kb, essentialGeneList);% 

% dataTableOfFeatures_WT_ab = dataTableOfFeatures;

%% Train classifier 
% [A, B, C] = trainClassifier(dataTableOfFeatures_WT) %USE BAGGED TREES !!!!
%%
% yfit = trainedModelWTStripped.predictFcn(dataTableOfFeatures);

% classifier = baggedTreesModel_220125;

%% FOR TRAINING REGRESSION PREDICTOR
fitnessData = readingFitnessData(standardGeneID); % Read in experimental fitness data for a set of genes
% fitnessData(essentialGeneList==1) = 0; %Set essential genes to have 0 fitness % Uncomment when using essential genes as well
dataTableOfFeatures.('readsPerTnMinTen') = readPerGenePerTnMinTen;
dataTableOfFeatures.('fitness') = fitnessData;
dataTableOfFeatures = removevars(dataTableOfFeatures,{'essentialGeneList'}); % remove essential gene list as we do not need it here
selection = (fitnessData>0); % Selection of all genes that have fitness data
% selection = (fitnessData+essentialGeneList)>0; % Selection of all genes that are essential or have fitness data % Uncomment when using essential genes as well
trainingTable = dataTableOfFeatures(selection,:);
% trainingTable = removevars(trainingTable,{'essentialGeneList'}); %Remove this variable as it completely maps all the 0 fitness values but we don't have properly for other datasets)


