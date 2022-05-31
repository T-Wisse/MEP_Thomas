%https://sites.google.com/site/satayusers/

clear all
close all

% I think I fixed chromosome order, but I should double check that because
% it did not seem to matter

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
% filebase = 'yWT04a_23_merged_cleaned_forward_reads_trimmed.sorted.bam'; %bem1bem3


%% Load satay data from transposonmapper & annotation file

% Should be able to do this with path and regexp...
perInsertionsFile = append(filebase,'_pergene_insertions.txt');
perGeneFile = append(filebase,'_pergene.txt');
essentialGeneFile = append(filebase, '_peressential.txt');
wigFile = append(filebase, '.wig');
bedFile = append(filebase,'.bed.txt');
[coordinatesGenesRaw, chromosomePerGeneRaw, chromosomeEndPosRaw, geneIDRaw, standardGeneIDRaw, geneSenseBinaryRaw] = readFileGFF(); %Read in the gff annotation file

% [bedFile, path]=uigetfile('Pick a file'); %<-point to the adjusted .bed
% file ALTERNATIVE METHOD


%pergene_insertions.txt: Gene name	Chromosome	Start location	End location  Insertion locations	Reads per insertion location

Table1Raw = readtable(perInsertionsFile); %load .insertions file
Table2Raw = readtable(perGeneFile); %load .pergene file
% Table3 = readtable(essentialGeneFile); %load peressential file % I DONT USE THIS
data = load(bedFile, '-ascii'); %.bed file: Chromosome name start location end location . score. Note: file must be changed manually to be numerical only
readDataRaw = loadWigfile(wigFile); %Takes a long time, can we shorten this somehow?
essentialGenes = readtable('Cerevisiae_AllEssentialGenes_List.txt'); %Read in the essential genes file % THIS IS BASED ON A COMBINED LIST TAKEN FROM: https://rdrr.io/bioc/SLGI/man/essglist.html and http://www.essentialgene.org/

data = data(:,[1,2]); %We are only interested in the first 2 entries for the .bed file.
data(:,1) = transformChromosomeIdentifiers(data(:,1)); % transformation of used chromosome identifiers to 1:17. 


%% reordering data to follow the same chromosome structure
% chr order per file ARE NOTED IN C:\Users\Thomas\Desktop\SATAY_DATA. They
% are as follows. % NOTE: we transfer to BedFile Order because the order makes the
% most sense.
% chrOrderBedFile = [1:17]; % ONLY BED FILE (data)

% chrOrderGFF = [1, 2, 3, 4, 9, 17, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]; %This is the order for TABLE 1, TABLE 2 and GFF file
chrGffToBed = [1:4,7:10,5,11:17,6];% mapping to the 1:17 order of the bedfile
chrWigToBed = [1, flip(2:17)]; % ONLY WIGFILE (readData). Is this mapping correct? should be

Table1 = table;
Table2 = table;
coordinatesGenes = [];
chromosomePerGene = [];
chromosomeEndPos = [];
geneID = [];
standardGeneID = [];
geneSenseBinary = [];

chrEntries = find(isnan(readDataRaw(:,1)));
chrEntries = [chrEntries; length(readDataRaw)];

for ii = 1:17 % Quick and dirty loop to map the chromosomal orders used in this data to the correct order
    Table1 = [Table1; Table1Raw((chromosomePerGeneRaw == ii),:)];
    Table2 = [Table2; Table2Raw((chromosomePerGeneRaw == ii),:)];
    coordinatesGenes = [coordinatesGenes; coordinatesGenesRaw((chromosomePerGeneRaw == ii),:)];
    chromosomePerGene = [chromosomePerGene; chromosomePerGeneRaw((chromosomePerGeneRaw == ii),:)];
    geneID = [geneID; geneIDRaw((chromosomePerGeneRaw == ii),:)];
    standardGeneID = [standardGeneID; standardGeneIDRaw((chromosomePerGeneRaw == ii),:)];
    geneSenseBinary = [geneSenseBinary; geneSenseBinaryRaw((chromosomePerGeneRaw == ii),:)];

    chromosomeEndPos = [chromosomeEndPos; chromosomeEndPosRaw(chrGffToBed(ii))];
end

% Quick and dirty loop to map the chromosomal orders used in this data to the correct order
readData = readDataRaw(chrEntries(1):chrEntries(2)-1,:); % part 1
readDataRaw = [readDataRaw; [80666, 3]]; % Add a false entry at the end to allow for easy looping. This entry is lost in the mapping loop so no worries
for ii = 1:16
    readData = [readData; readDataRaw((chrEntries(18-ii):chrEntries(19-ii)-1),:)];
end

clear Table1Raw Table2Raw coordinatesGenesRaw chromosomePerGeneRaw geneIDRaw standardGeneIDRaw geneSenseBinaryRaw chromosomeEndPosRaw readDataRaw;
%%
geneCount = height(Table1);
essentialGeneList = zeros(length(standardGeneID),1);
for ii = 1:height(standardGeneID)
    essentialGeneList(ii) = sum(strcmp(standardGeneID{ii},essentialGenes.AllEssentialGenesFoundInLists___Cerevisiae_EssentialGenes_List_)); %For each gene check if it is listed as essential
end

%% put all features and coordinates on one single huge concatenated chromosome
% Why do we have non-unique entries to tn locations? probably they go in reverse directions?

lengthPreviousChr = [0; cumsum(chromosomeEndPos)]; %get the length in bp of each chromosome. The first chromosome starts at 0
lengthPreviousChr = double(lengthPreviousChr(1:end-1)); % Remove the last entry as it does not act as a starting point to another chromosome
%force type double to allow for addition with other double type data. May
%be useful to force double immediately upon reading in, also for other data
tnCoordinates = data;

[tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat] = concatCoordinatesIntoSingleChr(Table1,tnCoordinates,geneCount,lengthPreviousChr);


%% reading in wig file
% NOTE: WIG FILE ORDER OF CHR DOES NOT MATCH OTHER FILES!!!

%Recalculate these for the reorderd data
chrEntries = find(isnan(readData(:,1)));
chrEntries = [chrEntries; length(readData)];

readDataChrCoords = readData;

for ii = 1:17 %ARE THEY IN THE SAME ORDER? OTHERWISE THIS DOES NOT MAKE SENSE. Readdata counts up 1133 - 1148. Mito last
    readDataChrCoords(chrEntries(ii):chrEntries(ii+1),3) = ii; % Add the chr. the coord belongs to
    readData(chrEntries(ii):chrEntries(ii+1),1) = readData(chrEntries(ii):chrEntries(ii+1),1) + lengthPreviousChr(ii); % Turn readdata from chr coordinates to genome coordinates
end

readData(isnan(readData(:,1)),:) = []; %Remove rows which contain a nan value
readDataChrCoords(isnan(readDataChrCoords(:,1)),:) = [];

%% Get number of reads & transposon per gene
tnPerGene = Table2.NumberOfTransposonsPerGene;
readPerGene = Table2.NumberOfReadsPerGene;

geneLength = Table1.EndLocation - Table1.StartLocation; %CAN USE GFF DATA FOR THIS
tnDensity = tnPerGene./geneLength; %define transposon density per gene

% Count number of transposon per gene minus end and beginning
[tnPerGeneMinTen,readPerGeneMinTen,tnDensityMinTen] = getTransposonsPerGeneMinTen(Table1,geneLength,geneCount); %the read count from table1 does not correspond with read count from table2. THIS IS INDEED INCONSISTENT????

% Find the longest transposon free interval within each gene
tnFreeInterval = getTnFreeIntervalPerGene(tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat,geneCount);

%% intergenic regions of 20kb 
% RUNS, BUT DOES IT actuaLLY WORK AS INTENDED? I think so, though I see at the start and end there are 2 equal entires for intergeneic length per gene, so maybe not... 

[intergenicLength20kbPerGene,intergenicTn20kbPerGene,intergenicTnDensity20kb] = getIntergenicTnDensity(geneLength,chromosomeEndPos,tnCoordinates,coordinatesGenes,chromosomePerGene);
[intergenicRead20kbPerGene, intergenicReadDensity20kb] = getIntergenicReadDensity(geneLength,chromosomeEndPos,readDataChrCoords,coordinatesGenes,chromosomePerGene);


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

%% 30/05/22 - Normalization
% Normalization of reads per gene:
% read count * 1/gene length * (10^6/total mapped reads entire genome) *
% (average read density chromosome / average read density local window)
readTotalChr = zeros(17,1); %17 = number of chromosomes
readsPerGeneNorm = zeros(height(Table1),1);

for ii = 1:17
    readTotalChr(ii) = sum(readDataChrCoords(readDataChrCoords(:,3) == ii,2)); % read count per chromosome
end
readDensityChr = readTotalChr./double(chromosomeEndPos); % reads per chr / length of chr
readCountTotal = sum(readData(:,2)); % Total number of reads in the genome
intergenicReadDensity20kAverage = (intergenicReadDensity20kb(:,1) + intergenicReadDensity20kb(:,2))/2;

for ii = 1:height(Table1)
    readsPerGeneNorm(ii) = (readPerGene(ii)/geneLength(ii))*(10^6/readCountTotal).*(readDensityChr(chromosomePerGene(ii))/(1+intergenicReadDensity20kAverage(ii))); % Added + 1 to intergenic read density to prevent inf values %removed *10^6 because it seems silly and has no reason?
end


%% Combine training/testing data into a single table
readnumb = 1; %TEMPORARY SO I CAN CONTINUE. VERY VERY WRONG! WHAT IS READNUMB? %READNUMB IS PROBABLY TOTAL READ COUNT
% 

readPerGenePerTn = readPerGene./(tnPerGene.*sum(readnumb)); %unused? % MAY WANT TO SWITCH TO NORMALIZED READS? BUT TN COUNT IS NOT NORMALIZED...
readPerGenePerTnMinTen = readPerGeneMinTen./(tnPerGeneMinTen.*sum(readnumb)); %unused?
tnFreeIntervalPerGeneLength = (double(tnFreeInterval)./geneLength);
intergenicTn20kbUp = intergenicTn20kbPerGene(:,1);
intergenicTn20kbDown = intergenicTn20kbPerGene(:,2);
NI20kb = tnDensity./((1+intergenicTnDensity20kb(:,1)+intergenicTnDensity20kb(:,2))/2); %IntergenicTnDensity20kb can be 0... --> inf values --> added 1
%NI20kb is compares the density of transposons in the local region to the3
%transposons in the gene
readsPerTnAverage = readPerGenePerTn;

featuresTable = table(readsPerGeneNorm, readsPerTnAverage, readPerGenePerTnMinTen, tnPerGeneMinTen, readPerGeneMinTen, geneLength, tnDensity, tnFreeIntervalPerGeneLength, ...
                        promotorTn, intergenicTn20kbUp, intergenicTn20kbDown, NI20kb, essentialGeneList);% 

% dataTableOfFeatures_WT_ab = dataTableOfFeatures;

%% Train classifier 
% [A, B, C] = trainClassifier(dataTableOfFeatures_WT) %USE BAGGED TREES !!!!
%%
% yfit = trainedModelWTStripped.predictFcn(dataTableOfFeatures);

% classifier = baggedTreesModel_220125;

%%
% imp = oobPermutedPredictorImportance() % FUCNTION TO DETERMINE FEATURE IMPORTANCE. USE: on predictor.ClassificationEnsemble
% labels = trainedModel_220531.ClassificationEnsemble.PredictorNames;

featureImportance = imp/sum(imp);

%% FOR TRAINING REGRESSION PREDICTOR
fitnessData = readingFitnessData(standardGeneID); % Read in experimental fitness data for a set of genes
% fitnessData(essentialGeneList==1) = 0; %Set essential genes to have 0 fitness % Uncomment when using essential genes as well
featuresTable.('readsPerTnMinTen') = readPerGenePerTnMinTen;
featuresTable.('fitness') = fitnessData;
% featuresTable = removevars(featuresTable,{'essentialGeneList'}); % remove essential gene list as we do not need it here
selection = (fitnessData>0); % Selection of all genes that have fitness data
% selection = (fitnessData+essentialGeneList)>0; % Selection of all genes that are essential or have fitness data % Uncomment when using essential genes as well
trainingTable = featuresTable(selection,:);
% trainingTable = removevars(trainingTable,{'essentialGeneList'}); %Remove this variable as it completely maps all the 0 fitness values but we don't have properly for other datasets)

%% Insertion count of gene vs gene length

scatter(geneLength,tnPerGene);
[rho,pval] = corr(geneLength,tnPerGene); %rho of ~0.56 shows a reasonably strong positive correlation

%% Getting read count per insertion according to TABLE 1. TABLE 1 doest not agree with other data afaik. Why?
readsPerInsertionConcat = [];
for ii = 1:6600
    readsPerInsertionConcat = [readsPerInsertionConcat str2num(Table1.ReadsPerInsertionLocation{ii})]; % is it faster to add first and perform str2num only once? probably
end
%% plot read count per insertion
maxk(readsPerInsertionConcat,100);
edges = linspace(0,100,101);

figure(1);
histogram(readsPerInsertionConcat,edges)
title('Reads per insertion in Wild-Type')
xlabel('Reads per insertion')
ylabel('Count')
set(gca,fontsize=20);
xlim([0,100]);

figure(2);
histogram(readsPerTnAverage,edges)
title('Reads per insertion averaged by gene in Wild-Type')
xlabel('Average reads per insertion per gene')
ylabel('Count')
set(gca,fontsize=20);
xlim([0,100]);

%% determine and plot the inter transposon distance

transposonDistance = tnCoordinatesConcat(2:end)- tnCoordinatesConcat(1:end-1);

figure(1);
edges = linspace(0,1000,1001);
histogram(transposonDistance,edges)
title('Insertion interval in Wild-Type')
xlabel('Insertion interval (bp)')
ylabel('Count')
set(gca,fontsize=20);
xlim([0,100]);

figure(2);
N = histcounts(transposonDistance,edges);
scalefactor = length(transposonDistance);

[lambdahat,lambdaci] = poissfit(transposonDistance);


plot(edges,poisspdf(edges,lambdahat)); % Basically the poisson distribution is a poor representation of the data. The tail results in shift to the right. I think it would be best to experimentally determine the 99% confidence interval of finding at least 1 transposon, which I think I can find from the histogram

%% testing
PD = fitdist(N','Poisson');
centerOfEdges = (edges(2:end)+edges(1:end-1))/2;
stem(centerOfEdges, pdf(PD,centerOfEdges)*numel(N));

%% experimental determination of resolution
cumulativePDF = cumsum(N)/sum(N);
confidence = 0.99;
bpConfidence = find(cumulativePDF(cumulativePDF <= confidence),1,'last'); % find the bp value for which the confidence value is reached. IS THIS REALLY THE RIGHT WAY TO DO THIS?


%%
% fitDist = fitdist(transposonDistance,'Poisson'); % lambda = 22.8234 for WT
fitDist = fitdist(transposonDistance); % lambda = 22.8234 for WT


