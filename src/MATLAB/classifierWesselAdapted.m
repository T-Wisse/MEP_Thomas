%https://sites.google.com/site/satayusers/

clear all 
close all

%ORDER OF FILES LOADED: BEDFILE, perInsertionsFile, perGeneFile, essentialGeneFile

%% Load insertion data (.bed file)

[bedfile, path]=uigetfile('Pick a file'); %<-point to the .bed file
cd(path)


%pergene_insertions.txt: Gene name	Chromosome	Start location	End location	Insertion locations	Reads per insertion location
%.bed file: Chromosome name start location end location . score

%We are only interested in the first 2 entries for the .bed file. The end location is
%just start loc + 1

data = load(bedfile, '-ascii'); %Load only works if file has no text. Note: automatically removes zeros (0) at the start of a number
data = data(:,[1,2]); %Select columns 1 and 2. It should be a score??
load('yeastGFF.mat') %<- both yeastGFF.mat should be in the same folder as the BAM file. Otherwise adapt here



%% Load insertion data per gene (.txt file)

% 
% perInsertionsFile = 'yLIC137_7_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt';
% perGeneFile = 'yLIC137_7_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt';
% essentialGeneFile = 'yLIC137_7_merged_cleaned_forward_reads_trimmed.sorted.bam_peressential.txt';

perInsertionsFile = 'yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt';
perGeneFile = 'yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt';
essentialGeneFile = 'yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam_peressential.txt';
wigFile = 'yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam.wig';

%Should be able to do this with path and regexp... but for now we cheat
% 
% perInsertionsFile = 'WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt';
% perGeneFile = 'WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt';
% essentialGeneFile = 'WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_peressential.txt';

% [perInsertionsFile, path]=uigetfile('Pick a file'); %<-point to the pergene_insertions.txt file
T = readtable(perInsertionsFile); %alternative to load, works with text.

% [perGeneFile, ~]=uigetfile('Pick a file'); %<-point to the pergene_insertions.txt file
Table2 = readtable(perGeneFile); %alternative to load, works with text.

% [essentialGeneFile, path]=uigetfile('Pick a file'); %<-point to the peressential.txt file
Table3 = readtable(essentialGeneFile); %alternative to load, works with text.
essentialGeneList = Table3.GeneName;

geneCount = height(T);


%% transformation of used chromosome identifiers to 1:17. Works but it's an ugly fix

chromosomeIdentifiers = unique(data(:,1)); %Get the 17 unique chromosome identifiers used in the data
chromosomeIdentifiers = chromosomeIdentifiers - 1132;
data(:,1) = data(:,1)-1132;
data(data(:,1)==92,1) = data((data(:,1)==92),1)-75;
if unique(data(:,1)) ~= [1:17]
    disp('An issue occured with transformation of chromosome identifiers')
end




%% NOTE: RUNS. BUT I DON't USE THIS FILE SINCE THE RUN WAS DONE WITH A DIFFERENT FILE...?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%prepare chromosomal features from gff file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


chromosomeEndPos = cell2mat(gff(find(strcmp(gff(:,3),'omosome')),5)); % find chromosome lengths


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create array for all bp positions per chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% res = 200;
% 
% 
% 
% for k=1:length(chromosomeEndPos)
%     my_field = strcat('ch',num2str(k));     % create variable with all bp positions for every chromosome
%     tnpos.(my_field) = zeros(chromosomeEndPos(k),1);
%     readpos.(my_field) = zeros(chromosomeEndPos(k),1);
%     
%     x=tnCoordinates(:,1)==k; %find all positions in tncoordinates that give tn locations for specific (k) chromosome
%     y=tnCoordinates(x,2); %takes coordinates of tns on chromosome k
%     tnpos.(my_field)(y) = 1; % gives a 1 for every bp position that has a transposon insertion
%     readpos.(my_field)(y) = readnumb(x); % puts number of reads on every bp position
%     
% end

 
%% Probably works now. Why do we have non-unique entries to tn locations? maybe they go in reverse directions?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%put all features and  coordinate on one single huge concatenated chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lengthPreviousChr = [0; cumsum(chromosomeEndPos)]; %get the length in bp of each chromosome. The first chromosome starts at 0
lengthPreviousChr = lengthPreviousChr(1:end-1); % Remove the last entry as it does not act as a starting point to another chromosome

tnCoordinates = data;
tnCoordinatesConcat = [];

geneChromosomeLoc = zeros(geneCount,1);

for ii = 1:geneCount
    geneChromosomeLoc(ii) = roman2num(T.Chromosome{ii});
end
geneChromosomeLoc(isnan(geneChromosomeLoc)) = 17; %Mitochondrial DNA (denoted as 'mito') is returned as NaN by roman2num. All this DNA is counted as chromosome 17
geneStartCoordinatesConcat = T.StartLocation;
geneEndCoordinatesConcat = T.EndLocation;

for ii=1:17
    %NOTE: tnCoordinatesConcat includes double locations. Presumably where
    %transposons hit either direction
    tnCoordinatesConcat= [tnCoordinatesConcat; tnCoordinates(tnCoordinates(:,1)==ii,2)+lengthPreviousChr(ii)]; %Translate all transposon coordinates to a total bp. 
    geneStartCoordinatesConcat(geneChromosomeLoc==ii) = geneStartCoordinatesConcat(geneChromosomeLoc==ii)+lengthPreviousChr(ii); %Translate all gene start coordinates to a total bp. 
    geneEndCoordinatesConcat(geneChromosomeLoc==ii) = geneEndCoordinatesConcat(geneChromosomeLoc==ii)+lengthPreviousChr(ii); %Translate all gene end coordinates to a total bp. 
end

%% reading in wig file

filename = 'D:\Thomas\Documents\Uni\MEP_Data\SATAY_Data\SATAY_Analysis\yTW001_4\yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam.wig';
readData = loadWigfile(filename); %CHANGE LOCATION! Location currently at D:\Thomas\Documents\Uni\MEP_Data\SATAY_Data\SATAY_Analysis\yTW001_4
chrEntries = find(isnan(readData(:,1)));
chrEntries = [chrEntries; length(readData)];

for ii = 1:17 %ARE THEY IN THE SAME ORDER? OTHERWISE THIS DOES NOT MAKE SENSE. Readdata counts up 1133 - 1148. Mito last
    readData(chrEntries(ii):chrEntries(ii+1),1) = readData(chrEntries(ii):chrEntries(ii+1),1) + lengthPreviousChr(ii); 
end
readData(isnan(readData(:,1)),:) = [];

%% Testing: Readpos, tnpos
% 
% readPosConcat = zeros(sum(chromosomeEndPos),1); %Array with the length of all chromosome combined
% 
% tnPosConcat = zeros(sum(chromosomeEndPos),1); %Array with the length of all chromosome combined
% tnPosConcat(tnCoordinatesConcat) = 1; %Set locations of transposons in the geneome to 1
% 
% readCountConcat = [];
% for ii = 1:geneCount
%     readCountConcat = [readCountConcat str2num(T.ReadsPerInsertionLocation{ii})];
% end
% 
% readPosConcat(logical(tnPosConcat)) =  readCountConcat';

%% brute forcing readpos?

% 
% kk = 0;
% for ii = 1:geneCount
%     myField = T.Chromosome{ii}; 
%     if ~exist(strcat('readCount.(',myField,')'))
%         readCount.(myField) = [];
%     end
%     for jj = 1:length(str2num(A{ii}))
%         kk = kk+1;
%         readCount.(myField) = [readCount.(myField) str2num(T.ReadsPerInsertionLocation{ii})];
%     end
% end
%% Get number of reads & transposon per gene
tnPerGene = Table2.NumberOfTransposonsPerGene;
readPerGene = Table2.NumberOfReadsPerGene;

geneLength = T.EndLocation - T.StartLocation;
% geneLength = genes.coordinates(:,2)-genes.coordinates(:,1); %define gene length USES a DIFFERENT LIST OF GENES THAN READ DATA HAS!
tnDensity = tnPerGene./geneLength; %define transposon density per gene


%% Count number of transposon per gene minus end and beginning


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

%% intergenic regions of 20kb % TO DO

wl = 20000; % define window length

roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "mt"];
t=0;
u=0;
intergeniclength20kbpergene = zeros(length(geneLength),2);


for k=1:length(chromosomeEndPos)
    my_field = strcat('ch',num2str(k));     % create variable with all bp positions for every chromosome
    genepos.(my_field) = zeros(chromosomeEndPos(k),1);
    
    for i=t+1:t+sum(strcmp(roman(k),genes.chr)) % create binary with 1 for every bp that is part of a gene for every chromosome
    genepos.(my_field)(genes.coordinates(i,1):genes.coordinates(i,2)) = 1; 
    t = i;
    end
    
    intergenictnpos.(my_field)=tnpos.(my_field);
    genepos.(my_field)=logical(genepos.(my_field));
    intergenictnpos.(my_field)(genepos.(my_field)) = 0;
    
    for j=u+1:u+sum(strcmp(roman(k),genes.chr)) % every gene of each chromosome
        if genes.coordinates(j,1)-wl<1
            intergeniclength20kbpergene(j,1) = sum(1-(genepos.(my_field)(1:genes.coordinates(j,1)))); % calculate intergenic length before gen in 20kb window
            intergenictn20kbpergene(j,1) = sum(intergenictnpos.(my_field)(1:genes.coordinates(j,1))); % calculates number of tn in this region
        else
            intergeniclength20kbpergene(j,1) = sum(1-(genepos.(my_field)(genes.coordinates(j,1)-wl:genes.coordinates(j,1)))); % calculate intergenic length before gen in 20kb window
            intergenictn20kbpergene(j,1) = sum(intergenictnpos.(my_field)(genes.coordinates(j,1)-wl:genes.coordinates(j,1)));
        end
        
        if genes.coordinates(j,2)+wl>chromosomeEndPos(k)
             intergeniclength20kbpergene(j,2) = sum(1-(genepos.(my_field)((genes.coordinates(j,2):chromosomeEndPos(k))))); % calculate intergenic length before gen in 20kb window
             intergenictn20kbpergene(j,2) = sum(intergenictnpos.(my_field)(genes.coordinates(j,2):chromosomeEndPos(k)));
        else
            intergeniclength20kbpergene(j,2) = sum(1-(genepos.(my_field)(genes.coordinates(j,2):genes.coordinates(j,2)+wl))); % calculate intergenic length before gen in 20kb window
            intergenictn20kbpergene(j,2) = sum(intergenictnpos.(my_field)(genes.coordinates(j,2):genes.coordinates(j,2)+wl));
        end
        u=j;
    end
    
end

intergenic_tndensity20kb = intergenictn20kbpergene./intergeniclength20kbpergene;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Longest transposon free interval within a gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:geneCount %Loop trough each gene. For each, find the longest interval between transposons in bp
    ww = tnCoordinatesConcat>=geneStartCoordinatesConcat(ii)&tnCoordinatesConcat<=geneEndCoordinatesConcat(ii);
    ww = tnCoordinatesConcat(ww);
    ww = [geneStartCoordinatesConcat(ii) ; ww ; geneEndCoordinatesConcat(ii)]; %Concat the transposon locations of gene ii with the start and end of the gene
    tnFreeInterval(ii,1) = max(diff(ww)); %determine the largest difference bertween previously concatenated locations to find longest region without a transposon in the gene 
end


%% creating a list with the number of transposons and reads in the 100bp 5' of each gene (promotor region)
% % TO DO ... This I need to use the genes file for...
% 
% genesense = cell2mat(gff(features.genes,7));
% genesensebinary = genesense=='+';
% 
% % creating a concatinated verson of all transposons and reads.
% tnpos_concat = [tnpos.ch1 ; tnpos.ch2 ; tnpos.ch3 ; tnpos.ch4 ; tnpos.ch5 ; tnpos.ch6 ; tnpos.ch7 ; tnpos.ch8 ; tnpos.ch9 ; tnpos.ch10 ; tnpos.ch11 ; tnpos.ch12 ; tnpos.ch13 ; tnpos.ch14 ; tnpos.ch15 ; tnpos.ch16 ; tnpos.ch17];
% readpos_concat = [readpos.ch1 ; readpos.ch2 ; readpos.ch3 ; readpos.ch4 ; readpos.ch5 ; readpos.ch6 ; readpos.ch7 ; readpos.ch8 ; readpos.ch9 ; readpos.ch10 ; readpos.ch11 ; readpos.ch12 ; readpos.ch13 ; readpos.ch14 ; readpos.ch15 ; readpos.ch16 ; readpos.ch17];
% 
% 
% for k=1:length(genes.coordinates_concat)
% 
% if genesensebinary(k)==1    
%     promotortn(k)=sum(tnpos_concat(genes.coordinates_concat(k,1)-100:genes.coordinates_concat(k,1)));
%     promotorread(k)=sum(readpos_concat(genes.coordinates_concat(k,1)-100:genes.coordinates_concat(k,1)));
% else
%     promotortn(k)=sum(tnpos_concat(genes.coordinates_concat(k,2):genes.coordinates_concat(k,2)+100));
%     promotorread(k)=sum(readpos_concat(genes.coordinates_concat(k,2):genes.coordinates_concat(k,2)+100));
% end
% end
% 
% promotorreadpertn=promotorread./promotortn;



%% Determine which genes are essential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Statistical Learning Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
essentiality = zeros(geneCount,1);

for ii = 1:geneCount
    
    if sum(strcmp(essentialGeneList, T.GeneName{ii})) == 1
        essentiality(ii) = 1;
    end
end

%% form table

tnpergeneT = tnPerGene';
readpergeneT = readPerGene';
tndensityT = tnDensity;

% essentiality = char(essentiality);
essentialyesno=essentiality;

readnumb = 1 %TEMPORARY SO I CAN CONTINUE. VERY VERY WRONG! %READNUMB IS PROBABLY THE READ COUNT PER BP

readpergenepertn = readPerGene./(tnPerGene.*sum(readnumb));
readpergenepertnminten = readPerGeneMinTen./(tnPerGeneMinTen.*sum(readnumb));
tnFreeIntervalPerGeneLength = (double(tnFreeInterval)./geneLength);
% intergenictn20kbup = intergenictn20kbpergene(:,1);
% intergenictn20kbdown = intergenictn20kbpergene(:,2);
% NI20kb = tndensityT./((intergenic_tndensity20kb(:,1)+intergenic_tndensity20kb(:,2))/2);

% Select which data set you are working with

%EssentialityTable_WTI = table(tnpergeneminten', readpergeneminten', genelength, tndensityminten', tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb, essentialyesno);
%EssentialityTable_WTII = table(tnpergeneminten', readpergeneminten', genelength, tndensityminten', tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb, essentialyesno);

% TestTable_DplI = table(tnPerGeneMinTen', readPerGeneMinTen', geneLength, tndensityT, tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb);% 

%TestTable_DplIPsdII = table(tnpergeneminten', readpergeneminten', genelength, tndensityT, tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb);% 

TestTable_DplI = table(tnPerGeneMinTen, readPerGeneMinTen, geneLength, tndensityT, tnFreeIntervalPerGeneLength, essentialyesno);% 


%% Normalization for selected datasets to combine them into a single set

%tntot1=410169; %WT 1
%readtot1=31794831; %WT 1
%tntot2=356875; %WT 2
%readtot2=1530285; %WT 2

% tnCountA = 100000; %Total number of transposons found in the experiment A
% readCountA = 1000000; %Total number of reads found in the experiment A
% 
% tableNormA = TestTable_DplI; %Select table A for normalization
% tableNormA{:,[1 4:5]} = tableNormA{:,[1, 4:5]}/tnCountA; %divide all entries based on transposon count (name includes tn) by the total number of transposons
% tableNormA{:,2} = tableNormA{:,2}/readCountA; %divide all entries based on read count by the total number of reads
% 
% 
% tnCountB = 100000; %Total number of transposons found in experiment B
% readCountB = 1000000; %Total number of reads found in experiment B
% 
% tableNormB = TestTable_DplI; %Select table B for normalization
% tableNormB{:,[1 4:5]} = tableNormB{:,[1, 4:5]}/tnCountB; %divide all entries based on transposon count (name includes tn) by the total number of transposons
% tableNormB{:,2} = tableNormB{:,2}/readCountB; %divide all entries based on read count by the total number of reads
% 
% EssentialityTable_A_B = [tableNormA ; tableNormB];

% Use Apps:  Classification Learner to train the classifier
% 
% %% Create training data set with equal number of ess and non-ess genes
% 
% A1=rand(6603,1)<=0.82;
% A2=logical(A1-essentialyesno);
% essentialyesnoreduced=essentialyesno;
% essentialyesnoreduced(A2)=[];
% tnpergenemintenred=tnPerGeneMinTen;
% tnpergenemintenred(A2)=[];
% readpergenemintenred=readPerGeneMinTen;
% readpergenemintenred(A2)=[];
% genelengthred=geneLength;
% genelengthred(A2)=[];
% tndensityred=tnDensity;
% tndensityred(A2)=[];
% tndensitymintenred=tnDensityMinTen;
% tndensitymintenred(A2)=[];
% %tndensityupred=tndensityup;
% %tndensityupred(A1)=[];
% %tndensitydownred=tndensitydown;
% %tndensitydownred(A1)=[];
% tnfreeintervalpergenelengthred=tnFreeIntervalPerGeneLength;
% tnfreeintervalpergenelengthred(A2)=[];
% promotortnred=promotortn;
% promotortnred(A2)=[];
% intergenictn20kbupred=intergenictn20kbup;
% intergenictn20kbupred(A2)=[];
% intergenictn20kbdownred=intergenictn20kbdown;
% intergenictn20kbdownred(A2)=[];
% NI20kbred=NI20kb;
% NI20kbred(A2)=[];
% 
% 
% %EssentialityTableGred = table(tnpergenemintenred', readpergenemintenred', genelengthred, tndensityred', tndensityupred, tndensitydownred, tnfreeintervalpergenelengthred', essentialyesnoreduced); %
% TestTable_DplI_red = table(tnpergenemintenred', readpergenemintenred', genelengthred, tndensityred', tnfreeintervalpergenelengthred', promotortnred', intergenictn20kbupred, intergenictn20kbdownred, NI20kbred);%
% %EssentialityTable_WTI_red = table(tnpergenemintenred', readpergenemintenred', genelengthred, tndensitymintenred', tnfreeintervalpergenelengthred', promotortnred', intergenictn20kbupred, intergenictn20kbdownred, NI20kbred, essentialyesnoreduced);
% 
% 
% %% Self prediction
% 
%  Predicted_WT2_WT1WT2training = SATAY_classifier_WT1_WT2_bagged.predictFcn(EssentialityTableK_WT2);
%  
%  Predicted_WT1_WT2training_bagged_10learners = SATAY_classifier_WTII_bagged_10learners.predictFcn(EssentialityTableK_WT1_norm);
%  sum(Predicted_WT1_WT2training_bagged_10learners==essentialyesno)
% 
% %% Test DplI data set
% 
% tntotDpl=590079;
% readtotDpl=15077158;
% 
% TestTable_DplI_norm = TestTable_DplI;
% 
% TestTable_DplI_norm{:,1} = TestTable_DplI{:,1}/tntotDpl;
% TestTable_DplI_norm{:,2} = TestTable_DplI{:,2}/readtotDpl;
% TestTable_DplI_norm{:,4:9} = TestTable_DplI{:,4:9}/tntotDpl;
% 
% Predicted_DplI_WT1WT2training = SATAY_classifier_WT1_WT2_bagged.predictFcn(TestTable_DplI_norm);
% %essdifference_DplI_WT = (Predicted_WT1_WT1WT2training == Predicted_DplI_WT1WT2training);
% 
% %% and DplIPsdII data set
% 
% tntotDplPsd=531369;
% readtotDplPsd=11649561;
% 
% TestTable_DplIPsdII_norm = TestTable_DplIPsdII;
% 
% TestTable_DplIPsdII_norm{:,1} = TestTable_DplIPsdII{:,1}/tntotDplPsd;
% TestTable_DplIPsdII_norm{:,2} = TestTable_DplIPsdII{:,2}/readtotDplPsd;
% TestTable_DplIPsdII_norm{:,4:9} = TestTable_DplIPsdII{:,4:9}/tntotDplPsd;
% 
% Predicted_DplIPsdII_WT1WT2training = SATAY_classifier_WT1_WT2_bagged.predictFcn(TestTable_DplIPsdII_norm);
% essdifference_DplIPsdII_WT = (Predicted_WT1_WT1WT2training == Predicted_DplIPsdII_WT1WT2training);
% 
% %% classification score
% 
% kNNMdl = SATAY_classifier_WT1_WT2_bagged.ClassificationEnsemble;
% [labels,score] = predict(kNNMdl,TestTable_DplI_norm);
% Class_score_DplI_WT1WT2 = score;
% 
% [labels,score] = predict(kNNMdl,TestTable_DplIPsdII_norm);
% Class_score_DplIPsdII_WT1WT2 = score;
% 
% %% check different essential genes in DplI and DplIPsdII with >95% accuracy
% 
% %sum(Predicted_DplIPsdII_WT1WT2training==Predicted_DplI_WT1WT2training)
% %Predicted_ess_in_DplI_not_in_DplIPsdII = find((Predicted_DplIPsdII_WT1WT2training-Predicted_DplI_WT1WT2training)==-1);
% 
% Predicted_Essential_DplI_95=(Class_score_DplI_WT1WT2(:,2)>=0.95);
% Predicted_Essential_DplIPsdII_95=(Class_score_DplIPsdII_WT1WT2(:,2)>=0.95);
% 
% Predicted_ess95_in_DplI_not_in_DplIPsdII = find((Predicted_Essential_DplIPsdII_95-Predicted_Essential_DplI_95)==-1);
% Predicted_ess95_in_DplIPsdII_not_in_DplI = find((Predicted_Essential_DplIPsdII_95-Predicted_Essential_DplI_95)==1);
% 
% % Histogram of classification scores
% 
% figure(1)
% hist(Class_score_DplI_WT1WT2(:,1),100)
% xlabel('Classification score')
% ylabel('#genes')
% title('Distribution of Classification scores in DplId')
% set (gca, 'Fontsize', 20)
% 
% %% test
% 
% TestTable_DplI_norm_red = TestTable_DplI_red;
% TestTable_DplI_norm_red{:,1} = TestTable_DplI_red{:,1}/tntotDpl;
% TestTable_DplI_norm_red{:,2} = TestTable_DplI_red{:,2}/readtotDpl;
% TestTable_DplI_norm_red{:,4:9} = TestTable_DplI_red{:,4:9}/tntotDpl;
% %%
% kNNMdl = SATAY_classifier_WT1_norm_bagged30.ClassificationEnsemble;
% [labels,score] = predict(kNNMdl,EssentialityTableK_WT2_norm);
% 
% Class_score_WT2_WT1tested = score;
% 
% figure(1)
% hist(Class_score_WT2_WT1tested(:,1),30)
% xlabel('Classification score')
% ylabel('#genes')
% title('Distribution of Classification scores in DplId')
% set (gca, 'Fontsize', 20)
% 
% %%
% kNNMdl = SATAY_Classifier_WT1_100bagged.ClassificationEnsemble;
% [labels,score] = predict(kNNMdl,EssentialityTableK_WT2_norm);
% 
% Class_score_WT2_WT1tested = score;
