clear
clc

model = readCbModel('iCS1224.mat');
model = changeObjective(model,'rxnBIOMASS');

%% Define minimal media boundary conditions

[selExc, selUpt] = findExcRxns(model,0,0);
excRxns = model.rxns(selExc);
model = changeRxnBounds(model, excRxns, 0, 'l');

%List of exchange reactions for compounds present in UMS media
cpdsUMS = {'EX_cpd00305','EX_cpd00104','EX_cpd10516',...
    'EX_cpd00099','EX_cpd00205','EX_cpd00254','EX_cpd00048',...
    'EX_cpd00030','EX_cpd00971','EX_cpd00034','EX_cpd11574','EX_cpd00058',...
    'EX_cpd00149','EX_cpd00644','EX_cpd00009','EX_cpd00063','EX_cpd00007',...
    'EX_cpd00067'};

succinate = 'EX_cpd00036';
NSource = 'EX_cpd00013';

model = changeRxnBounds(model,cpdsUMS,-1000,'l');
model = changeRxnBounds(model,succinate,-1000,'l');
model = changeRxnBounds(model,NSource,-1000,'l');

%essential genes based on DOI: 10.1128/JB.00572-16 and DOI: 10.1073/pnas.2009094117
esGenes = readtable('Data/essentialGenes.csv','Format','%s%s','Delimiter',',');
esArray = table2array(esGenes);
esArray = esArray(:,1);

[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model,'MOMA');


esPredicted = cell2table(cell(0,2));
TP = string.empty();
TN = string.empty();
FP = string.empty();
FN = string.empty();
growthThreshold = 0.5;

for i = 1:size(model.genes,1)
    esPredicted = [esPredicted; cell2table({model.genes(i),grRatio(i)})];
end


for i = 1:size(model.genes,1)
    if esPredicted{i,2} < growthThreshold || isnan(esPredicted{i,2})
        if any(strcmp(esArray,char(esPredicted{i,1})))
            TP = [TP; char(esPredicted{i,1})];
        else
            FP = [FP; char(esPredicted{i,1})];      
        end
    elseif esPredicted{i,2} > growthThreshold 
        if ~any(strcmp(esArray,char(esPredicted{i,1})))
            TN = [TN; char(esPredicted{i,1})];
        else   
            FN = [FN; char(esPredicted{i,1})];
        end
    end
end

FP(FP == 'Unknown') = [];
FP(FP == 'Spontaneous') = [];

sprintf('True positives: %d\nTrue negatives: %d\nFalse positives: %d\nFalse negatives: %d',...
    size(TP,1),size(TN,1),size(FP,1),size(FN,1))

accuracy = (size(TP,1)+size(TN,1))/(size(TP,1)+size(FP,1)+size(TN,1)+size(FN,1));
precision = size(TP,1)/(size(TP,1)+size(FP,1));
recall = size(TP,1)/(size(TP,1)+size(FN,1));

sprintf('Accuracy: %.4f\nPrecision: %.4f\nRecall: %.4f',accuracy,precision,recall)
