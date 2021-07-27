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
    'EX_cpd00067','EX_cpd00013'};

biologData = readtable('Data/BiologData.csv');

model = changeRxnBounds(model,cpdsUMS,-1000,'l');
matched = string.empty();
non_matched = string.empty();

for i = 1:size(biologData,1)
    if any(strcmp(model.mets,strcat(biologData{i,2},'[c0]'))) || any(strcmp(model.mets,strcat(biologData{i,2},'[e0]')))
        matched = [matched; biologData{i,2}];
        if any(strcmp(model.rxns,strcat('EX_',char(biologData{i,2}))))
            excRxn = strcat('EX_',char(biologData{i,2}));
            model = changeRxnBounds(model,excRxn,-1000,'l');
            FBA = optimizeCbModel(model,'max','one');
            if FBA.f > 0
                biologData{i,4} = 1;
            else
                biologData{i,4} = 0;
            end
            model = changeRxnBounds(model,excRxn,0,'l');
            clear FBA
        else
            model = addSinkReactions(model,strcat(biologData{i,2},'[e0]'));
            model = changeRxnBounds(model,strcat('sink_',biologData{i,2},'[e0]'),-1000,'l');
            model = changeRxnBounds(model,strcat('sink_',biologData{i,2},'[e0]'),0,'u');
            FBA = optimizeCbModel(model,'max','one');
            if FBA.f > 0
                biologData{i,4} = 1;
            else
                biologData{i,4} = 0;
            end
            clear FBA
            model = removeRxns(model,strcat('sink_',biologData{i,2},'[e0]'));
        end
        
    else
        non_matched = [non_matched; biologData{i,2}];
    end    
end

totalMatched = size(matched,1);
TP = 0;
FP = 0;
TN = 0;
FN = 0;

for k = 1:size(biologData,1)
    if biologData{k,3} == 1
        if biologData{k,4} == 1
            TP = TP + 1;
        else
            FN = FN + 1;
        end
    elseif biologData{k,3} == 0
        if biologData{k,4} == 0
            TN = TN + 1;   
        elseif biologData{k,4} == 1  
            FP = FP + 1;
        end
    end
end

sprintf('True positives: %d\nTrue negatives: %d\nFalse positives: %d\nFalse negatives: %d',...
    TP,TN,FP,FN)

accuracy = (TP+TN)/totalMatched;
precision = TP/(TP+FP);
recall = TP/(TP+FN);

sprintf('Accuracy: %.4f\nPrecision: %.4f\nRecall: %.4f',accuracy,precision,recall)
