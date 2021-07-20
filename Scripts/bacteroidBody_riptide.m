clear
clc

model = readCbModel('iCS1224.mat');
model = changeObjective(model,'rxn06874');

%% Define bacteroid boundary conditions

[selExc, selUpt] = findExcRxns(model,0,0);
excRxns = model.rxns(selExc);
model = changeRxnBounds(model, excRxns, 0, 'l');

%add demand 
model = addDemandReaction(model,'cpd11463[c0]');

diffRxns = {'rxn05202','rxn05469','rxn05488','rxn05561','rxn05634','rxn05682',...
    'rxn08423','rxn08525','rxn08651','rxn08863','rxn09256','rxn09269','rxn10153',...
    'rxn12569','rxnTDIFFcpd00441','rxnTPScpd01217','rxn12671','rxn12591'};

model = changeRxnBounds(model,diffRxns,0,'l');

% Add demand reactions for protein and lipids
model = changeRxnBounds(model,'DM_cpd11463[c0]',1,'l');
model = changeRxnBounds(model,'DM_cpd01080',1,'l');
model = changeRxnBounds(model,'DM_cpd25615',1,'l');

%prevent stoichiometrically balanced cycles through these reactions
model = changeRxnBounds(model,'rxn06033',0,'b');
model = changeRxnBounds(model,'rxn00669',0,'b');
model = changeRxnBounds(model,'rxn00868',0,'b');
model = changeRxnBounds(model,'rxn00872',0,'b');

[names, IDs] = textread('Data/noduleCpds_biosensor.txt','%s%s');

for i = 1:size(IDs,1)
    if any(strcmp(model.rxns,char(strcat('EX_',IDs(i)))))
        excRxn = strcat('EX_',IDs(i));
        model = changeRxnBounds(model,excRxn,-1000,'l');
    else
        model = addSinkReactions(model,strcat(IDs(i),'[c0]'));
        model = changeRxnBounds(model,strcat('sink_',IDs(i),'[c0]'),-1000,'l');
        model = changeRxnBounds(model,strcat('sink_',IDs(i),'[c0]'),0,'u');
    end
end

Nsecretioncpds = {'cpd00013','cpd00035','cpd00041'};

for i = 1:numel(Nsecretioncpds)
    excRxn = strcat('EX_',Nsecretioncpds(i));
    model = changeRxnBounds(model,excRxn,0,'l');
    model = changeRxnBounds(model,excRxn,1000,'u');
end

%set oxygen limitation if required
%model = changeRxnBounds(model,'EX_cpd00007',-100,'l');

%allow glycolate uptake and inactivate isocitrate lyase
model = changeRxnBounds(model,'EX_cpd00139',-1000,'l');
model = changeRxnBounds(model,'rxn00336',0,'b');


FBA = optimizeCbModel(model,'max','one')
writeCbModel(model,'bacteroidBody.xml');
