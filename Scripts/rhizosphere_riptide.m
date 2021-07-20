clear
clc

model = readCbModel('iCS1224.mat');
model = changeObjective(model,'rxnBIOMASS');

[selExc, selUpt] = findExcRxns(model,0,0);
excRxns = model.rxns(selExc);
model = changeRxnBounds(model, excRxns, 0, 'l');
model = changeRxnBounds(model,'EX_cpd00013',0,'b');
model = changeRxnBounds(model,'EX_cpd00035',0,'u');
model = changeRxnBounds(model,'EX_cpd00041',0,'u');

%set demand for Nod factor, EPS and LPS
model = changeRxnBounds(model,'DM_cpdNF',1,'l');
model = changeRxnBounds(model,'DM_cpdEPS',1,'l');
model = changeRxnBounds(model,'DM_cpdLPS',1,'l');

%set lower bound of diffusion reactions to 0
diffRxns = {'rxn05202','rxn05469','rxn05488','rxn05561','rxn05634','rxn05682',...
    'rxn08423','rxn08525','rxn08651','rxn08863','rxn09256','rxn09269','rxn10153',...
    'rxn12569','rxnTDIFFcpd00441','rxnTPScpd01217','rxn12671','rxn12591'};

model = changeRxnBounds(model,diffRxns,0,'l');

%[names, IDs] = textread('rhizosphereCpds.txt','%s%s');
[names, IDs] = textread('Data/rhizosphereCpds.txt','%s%s');

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

FBA = optimizeCbModel(model,'max','one')
writeCbModel(model,'rhizo_riptide.xml');        