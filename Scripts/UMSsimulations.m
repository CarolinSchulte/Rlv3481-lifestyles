clear
clc

model = readCbModel('iCS1224.mat');
model = changeObjective(model,'rxnBIOMASS');

%% Define minimal media boundary conditions

[selExc, selUpt] = findExcRxns(model,0,0);
excRxns = model.rxns(selExc);
model = changeRxnBounds(model, excRxns, 0, 'l');

%list of exchange reactions for compounds present in UMS media
cpdsUMS = {'EX_cpd00305','EX_cpd00104','EX_cpd10516',...
    'EX_cpd00099','EX_cpd00205','EX_cpd00254','EX_cpd00048',...
    'EX_cpd00030','EX_cpd00971','EX_cpd00034','EX_cpd11574','EX_cpd00058',...
    'EX_cpd00149','EX_cpd00644','EX_cpd00009','EX_cpd00063','EX_cpd00007',...
    'EX_cpd00067','EX_cpd00013'};

model = changeRxnBounds(model,cpdsUMS,-1000,'l');

%read C13 flux data from Terpolilli et al. (DOI: 10.1128/JB.00451-16)
fluxData = readtable('Data/C13data.csv');

%succinate as carbon source, set uptake to one
model = changeRxnBounds(model,'EX_cpd00036',-1,'b');
FBA = optimizeCbModel(model,'max','one');

for i = 1:size(fluxData,1)
    flux = FBA.v(findRxnIDs(model,fluxData{i,1}));
    fluxData{i,5} = flux;
end

[minFlux, maxFlux] = fluxVariability(model,95,'max','allowLoops',0);

for i = 1:size(fluxData,1)
    lb = minFlux(findRxnIDs(model,fluxData{i,1}));
    ub = maxFlux(findRxnIDs(model,fluxData{i,1}));
    fluxData{i,6} = lb;
    fluxData{i,7} = ub;
end
