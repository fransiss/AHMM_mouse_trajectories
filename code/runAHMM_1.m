
%cd C:\FullBNT-1.0.4  % modify as needed - where the bnt code is - from https://github.com/bayesnet/bnt
%addpath(genpathKPM(pwd))

%Read in the entire data
addpath(genpath(pwd));

% Read the data from the specified subfolder
table = readtable('data/clean_weight.xlsx'); 

[ids_ok, ids_Animal] = findgroups(table.AnimalID);
save('results/ids_Animal.mat', 'ids_Animal');  

[genotypes_ok, genotypesID] = findgroups(table.Genotype);
save('results/genotypesID.mat', 'genotypesID');

[sex_ok, sedxID] = findgroups(table.Sex);
save('results/sedxID.mat', 'sedxID');  

dat = [ids_ok, genotypes_ok, sex_ok, table.AgeInMonth, table.Weight];


[m,n]=size(dat);%get number of vars
h=10;%hidden states

MTSids=dat(:,1);
MTS=dat(:,2:n);


%Standardise the time-series
%sMTS=standardize(MTS);
%sMTS(1,:)=MTS(1,:)+1;
%build the cell array and train an ARHMM

[~,~,MTSids_uni]=unique(MTSids);
A = splitapply(@(v) {v}, MTS, MTSids_uni);
PMTS=[];
for i=1:length(A)
     PMTS{i}=num2cell(A{i}');
end
m=max(MTSids);

niters=10;
repetition_cell={};
for iter = 1:niters
    [bnet2, LLtrace] = TrainHMM(PMTS, 2, [], h);
    save(['results/bnet2_', num2str(iter), '.mat'], 'bnet2');  
    save(['results/LLtrace_', num2str(iter), '.mat'], 'LLtrace');  
    repetition_cell{iter, 1} = LLtrace;
    repetition_cell{iter, 2} = bnet2;
    
    % Plot the learning curve
    figure;
    plot(LLtrace(2:length(LLtrace)));
    
    % Explore the hidden state transitions
    get_field(bnet2.CPD{6}, 'cpt');
    
    % Use inference on the training data to infer the hidden states
    out = cell(1, m);
    for i = 1:m
        [mm, nn] = size(PMTS{i});
        % Original data with hidden state sample trajectories
        engine = jtree_unrolled_dbn_inf_engine(bnet2, nn);
        ev = [cell(1, nn); PMTS{i}];
        [engine, ll] = enter_evidence(engine, ev);
        for t = 1:nn
            marg = marginal_nodes(engine, 1, t);
            ev{1, t} = argmax(marg.T);
        end
        out{i} = ev;
    end
    save(['results/out5hidden_correct_', num2str(iter), '.mat'], 'out'); 
end

save('results/repetition_cell.mat', 'repetition_cell');  