clear all
close all
niters=10;
nstatesi=10;
nfig=1;


% Load data
table = readtable('data/clean_weight.xlsx');
table_f = table(strcmp(table.Sex, 'Female'), :);
table_m = table(strcmp(table.Sex, 'Male'), :);

real_genotype={};
real_sex={};
all_ids=unique(table.AnimalID);
for u =1:length(all_ids)
    index_ids=find(table.AnimalID==all_ids(u));
    genotype_u=table.Genotype(index_ids);
    sex_u=table.Sex(index_ids);
    real_genotype(end+1) = unique(genotype_u);
    real_sex(end+1) = unique(sex_u);
end 

figure(nfig)
pie(categorical(real_genotype));
nfig=nfig+1;
figure(nfig)
pie(categorical(real_sex));
nfig=nfig+1;

last_trace_val=[];
for iter = 1:niters
    load(strcat('results/trajectories10states/LLtrace_', num2str(iter), '.mat'))
    last_trace_val(iter) = x(end);
end

% select the results with the max trace of the 10 iterations that were run
[mini ind]=max(last_trace_val);

% Load best model output
load(strcat('results/trajectories10states/out5hidden_correct_', num2str(ind), '.mat'))
out = x;
load(strcat('results/trajectories10states/bnet2_', num2str(ind), '.mat'))
bnet2 = x;

trans_matrix=get_field(bnet2.CPD{6},'cpt');
trans_matrix=round(trans_matrix,5);

% check ending trajectories >0.1 probabilities
trans_matrix_masked=trans_matrix.*(trans_matrix>0.1);
trans_matrix_masked=round(trans_matrix_masked,5);

% check ending trajectories >0.3 probabilities
trans_matrix_masked3=trans_matrix.*(trans_matrix>0.3);
trans_matrix_masked3=round(trans_matrix_masked3,5);

figure(nfig)
imagesc(trans_matrix);
colormap(flipud(bone));
colorbar;
axis square
h = gca;
h.XTick = 1:length(trans_matrix);
h.YTick = 1:length(trans_matrix);
set(gca,'Yticklabel',{'A', 'B', 'C','D','E','F','G','H','I','J'},'fontsize',10,'FontWeight','bold')
set(gca,'xticklabel',{'A', 'B', 'C','D','E','F','G','H','I','J'},'fontsize',10,'FontWeight','bold')
title 'Transition Matrix Heatmap';
nfig=nfig+1;


stateNames = ["A" "B" "C" "D" "E" "F" "G" "H" "I" "J"];
mc = dtmc(trans_matrix,'StateNames',stateNames);

figure(nfig)
graphplot(mc,'ColorEdges',true,'ColorNodes',true);


figure(nfig)
imagesc(trans_matrix_masked3);
colormap(flipud(bone));
colorbar;
axis square
h = gca;
h.XTick = 1:length(trans_matrix_masked3);
h.YTick = 1:length(trans_matrix_masked3);
set(gca,'Yticklabel',{'A', 'B', 'C','D','E','F','G','H','I','J'},'fontsize',10,'FontWeight','bold')
set(gca,'xticklabel',{'A', 'B', 'C','D','E','F','G','H','I','J'},'fontsize',10,'FontWeight','bold')
title 'Transition Matrix Heatmap';
nfig=nfig+1;


figure(nfig)
imagesc(trans_matrix_masked);
colormap(flipud(bone));
colorbar;
axis square
h = gca;
h.XTick = 1:length(trans_matrix_masked);
h.YTick = 1:length(trans_matrix_masked);
set(gca,'Yticklabel',{'A', 'B', 'C','D','E','F','G','H','I','J'},'fontsize',10,'FontWeight','bold')
set(gca,'xticklabel',{'A', 'B', 'C','D','E','F','G','H','I','J'},'fontsize',10,'FontWeight','bold')
title 'Transition Matrix Heatmap';
nfig=nfig+1;


%Read in the entire data
[ids_ok, ids_order] = findgroups(table.AnimalID);
[genotypes_ok, geno_order]= findgroups(table.Genotype);
[sex_ok, sex_order] = findgroups(table.Sex);
dat=[ids_ok genotypes_ok sex_ok table.AgeInMonth table.Weight];
n=5;%6

data2write=[];
for i=1:length(out)
    celldata=out{i};
    hiddenstates=celldata(1,:);
    v_ids =zeros(length(celldata(1,:)), 1);
    v_ids(1:length(celldata(1,:)), 1) = i;
    state2write=cell2mat(hiddenstates)';
    sex_out = cell2mat(celldata(3,:))';
    genotype_out = cell2mat(celldata(2,:))';
    age_out = cell2mat(celldata(4,:))';
    weights_out = cell2mat(celldata(5,:))';
   
    mat=[v_ids state2write sex_out age_out genotype_out weights_out];
    data2write = vertcat(data2write,mat);
end

dataTable_t = array2table(data2write, 'VariableNames', {'IDs', 'States', 'Sex1', 'Age1', 'Genotype1', 'Weight1'});
filename = 'results/trajectories10states/dataTable_TransitionStates.xlsx';
writetable(dataTable_t, filename, 'Sheet', 'MyNewSheet');

% Save marginal states
fid = fopen('results/trajectories10states/states2marg.txt', 'w');


engine = jtree_inf_engine(bnet2);
evidence = cell(5, 2);
%look for each state 1 to 10
states2marg = {};

for i=1:nstatesi
    evidence{1,2} = i;
    [engine, loglik] = enter_evidence(engine, evidence);
    %look genotype do it 2 3 4 5
    for jj=2:3
        marg = marginal_nodes(engine, jj);
        p = marg.T;
        states2marg{i,jj-1}=p;
        str2w='';
        for oo=1:length(p)
            if oo==1
                str2w=num2str(p(oo));
            else
                str2w=strcat(str2w,',',num2str(p(oo)));
            end
        end
        fprintf(fid, '%s\t', str2w);
    end
    for jj=4:5
        marg = marginal_nodes(engine, jj);
        p = marg.T;
        states2marg{i,jj-1}=[marg.mu marg.Sigma];
        vect2w=[marg.mu marg.Sigma];
        str2w='';
        for oo=1:length(vect2w)
            if oo==1
                str2w=num2str(vect2w(oo));
            else
                str2w=strcat(str2w,',',num2str(vect2w(oo)));
            end
        end
        fprintf(fid, '%s\t', str2w);
    end
    fprintf(fid, '\n');
end

fclose(fid);
save('results/trajectories10states/states2marg')


%B end
B_end={};
%C end
C_end={};
%F end
F_end={};
%I end
I_end={};

%J end
J_end={};
data2write=[];
for i=1:length(out)
    celldata=out{i};
    hiddenstates=celldata(1,:);
    v_ids =zeros(length(celldata(1,:)), 1);
    v_ids(1:length(celldata(1,:)), 1) = i;
    
    state2write=cell2mat(hiddenstates)';
    sex_out = cell2mat(celldata(3,:))';
    genotype_out = cell2mat(celldata(2,:))';
    age_out = cell2mat(celldata(4,:))';
    weights_out = cell2mat(celldata(5,:))';
    
    mat=[v_ids state2write sex_out age_out genotype_out weights_out];
    data2write = vertcat(data2write,mat);
  
    state2write_orig=state2write;
    if length(unique(state2write_orig))==1
        state_i_1=0;
    else
        last = state2write_orig(end);
        state2write_orig(end)=[];

        for o=1:length(state2write_orig)
            rev_v=flip(state2write_orig);
            if rev_v(o)~=last
                res=o;
                break
            end
        end
        state_i_1=state2write(length(state2write)-res);
    end  
    % B end
    if last==2
        endok=size(B_end,1);
        B_end{endok+1,1} = unique(genotype_out);
        B_end{endok+1,2} = unique(sex_out);
        B_end{endok+1,3} = age_out;
        B_end{endok+1,4} = weights_out;
 
    end
    %end state C
    if last==3
        endok=size(C_end,1);
        C_end{endok+1,1} = unique(genotype_out);
        C_end{endok+1,2} = unique(sex_out);
        C_end{endok+1,3} = age_out;
        C_end{endok+1,4} = weights_out;
     
    end
    %end state I
    if last==9
        endok=size(I_end,1);
        I_end{endok+1,1} = unique(genotype_out);
        I_end{endok+1,2} = unique(sex_out);
        I_end{endok+1,3} = age_out;
        I_end{endok+1,4} = weights_out;
    end
    %end State F
    if last==6
        endok=size(F_end,1);
        F_end{endok+1,1} = unique(genotype_out);
        F_end{endok+1,2} = unique(sex_out);
        F_end{endok+1,3} = age_out;
        F_end{endok+1,4} = weights_out;
    end
       %end State J
    if last==10
        endok=size(J_end,1);
        J_end{endok+1,1} = unique(genotype_out);
        J_end{endok+1,2} = unique(sex_out);
        J_end{endok+1,3} = age_out;
        J_end{endok+1,4} = weights_out
    end
  
end

% state B
figure(nfig)
subplot(2,3,1)
histogram(categorical(cell2mat(B_end(:,1))));
title('end B genotype')

subplot(2,3,4)
histogram(categorical(cell2mat(B_end(:,2))));
title('state B sex')
nfig=nfig+1;

% state C
figure(nfig)
subplot(2,3,1)
histogram(categorical(cell2mat(C_end(:,1))));
title('end C genotype')

subplot(2,3,4)
histogram(categorical(cell2mat(C_end(:,2))));
title('state C sex')

nfig=nfig+1;

% state F
figure(nfig)
subplot(2,4,1)
histogram(categorical(cell2mat(F_end(:,1))));
title('end F genotype')

subplot(2,4,5)
histogram(categorical(cell2mat(F_end(:,2))));
title('state F sex')
nfig=nfig+1;


% state I
figure(nfig)
subplot(2,3,1)
histogram(categorical(cell2mat(I_end(:,1))));
title('end J genotype')

subplot(2,3,4)
histogram(categorical(cell2mat(JI_end(:,2))));
title('state I sex')
nfig=nfig+1;

% state J
figure(nfig)
subplot(2,3,1)
histogram(categorical(cell2mat(J_end(:,1))));
title('end J genotype')

subplot(2,3,4)
histogram(categorical(cell2mat(J_end(:,2))));
title('state J sex')
nfig=nfig+1;

 
data2write2=[];
v_ii=0;
[data2write2,v_ii]=writeDataMatrix(v_ii,data2write2,2,B_end);
[data2write2,v_ii]=writeDataMatrix(v_ii,data2write2,3,C_end);
[data2write2,v_ii]=writeDataMatrix(v_ii,data2write2,6,F_end);
[data2write2,v_ii]=writeDataMatrix(v_ii,data2write2,9,I_end);
[data2write2,v_ii]=writeDataMatrix(v_ii,data2write2,10,J_end);


dataTable = array2table(data2write, 'VariableNames', {'IDs', 'States', 'Sex1', 'Age1', 'Genotype1', 'Weight1'});
filename = 'results/trajectories10states/dataTable_all.xlsx';
writetable(dataTable, filename, 'Sheet', 'MyNewSheet');
