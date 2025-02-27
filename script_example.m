% This script shows how to run the MATLAB_IterativeInputSelection_with_Rtree-c 
% toolbox for a test example. The algorithm ranks and selects the most relevant 
% variables with respect to an output of interest. 
% 
% For the script to work, the Regression tree package source code must have 
% been compiled, and the resulting mex files copied to the main directory. 
% This directory must then be added to the MATLAB PATH. Refer to INSTALL.txt 
% for further information on the steps to follow. 
%
%
%
% Copyright 2014 Stefano Galelli and Riccardo Taormina
%
% Prof. Galelli is Assistant Professor, Singapore University of Technology and Design
% stefano_galelli@sutd.edu.sg
% http://people.sutd.edu.sg/~stefano_galelli/index.html
%
% Riccardo Taormina is a Ph.D. candidate at the Hong Kong Polytechnic University
% riccardo.taormina@connect.polyu.hk
%
% This file is part of MATLAB_IterativeInputSelection_with_RTree-c.
% 
%     MATLAB_IterativeInputSelection_with_RTree-c is free software: you can redistribute 
%     it and/or modify it under the terms of the GNU General Public License 
%     as published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.     
% 
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with MATLAB_IterativeInputSelection_with_RTree-c.  
%     If not, see <http://www.gnu.org/licenses/>.
% 

%% Set workspace
clear
clc


%% Load and prepare data

% load data
load -ascii Friedman_dataset.txt

% rename data
data = Friedman_dataset;
clear Friedman_dataset

% definition of the calibration and validation data-set
subset_cal = data(1:180,:);      
subset_val = data(181:end,:);

% Set the parameters for the Extra-Trees
M    = 500; % number of extra trees in the forest
nmin = 5;   % number of points per leaf
k    = 10;  % Number of random cuts

% Create struct of parameters for rtree-c Extra-Trees
rtensparam                     = init_extra_trees();
rtensparam.nbterms             = M;
rtensparam.rtparam.nmin        = nmin;
rtensparam.rtparam.extratreesk = k;

% Build an ensemble of Extra-Trees and return the predictions on the
% training and test datasets
X1  = single(subset_cal(:,1:end-1));
Y1  = single(subset_cal(:,end));
ls  = int32(1:size(subset_cal,1));    
X2  = single(subset_val(:,1:end-1)); 

[finalResult_val var_imp ensemble finalResult_cal] =...
    rtenslearn_c(X1,Y1,ls,[],rtensparam,X2,0);

% Calculate the model performance in calibration and validation
R2 = Rt2_fit(subset_cal(:,end),finalResult_cal)  % 
R2 = Rt2_fit(subset_val(:,end),finalResult_val)  % 

% Graphical analysis
figure;
subplot(221)
plot(subset_cal(:,end),'.-'); hold on; plot(finalResult_cal,'.-r'); grid on;
axis([1 length(subset_cal) min(subset_cal(:,end)) max(subset_cal(:,end))]);
xlabel('time'); ylabel('output');
legend('measured','predicted');
title('calibration - trajectory');
subplot(222)
plot(finalResult_cal,subset_cal(:,end),'.'); grid on
axis([min(subset_cal(:,end)) max(subset_cal(:,end)) min(subset_cal(:,end)) max(subset_cal(:,end))]); 
xlabel('measured'); ylabel('predicted');
title('calibration - scatter plot');
subplot(223)
plot(subset_val(:,end),'.-'); hold on; plot(finalResult_val,'.-r'); grid on;
axis([1 length(subset_val) min(subset_val(:,end)) max(subset_val(:,end))]);
xlabel('time'); ylabel('output');
legend('measured','predicted');
title('validation - trajectory');
subplot(224)
plot(finalResult_val,subset_val(:,end),'.'); grid on
axis([min(subset_val(:,end)) max(subset_val(:,end)) min(subset_val(:,end)) max(subset_val(:,end))]);
xlabel('measured'); ylabel('predicted');
title('validation - scatter plot');


%% k-fold cross-validation

% Define the parameters for the cross-validation
ns   = 5; % number of folds
flag = 1; % if flag == 1, an ensemble is built on the whole dataset at the end of the cross-validation. 
          % Otherwise (flag == 0), such model is not built.

% Shuffle the data
data_sh = shuffle_data(data);

% Run the cross-validation
[model] = crossvalidation_extra_tree_ensemble(data_sh,M,k,nmin,ns,flag);

% Model performance in calibration and validation
model.cross_validation.performance.Rt2_cal_pred_mean  % 
model.cross_validation.performance.Rt2_val_pred_mean  % 

% Graphical analysis
figure;
subplot(211)
plot(data_sh(:,end),'.-'); hold on; plot(model.complete_model.trajectories,'.-r'); grid on;
axis([1 length(data_sh) min(data_sh(:,end)) max(data_sh(:,end))]);
xlabel('time'); ylabel('output');
legend('measured','predicted');
title('calibration - trajectory');
subplot(212)
plot(model.complete_model.trajectories,data_sh(:,end),'.'); grid on
axis([min(data_sh(:,end)) max(data_sh(:,end)) min(data_sh(:,end)) max(data_sh(:,end))]); 
xlabel('measured'); ylabel('predicted');
title('calibration - scatter plot');


%% Input ranking

% Shuffle the data
data_sh = shuffle_data(data);

% Run the ranking algorithm
[result_rank] = input_ranking(data_sh,M,k,nmin);

% Graphical analysis

% sort variables for bar plot
[temp,ixes] = sort(result_rank(:,2))
figure;
bar(result_rank(ixes,1));
xlabel('variable'); 
ylabel('normalized variable importance');
title('variable ranking - bar plot');


%% Iterative input selection
% Set the parameters for IIS
rpar.ns  = 5;  % number of folds for the cross-validation
rpar.p   = 5;  % number of SISO models evaluated at each iteration
rpar.epsilon = 0;  % tolerance
rpar.max_iter = 6;  % maximum number of iterations
verbose  = 1;  % 0 for silent run / 1 for verbose mode 
% set the parameters for the underlying Extra-trree model
rpar.M = 500;  % number of tree in the forest
rpar.nmin= 5;  % number of points per leaf
rpar.k =size(data,2)-1; %number of random cut = number of input signals

% Launch the IIS
result_iis = perform_IIS(data,rpar,flag,verbose)

% Report exit condition
disp(result_iis.exit_condition);

% Selected variables (by iteration)
% Determine the number of selected variables
if strcmp(result_iis.exit_condition,...
        'An input variable was selected twice') == 1
    nVariables = length(fieldnames(result_iis)) - 3;
else if strcmp(result_iis.exit_condition,...
        'The maximum number of iterations was reached') == 1
    nVariables = max_iter;
    else
        nVariables = length(fieldnames(result_iis)) - 2;
    end
end
% Selected variables, Cumulated R2 of the MISO model, R2 of the selected
% SISO model on the residual
[X,R2,R2_res] = summarize_IIS_result( result_iis);
deltaR2 = [R2(1) ; diff(R2)];
nVariables = length(X);

% Plotting
figure; 
bar(1:nVariables,deltaR2,'FaceColor','b'); hold on;
plot(R2,'o-','Color','k','LineWidth',...
    1, 'MarkerSize', 8, 'MarkerFaceColor', 'w',...
    'MarkerEdgeColor', 'k'); grid on;
axis([0.5 5.5 0 1.0]);
set(gca,'XTick',1:nVariables); 
set(gca,'XTickLabel', string(X),'Ylim',[0.00 1.0]);
xlabel('selected variables'); ylabel('R^2');
title('IIS');

%% Multiple runs of the IIS algorithm (with different shuffled datasets)- standard mode

% Define the parameters as before plus 
rpar.mult_runs = 10; % number of runs for the IIS algorithm               

% Run the IIS algorithm
for i = 1:rpar.mult_runs    
    fprintf( 'Run #%d\n',i );
    % Shuffle the data
    data_sh = shuffle_data(data);
    % run iis withouth fixed prior order
    results_iis_n(i) = iterative_input_selection(data_sh,rpar, 1, []);
    clear data_sh   
end

% Plot the results
[X,R2,R2_res] = summarize_IIS_result( results_iis_n);
figure, draw_R2( R2, gca);

draw_colorMap( X, R2);


%% Multiple runs of the IIS algorithm (with different shuffled datasets) - residual mode

% Define the parameters as before plus 
rpar.mult_runs = 10; % number of runs for the IIS algorithm               

% Run the IIS algorithm
for i = 1:rpar.mult_runs    
    fprintf( 'Run #%d\n',i );
    % Shuffle the data
    data_sh = shuffle_data(data);
    % run iis withouth by fixing the first q variables in the ranking (useful
    % to force to remove correlated variables)
    results_iis_n(i) = iterative_input_selection(data_sh,rpar, 1, [5]);
    clear data_sh   
end

% Plot the results
[X,R2,R2_res] = summarize_IIS_result( results_iis_n);
figure, draw_R2( R2, gca); 
xticklabels(string(X) );

draw_colorMap( X, R2);

% This code has been written by Stefano Galelli and Riccardo Taormina.
% Last updated by Dennis Zanutto