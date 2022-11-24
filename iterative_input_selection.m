
function [result] = iterative_input_selection(subset,iispar,Vflag,fxV,varargin)
    
    % This function implements the IIS algorithm
    %
    % subset   = observations
    % iispar   = struct containing the following parameters:
    %   M        = number of trees in the ensemble
    %   nmin     = minimum number of points per leaf
    %   ns       = number of folds in the k-fold cross-validation process
    %   p        = number of SISO models (it must be smaller than the number of
    %              candidate inputs).
    %   k        = number of random cuts, if empty set to the number of
    %              candidate variables.
    %   epsilon  = tolerance
    %   max_iter = maximum number of iterations
    % Vflag     = selection of the type of validation:
    %               1 = k-fold(default)
    %               2= repeated random sub-sampling
    % fxV      = column index of the subset containg the Variable from which
    %           you have to start the iterative procedure. During the first q
    %           iterations the variable is FiXed, first q elements of the miso
    %           model.
    %
    %
    % Output result   = structure containing the result for each iteration
    %
    %
    % Copyright 2014 Stefano Galelli and Matteo Giuliani
    % Assistant Professor, Singapore University of Technology and Design
    % stefano_galelli@sutd.edu.sg
    % http://people.sutd.edu.sg/~stefano_galelli/index.html
    % Research Fellow, Politecnico di Milano
    % matteo.giuliani@polimi.it
    % http://giuliani.faculty.polimi.it
    %
    % Please refer to README.txt for further information.
    %
    %
    % This file is part of MATLAB_IterativeInputSelection.
    %
    %     MATLAB_IterativeInputSelection is free software: you can redistribute
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
    %     along with MATLAB_IterativeInputSelection.
    %     If not, see <http://www.gnu.org/licenses/>.
    %
    
    %% 0) SET THE PARAMETERS
    M       = iispar.M;
    nmin    = iispar.nmin;
    if isempty( iispar.k )
        k   = size(subset,2)-1;
    else
        k   = iispar.k;
    end
    
    ns      = iispar.ns;
    p       = iispar.p;
    epsilon = iispar.epsilon;
    max_iter= iispar.max_iter;
    
    iP = inputParser;
    
    %This two variables are optional, so it works even if you don't insert
    %them, on the other hand the order is fixed, so 3rd element is alwasy Vflag
    %           'variable', 'default', 'check function'
    addOptional(iP,'Vflag',  1,     @(x) isnumeric(x) && isscalar(x) ); 
    addOptional(iP,'fxV',   [],     @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x)<= max_iter );
    
    % this is a couple Name Value for the function, use if during the procedure
    % you are intereseted in seeing the name of the variable rather than the
    % number representing it.
    addParameter(iP, 'Name', string(1:size(subset,2)-1), @(x) (isstring(x) || iscellstr(x)) && length(x) == size(subset,2)-1 );
    
    parse( iP, Vflag, fxV, varargin{:} );
    
    f           = iP.Results.Vflag;
    fxV         = iP.Results.fxV;
    listNames   = string(iP.Results.Name);
    
    % number of forced variable selected
    q       = length( fxV );
    
    % Initialize the counter and the exit condition flag
    iterC     = 1;    % iterations counter
    diff     = 1;    % exit flag
    
    % Re-define the subset matrix
    l = floor(size(subset,1)/ns);
    sel_row = 1:l*ns;
    
    % Define the MISO model output
    miso_output = subset(sel_row,end);
    miso_model  = [];
    
    % Define the set of candidate input variables
    input  = subset(sel_row,1:end-1);
    
    % Other variables to be initialized
    miso_input = [];  % initialize an empty set to store the input variables to be selected%
    
    
    %% 1) IIS ALGORITHM
    while (diff > epsilon) && (iterC <= max_iter)
        
        % Visualize the iteration
        fprintf('ITERATION:\n\t%d\n', iterC);
        
        % Define the output variable to be used during the ranking
        if iterC == 1
            rank_output = miso_output; % at the first iteration the MISO model output and ranking output are the same variable%
        else
            rank_output = residual;    % at the other iterations, the ranking output is the residual of the previous MISO model%
        end
        
        if iterC <= q
            % the first q element are not set according to the input ranking algorithm,
            % but from the user
            % ranking is fixed from the user, the selected feature is
            features = fxV(iterC);
            list = 1:size(input,2);
            [X, I] = sort( list'==features, 1, 'descend' );
            ranking = [X,I];
            result.iter(iterC).ranking = ranking;
            fprintf('Evaluating SISO model:\n\t%s\n', listNames(features));
            %evaluate the siso model, for consistency; k = 1 (SISO)
            if f == 1
                [siso_model] = crossvalidation_extra_tree_ensemble([input(:,features) rank_output],M,1,nmin,ns,0);
            else
                [siso_model] = repeatedRandomSubSamplingValidation_extra_tree_ensemble([input(:,features) rank_output],M,1,nmin,ns,0);
            end
            performance = siso_model.cross_validation.performance.Rt2_val_pred_mean;
            
            result.iter(iterC).SISO = [features, performance];
           
            % Choose the SISO model with the best performance
            val = performance;
            best_siso_input = features;
            
        else
            % Define the ranking matrix
            matrix_ranking = [input rank_output];
            
            % Run the feature ranking; k = number of cand var 
            [ranking] = input_ranking(matrix_ranking,M,k,nmin);
            result.iter(iterC).ranking = ranking;
            
            % Select and cross-validate p SISO models (the first p-ranked models); k = 1 (SISO)
            features = ranking(1:p,2);                             % p features to be considered
            performance = zeros(p,1);	                           % initialize a vector for the performance of the p SISO models%
            fprintf('Evaluating SISO models:\n');
            fprintf('\t%s\n', listNames(features) );
            for i = 1:p
                % to print at runtime the number of the model(this for cicle is the one
                % that takes the most time without any info)
                lineLenght = fprintf('%d/%d\n', i, p);
                if f == 1
                    [siso_model] = crossvalidation_extra_tree_ensemble([input(:,features(i)) rank_output],M,1,nmin,ns,0);
                else
                    [siso_model] = repeatedRandomSubSamplingValidation_extra_tree_ensemble([input(:,features(i)) rank_output],M,1,nmin,ns,0);
                end
                performance(i) = siso_model.cross_validation.performance.Rt2_val_pred_mean;
                
                % to avoid using multiple lines
                fprintf(repmat('\b', 1,lineLenght));
            end
            result.iter(iterC).SISO = [features, performance];
            
            % Choose the SISO model with the best performance
            [val,idx_siso] = max(performance);
            best_siso_input = features(idx_siso);
            
        end
        
        result.iter(iterC).best_SISO = [best_siso_input val];
        fprintf('Selected variable:\n\t%d %4.2f %s\n\n', best_siso_input, val, listNames(best_siso_input) );
        
        % Check the exit condition
        if any( miso_input == best_siso_input )
            result.exit_condition = 'An input variable was selected twice';
            result.iters_done = iterC;
            result.iters_valid = iterC-1;
            return
        end
        
        % Build a MISO model with the selected inputs; k = iterC (MISO, the number of variables in the matrix is equal to the number of iterations done)
        % Save old for performance comparison
        miso_model_old = miso_model;
        fprintf('Evaluating MISO model:\n');
        miso_input = [miso_input best_siso_input]; %#ok<AGROW>
        if f==1
            [miso_model] = crossvalidation_extra_tree_ensemble([input(:,miso_input) miso_output],M,iterC,nmin,ns,1);
        else
            [miso_model] = repeatedRandomSubSamplingValidation_extra_tree_ensemble([input(:,miso_input) miso_output],M,iterC,nmin,ns,1);
        end
        result.iter(iterC).MISO = miso_model;
        fprintf('\t%4.2f\n\n', miso_model.cross_validation.performance.Rt2_val_pred_mean);
        
        % Evaluate the performance of the MISO model and calculate the
        % difference with respect to the previous MISO model
        if iterC == 1   % at the first iteration, use a default value
            diff = 1;
        else
            diff = miso_model.cross_validation.performance.Rt2_val_pred_mean - miso_model_old.cross_validation.performance.Rt2_val_pred_mean;
        end
        
        % Compute the MISO model residual
        residual = miso_output - miso_model.complete_model.trajectories;
        
        
        % Check the exit condition
        if iterC == max_iter
            result.exit_condition = 'The maximum number of iterations was reached';
            result.iters_done = iterC; 
            result.iters_valid = iterC;
        end
        if diff <= epsilon
            result.exit_condition = 'The tolerance epsilon was reached';
            result.iters_done = iterC;
            result.iters_valid = iterC-1;
        end
        
        % Update the counter at the end! 
        iterC = iterC + 1;
        
    end
    
end
% This code has been written by Stefano Galelli, Matteo Giuliani
% Updated by Dennis Zanutto