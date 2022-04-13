
function [result] = iterative_input_selection(subset,iispar,Vflag,fxV,varargin)
    
    % This function implements the IIS algorithm
    %
    % subset   = observations
    % iispar   = struct containing the following parameters:
    % M        = number of trees in the ensemble
    % nmin     = minimum number of points per leaf
    % ns       = number of folds in the k-fold cross-validation process
    % p        = number of SISO models (it must be smaller than the number of
    %            candidate inputs).
    % epsilon  = tolerance
    % max_iter = maximum number of iterations
    % Vflag     = selection of the type of validation,
    %               1 = k-fold(default)
    %               2= repeated random sub-sampling
    % fxV      = column index of the subset containg the Variable from which
    %           you have to start. First q iterations the variable is FiXed,
    %           first q elements of the miso model.
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
    
    % 0) SET THE PARAMETERS
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
    
    addOptional(iP,'Vflag', 1, @(x) isnumeric(x) && isscalar(x) ); 
    addOptional(iP,'fxV',[], @(x) isnumeric(x) && isvector(x) && length(x)<= max_iter );
    
    addParameter(iP, 'Name', string(1:size(subset,2)-1), @(x) (isstring(x) || iscellstr(x)) && length(x) == size(subset,2)-1 );
    
    parse( iP, Vflag, fxV, varargin{:} );
    
    f           = iP.Results.Vflag;
    fxV         = iP.Results.fxV;
    listNames   = string(iP.Results.Name);
    
    % number of forced variable selected
    q       = length( fxV );
    
    % Initialize the counter and the exit condition flag
    iter     = 1;    % iterations counter
    diff     = 1;    % exit flag
    
    % Re-define the subset matrix
    l = floor(size(subset,1)/ns);
    subset = subset(1:l*ns,:);
    
    % Define the MISO model output
    miso_output = subset(:,end);
    miso_model  = [];
    
    % Define the set of candidate input variables
    input  = subset(:,1:end-1);
    
    % Other variables to be initialized
    miso_input = [];  % initialize an empty set to store the input variables to be selected%
    
    
    % 1) IIS ALGORITHM
    while (diff > epsilon) && (iter <= max_iter)
        
        % Visualize the iteration
        disp('ITERATION:'); disp(iter);
        iterS = ['iter_', num2str(iter)];
        
        % Define the output variable to be used during the ranking
        if iter == 1
            rank_output = miso_output; % at the first iteration the MISO model output and ranking output are the same variable%
        else
            rank_output = residual;    % at the other iterations, the ranking output is the residual of the previous MISO model%
        end
        
        if iter <= q
            % the first q element are not set according to the input ranking algorithm,
            % but from the user
            % ranking is fixed from the user, the selected feature is
            features = fxV(iter);
            list = 1:size(input,2);
            [X, I] = sort( list'==features, 1, 'descend' );
            ranking = [X,I];
            result.(iterS).ranking = ranking;
            %disp(ranking);
            
            %evaluate the siso model, for consistency
            if f == 1
                [siso_model] = crossvalidation_extra_tree_ensemble([subset(:,features) rank_output],M,1,nmin,ns,0);
            else
                [siso_model] = repeatedRandomSubSamplingValidation_extra_tree_ensemble([subset(:,features) rank_output],M,1,nmin,ns,0);
            end
            performance = siso_model.cross_validation.performance.Rt2_val_pred_mean;
            
            result.(iterS).SISO = [features, performance];
            disp('Tested SISO:');
            disp([features,  performance, listNames(features)]);
            
            % Choose the SISO model with the best performance
            val = performance;
            best_siso_input = features;
            result.(iterS).best_SISO = [best_siso_input val];
            disp('Select variable:'); disp([best_siso_input, listNames(best_siso_input)]);
            
        else
            % Define the ranking matrix
            matrix_ranking = [input rank_output];
            
            % Run the feature ranking
            %disp('Ranking:');
            [ranking] = input_ranking(matrix_ranking,M,k,nmin);
            result.(iterS).ranking = ranking;
            %disp(ranking);
            
            % Select and cross-validate p SISO models (the first p-ranked models)
            disp('Evaluating SISO models:');
            features = ranking(1:p,2);                             % p features to be considered
            performance = zeros(p,1);	                           % initialize a vector for the performance of the p SISO models%
            for i = 1:p
                lineLenght = fprintf('%d/%d\n', i, p);
                if f == 1
                    [siso_model] = crossvalidation_extra_tree_ensemble([subset(:,features(i)) rank_output],M,1,nmin,ns,0);
                else
                    [siso_model] = repeatedRandomSubSamplingValidation_extra_tree_ensemble([subset(:,features(i)) rank_output],M,1,nmin,ns,0);
                end
                performance(i) = siso_model.cross_validation.performance.Rt2_val_pred_mean; 
                fprintf(repmat('\b', 1,lineLenght));
            end
            result.(iterS).SISO = [features, performance];
            disp('Tested SISO:');
            disp([features,  performance, listNames(features)]);
            
            % Choose the SISO model with the best performance
            [val,idx_siso] = max(performance);
            best_siso_input = features(idx_siso);
            result.(iterS).best_SISO = [best_siso_input val];
            disp('Select variable:'); disp([best_siso_input, listNames(best_siso_input)]);
        end
        
        % Check the exit condition
        if any( miso_input == best_siso_input )
            result.exit_condition = 'An input variable was selected twice';
            result.iters_done = iter;
            return
        end
        
        % Build a MISO model with the selected inputs
        % Save old for performance comparison
        miso_model_old = miso_model;
        disp('Evaluating MISO model:');
        miso_input = [miso_input best_siso_input]; %#ok<AGROW>
        k = length(miso_input);
        if f==1
            [miso_model] = crossvalidation_extra_tree_ensemble([subset(:,miso_input) miso_output],M,k,nmin,ns,1);
        else
            [miso_model] = repeatedRandomSubSamplingValidation_extra_tree_ensemble([subset(:,miso_input) miso_output],M,k,nmin,ns,1);
        end
        result.(iterS).MISO = miso_model;
        disp(miso_model.cross_validation.performance.Rt2_val_pred_mean);
        
        % Evaluate the performance of the MISO model and calculate the
        % difference with respect to the previous MISO model
        if iter == 1   % at the first iteration, use a default value
            diff = 1;
        else
            diff = miso_model.cross_validation.performance.Rt2_val_pred_mean - miso_model_old.cross_validation.performance.Rt2_val_pred_mean;
        end
        
        % Compute the MISO model residual
        residual = miso_output - miso_model.complete_model.trajectories;
        
        % Update the counter
        iter = iter + 1;
        
        % Check the exit condition
        if iter > max_iter
            result.exit_condition = 'The maximum number of iterations was reached';
            result.iters_done = iter;
        end
        if diff <= epsilon
            result.exit_condition = 'The tolerance epsilon was reached';
            result.iters_done = iter;
        end
        
    end
    
end
% This code has been written by Stefano Galelli, Matteo Giuliani
% Updated by Dennis Zanutto