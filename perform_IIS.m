function result = perform_IIS(data,iispar,Vflag,verbose)

    % This function is a wrapper around iterative_input_selection.m,
    % the function implementing the IIS technique using Extra-Trees.
    %
    %
    % data              = dataset;
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
    % verbose           = 0 for silent run. 1 for verbose mode
    %
    % Outputs
    % result   = structure containing the result for each iteration
    

% 0) check if p <= number of attributes
natt = size(data,2)-1;

if iispar.p > natt
    error(['The number of SISO models evaluated',...
        'has to be < number of candidate inputs'])
end

% 1)  Launch IIS algorithm

% Shuffle the data
data_sh = shuffle_data(data);

% Run the IIS algorithm
if verbose == 0
    evalc('result = iterative_input_selection(data_sh,iispar,Vflag,[])');
else
    result = iterative_input_selection(data_sh,iispar,Vflag,[]);
end


% This code has been written by Riccardo Taormina.
% Updated by Dennis Zanutto
