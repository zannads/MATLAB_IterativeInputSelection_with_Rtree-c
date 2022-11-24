function [X, R2, R2_res] = summarize_IIS_result( iis_res, Nvar )
    
    % [X, R2, R2_res] = visualize_inputSel( iis_res, Nvar )
    %
    % This function summarizes the results of multiple runs of the IIS algorithm to
    % analyze the input ranking, the variability of the final MISO models'
    % performance and the variability during the process over different runs.
    %
    % input:
    %   iis_res     = cell array containing the output of the
    %               iterative_input_selection funtion
    %   Nvar        = number of candidate input variables (in theory this is
    %               equal to the number of columns of the dataset, in practice
    %               it can be reduced to max_iter or not inserted and
    %               calculated automatically from iis_res)
    % output:
    %   X           = matrix containing the indices of the selected input for
    %               each run (on the columns)
    %   R2          = corresponding model performance
    %   R2_res      = corresponding SISO model performance on the residual of
    %               the previous iteration.
    %
    %
    % Copyright 2014 Stefano Galelli and Matteo Giuliani
    % Assistant Professor, Singapore University of Technology and Design
    % stefano_galelli@sutd.edu.sg
    % http://people.sutd.edu.sg/~stefano_galelli/index.html
    % Research Fellow, Politecnico di Milano
    % matteo.giuliani@polimi.it
    % http://home.deib.polimi.it/giuliani/
    %
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
    
    if nargin == 1
        Nvar = max( cat(1, iis_res.iters_valid) );
    end
    Nrun = length( iis_res );
    
    % extract results
    X           = nan(Nvar,Nrun);
    R2          = nan(Nvar,Nrun);
    R2_res      = nan(Nvar,Nrun);
    
    for i = 1:Nrun
      
        for j = 1: iis_res(i).iters_valid 
            X(j,i)      = iis_res(i).iter(j).best_SISO(1);
            R2(j,i)     = iis_res(i).iter(j).MISO.cross_validation.performance.Rt2_val_pred_mean;
            R2_res(j,i) = iis_res(i).iter(j).best_SISO(2);
        end
    end
    
end
% This code has been written by Matteo Giuliani.
% Updated by Dennis Zanutto