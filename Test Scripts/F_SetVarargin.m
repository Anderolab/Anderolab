function [SpikeThreshold] = F_SetVarargin(FilteredTrace, varargin)
% C_SpikeThreshold
    % AIM - Compute the threshold for peak detection of a denoised signal.
    % OPTIONAL INPUTS
        % 'Method' - How the threshold is computed
            % 'CI' (DEFAULT) - Using the 95% confidence interval
            % 'Median' - Number of SDs away from the median
                % OPTIONAL INPUT = 'Threshold'
                    % Number of SDs away from the median
                    % Default = 3;
            % 'Mean' - Number of SDs away from the mean
                % OPTIONAL INPUT = 'Threshold'
                    % Number of SDs away from the median
                    % Default = 3;
            % 'Quian'
        % 'Visualise' - Whether the figure is plotted or not
            % True or False
    %% Within function vars
    % Set the name for the optional arguments
    OpArgs = ["Method", "Visualise", "Threshold"];
    
    % Then set the expected imput
        % When a specific string is needed write the list of option
    ValsOptions = {["CI","Median", "Mean", "Quian"]; "logical"; "double"};

    % Set the default values if the user does not provide any
    % specifications
    ValsDefault = {'CI', false, 3};

    %% Storage of variables
    ArgLoco = repelem(1, length(OpArgs)); % Order of specified variables
    % Storage of the specified options as a dictionary for ease of use
    ValsOptions = dictionary(OpArgs.', ValsOptions);
            
    %% Identifying errors
    % If there are optional arguments
    if ~isempty(varargin)

    % Checking correct number of arguments
        if rem(length(varargin), 2)~=0
            error('Incorrect number of input arguments.')
        end

    % Reading the inputs
        KeyWords = string(varargin(1:2:end))
        OpArgs

    % Verifying that the argument names are correct
        if sum(~contains(KeyWords, OpArgs))> 0
            error(strcat('Unrecognized input ', ...
                string(KeyWords(find(contains(KeyWords, OpArgs)== 0, 1)))))
        end
        
    %% Determining the location of each argument
        % Counter to determine location of argument
        c = 2;

        % Iterating per argument name
        for Word = KeyWords

            % And storing the answer location
            ArgLoco(OpArgs == Word) = c;
            c = c+2;
        end
        
        % Extracting user answers according to input order
        ArgValues = varargin(ArgLoco);

        % Replacing non-defined arguments with the default value
        ArgValues(ArgLoco == 1) = ValsDefault(ArgLoco == 1);

        % Converting to dictionary format for ease of use
        ArgValues = dictionary(OpArgs, ArgValues)

        % Checking if the argument given meets the specifications
        for Arg = OpArgs
            Ops = ValsOptions(Arg);     Ops = Ops{:};
            if class(Ops) == 'string'
                ArgValues(Arg)
                if ~contains(ArgValues(Arg), Ops)
                    error(strcat("Unrecognised value ", ArgValues(Arg), ...
                        " for ", Arg, '.'))
                end
            else class(Ops) == 'double';
                if length(Ops) > 1
                    Val = ArgValues(Arg);       Val = Val{:}
                    if string(class(Val)) ~= "logical"
                        error(strcat("Expected logical for ", Arg, '.'))
                    end
                end
            end
        end
        
    else
        ArgValues = ValsDefault;
    end
    
end
