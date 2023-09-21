function GUI_CreateTrialType

%% Saving colours for visualisations
palette = [88, 0, 0;
    121, 38, 11;
    156, 69, 17;
    189, 101, 26;
    221, 134, 41;
    245, 173, 82;
    254, 214, 147;
    255, 255, 224;
    187, 228, 207;
    118, 199, 190;
    62, 168, 166;
    32, 136, 136;
    7, 103, 105;
    0, 75, 75;
    0, 44, 45]./255;

contrast = palette([14, 6, 11, 4, 7, 10], :);
c_ix = repmat(1:6, 1, 10);
%% GENERATING EMPTY DATASETS TO STORE THE RESULTS
% Generating the save struct
Task.Titles = {};
Task.Lengths = {};
Task.FPS = [];

Loop.Titles = {};
Loop.Lengths = {};
Loop.Index = {};

% Generating the momentary storage place for epochs
epoch_title = [];
epoch_length = [];

% Generating momentary storage place for loop parameters
loop_name = [];
loop_repeats = [];
file_name = [];
save_loc = [];
save_path = [];

% Counter for determining in which stage of the trial we are
counter_loop = 1;

%% CREATING THE START WINDOW
% Will ask user to develop a protocol for the analyses
fig = uifigure;
fig.Name = "Trial";

% Indicating the first step
uicontrol("Parent", fig, "Style", "text", "String", ...
    "1 - Select frame rate", "FontWeight","bold", "Position", ...
    [10, 350, 120, 20]);

% Asking for frame rate
uicontrol("Parent", fig, 'Style', 'text', "String", ...
    "FPS:", "FontWeight","bold", "Position", ...
    [10, 327, 30, 20]);
uicontrol("Parent", fig, 'Style', 'edit', "Position", [50, 330, 90, 20], ...
    "Callback", @save_FPS)

% Indicating the second step
uicontrol("Parent", fig, "Style", "text", "String", ...
    "2 - Add epochs", "FontWeight","bold", "Position", ...
    [10, 280, 120, 20]);

% Creating the button to ask user for single epoch input
uicontrol('Parent', fig, "String", "Add new epoch", ...
    "Position", [10, 260, 120, 20], "Callback", @ask_epoch_input);

% Creating the button to ask the user for loop input
uicontrol('Parent', fig, "String", "Add new epoch loop", ...
    "Position", [10, 240, 120, 20], "Callback", @ask_loop_input);

% Creating the button to save
uicontrol("Parent", fig, "Style", "text", "String", ...
    "3 - Remember to save!", "FontWeight","bold", "Position", ...
    [10, 190, 120, 20]);
uicontrol("Parent", fig, "String", "Save!", "Style","pushbutton", ...
    "Position", [10, 170, 120, 20], "Callback", @saveTrial);


% Creating the empty axis
ax = uiaxes(fig, "Position", [150, 10, 350, 400]);

% Creating the cancel button
uicontrol("Parent", fig, "String", "Cancel", "Style", "pushbutton", ...
    "Position", [10, 10, 120, 20], "Callback", @cancel_design);

%% CREATING THE EPOCH WINDOW
input_epoch = dialog('Visible','off','Position',[500,500,250,150]);

% Add two input fields to the input window
uicontrol('Parent',input_epoch,'Style','edit','Position',[10,95,230,20], ...
    "Callback", @update_title);
uicontrol('Parent',input_epoch,'Style','edit','Position',[10,40,230,20], ...
    "Callback", @update_length);
uicontrol('Parent', input_epoch, 'Style', 'text', "String", ...
    "Enter epoch title", 'Position', [10, 120, 230, 20]);
uicontrol('Parent', input_epoch, 'Style', 'text', "String", ...
    "Enter epoch length in seconds", 'Position', [10, 65, 230, 20]);

% Add a "Save" button to the input window
uicontrol('Parent',input_epoch,'Style','pushbutton','String','Cancel', ...
    'Position',[160,10,70,25], "Callback", @cancel_epoch);

% Add a "Cancel" button to the input window
uicontrol('Parent',input_epoch,'Style','pushbutton','String','Save', ...
    'Position',[90,10,70,25], "Callback", @save_epoch);

%% CREATING THE EPOCH LOOP WINDOW
input_loop = dialog("Visible", "off", "Position", [500,500,270,100]);

% For interacting with the user
uicontrol("Parent", input_loop, "String", "Enter loop parameters", ...
    "Style", "text", "Position", [10, 67, 270, 30], "FontSize", 12, ...
    "FontWeight","bold");
uicontrol("Parent", input_loop, "String", "Enter loop name:", "Style", ...
    "text", "Position", [10, 55, 120, 20]);
uicontrol("Parent", input_loop, "String", "Number of repeats:", "Style", ...
    "text", "Position", [10, 35, 120, 20]);

% For user input
uicontrol("Parent", input_loop, "Style", 'edit', "Position", ...
    [140, 55, 130, 20], "Callback", @save_loop_name);
uicontrol("Parent", input_loop, "Style", 'edit', "Position", ...
    [140, 35, 130, 20], "Callback", @save_loop_repeats);

% Buttons
uicontrol("Parent", input_loop, "Style", "pushbutton", "String", ...
    "Add epoch", "Position", [10, 10, 80, 20], "Callback", @add_epoch);
uicontrol("Parent", input_loop, "Style", "pushbutton", "String", ...
    "Save", "Position", [100, 10, 80, 20], "Callback", @val_loop);
uicontrol("Parent", input_loop, "Style", "pushbutton", "String", ...
    "Cancel", "Position", [190, 10, 80, 20], "Callback", @cancel_loop);

    %% CREATING THE EPOCH IN LOOP INPUT WINDOW
    ep_loop_window = dialog('Visible','off','Position',[500,300,350,150]);

    % Add two input fields to the input window
    uicontrol('Parent', ep_loop_window,'Style','edit','Position',[10,95,230,20], ...
        "Callback", @update_title);
    uicontrol('Parent', ep_loop_window,'Style','edit','Position',[10,40,230,20], ...
        "Callback", @update_length);
    uicontrol('Parent', ep_loop_window, 'Style', 'text', "String", ...
        "Enter epoch title", 'Position', [10, 120, 230, 20]);
    uicontrol('Parent', ep_loop_window, 'Style', 'text', "String", ...
        "Enter epoch length in seconds", 'Position', [10, 65, 230, 20]);
    
    % Add a "Save" button to the input window
    uicontrol('Parent', ep_loop_window,'Style','pushbutton','String','Cancel', ...
        'Position',[160,10,70,25], "Callback", @cancel_epoch_loop);
    
    % Add a "Cancel" button to the input window
    uicontrol('Parent', ep_loop_window,'Style','pushbutton','String','Save', ...
        'Position',[90,10,70,25], "Callback", @save_epoch_loop);

    %% CREATING LOOP VALIDATION WINDOW
    val_loop_window = dialog("Visible", 'off', "Position", ...
        [800, 300, 150, 300]);

    val_end = uicontrol("Parent", val_loop_window, "Style", ...
        "pushbutton", "String", "Continue", "Position", [10, 10, 50, 20], ...
        "Callback", @close_val);

    uicontrol("Parent", val_loop_window, "Style", "pushbutton", ...
        "String", "Delete", "Position", [90, 10, 50, 20], ...
        "Callback", @erase_last);
    
    %% CREATING ERROR WINDOWS
    ERR_NoFPS = dialog("Visible", "off", "Position", ...
                [800, 300, 150, 150]);

    % Close button
    uicontrol("Parent", ERR_NoFPS, "Style","pushbutton", ...
        "String", "Continue", "Callback", @close_FPS_error, ...
        "Position", [50, 10, 50, 20]);
    
    % Text notifyong of the error
    uicontrol("Parent", ERR_NoFPS, "Style", 'text', "String", ...
        "WARNING - Input FPS before saving", "FontWeight", ...
        "bold", "Position", [10, 50, 130, 60]);

    %% CREATING THE SAVE WINDOW
    % Creating a window to ask for save inputs
    save_input_w = dialog("Visible", 'OFF', 'Position', ...
        [800, 300, 300, 150]);

    % Asking for save name
    uicontrol("Parent", save_input_w, "Style", 'text', 'String', ...
        "Enter task name (no spaces admitted)", ...
        "Position", [10, 110, 280, 20]);
    uicontrol("Parent", save_input_w, "Style", 'edit', "Position", ...
        [10, 95, 280, 20], "Callback", @save_taskname);
    
    % Asking for save path
    uicontrol("Parent", save_input_w, "Style", "text", "String", ...
        "Select path to save file", "Position", [10, 75, 280, 20], ...
        "Callback", @select_folder);
    uicontrol("Parent", save_input_w, "Style", "pushbutton", "String", ...
        "Select", "Position", [125 , 55, 50, 20], ...
        "Callback", @set_save_path);
    path_notif = uicontrol("Parent", save_input_w, "Style", "text", ...
        "Position", [10, 30, 280, 20], "Callback", @select_folder);
    % Create the final save button
    uicontrol("Parent", save_input_w, "Style", "pushbutton", "String", ...
        "Save", "Position", [125 , 10, 50, 20], "Callback", @save_task);
    % Create a cancel button
    uicontrol("Parent", save_input_w, "Style", "pushbutton", "String", ...
        "Cancel", "Position", [240 , 10, 50, 20], "Callback", ...
        @cancel_save);
    %% FUNCTIONS FOR THE MAIN WINDOW
    % Creating the function for opening the pop up window for epoch inputs
    function ask_epoch_input(~,~)
        input_epoch.Visible = 'On';
    end
    
    % ... and for the loop epoch inputs
    function ask_loop_input(~, ~)
        input_loop.Visible = 'On';
        Loop.Titles = {};
        Loop.Lengths = {};
        Loop.Index = {};
        counter_loop = 1;
    end

    % Visualising the task as it goes
    function view_task(~, ~)
        
        Lengths = str2double(string(Task.Lengths));
        % Setting the start and end of the trial
        Task.Start = [1, cumsum(Lengths(1, 1:end-1))+1];
        Task.End = cumsum(Lengths);
        base = zeros(length(Lengths), Task.End(1, length(Task.End)));
        for ix = 1:length(Lengths)
            base(ix, Task.Start(1, ix):Task.End(1, ix)) = 1;
            area(ax, base(ix, :), 'FaceColor', ...
                contrast(c_ix(1, ix), :))
            hold(ax, 'on')
        end
        legend(ax, string(Task.Titles))
        hold(ax, 'off')
        xlabel(ax, "Seconds")
        ylim(ax, [0, 1.3])
        set(ax,'ytick',[])

    end

    % Saving the FPSs
    function save_FPS(src, ~)
        Task.FPS = str2double(src.String);
    end
    
    % Saving the trial
    function saveTrial(~, ~)

        % Figuring out if the frame rate is empty
        if isempty(Task.FPS)
            ERR_NoFPS.Visible = 'on';

        else
            % Asking for the save parameters
            save_input_w.Visible = 'on';
        end
    end
    
    % Close error
    function close_FPS_error(~, ~)
        ERR_NoFPS.Visible = 'Off';
    end

    function cancel_design(~, ~)
        fig.Visible = 'off';
        closereq()
    end
    %% FUNCTIONS FOR EPOCH INPUTS
    % Function to save the title of the epoch
    function update_title(src, ~)
        epoch_title = src.String;

    end
    
    % Function to save the length of the epoch
    function update_length(src, ~)
        epoch_length = src.String;
    end

    
    % Function to save the parameters
    function save_epoch(~, ~)

        % Saving the title and length
        Task.Titles{length(Task.Titles)+1} = epoch_title;
        Task.Lengths{length(Task.Lengths)+1} = epoch_length;

        % Counter updating
        view_task

        % Closing the epoch input window
        input_epoch.Visible = 'Off';
    end

    % Function to cancel the epoch
    function cancel_epoch(~, ~)
        input_epoch.Visible = 'Off';
    end
%% CREATING FUNCTIONS FOR THE MAIN WINDOW OF THE LOOP GENERATION WINDOW
    % Momentary storage of the loop name
    function save_loop_name(src, ~)
        loop_name = src.String;
    end
    
    %... and of the number of repeats
    function save_loop_repeats(src, ~)
        loop_repeats = src.String;
    end
    
    function add_epoch(~, ~)
        ep_loop_window.Visible = 'On';
    end

    % Validating the loop output
    function val_loop(~, ~)
        val_loop_window.Visible = 'On';
        counter_loop = counter_loop - 1;
        update_validation
        
        % Changing the parameters of the close button
        val_end.String = "Cancel";
        val_end.Callback = @to_menu;

        % Creating a save button
        save_loop_btn = uicontrol("Parent", val_loop_window, "Style", ...
            "pushbutton", "String", "Save!", "Position", ...
            [90, 10, 50, 20], "Callback", @save_loop);
    end

    % Saving the loop
    function save_loop(~, ~)
        val_loop_window.Visible = 'Off';

        % Saving the lengths
        Task.Lengths = cat(2, Task.Lengths, ...
            repmat(Loop.Lengths, 1, str2double(loop_repeats)));
        
        % Saving the titles for each epoch
        Titles = string(repmat(Loop.Titles, 1, str2double(loop_repeats)));

        % Creating a loop constant part
        LoopK = repmat(convertCharsToStrings(loop_name), 1, ...
            str2double(loop_repeats)*length(Loop.Lengths));

        % Creating loop repeat index
        LoopIx = reshape((ones(str2double(loop_repeats), ...
            length(Loop.Lengths)).* [1:str2double(loop_repeats)].').', 1, []);
        Task.Titles = cat(2, Task.Titles, ...
            num2cell(Titles + " " + LoopK + string(LoopIx)));
        view_task

        input_loop.Visible = 'Off';

    end

    % Cancelling the loop
    function cancel_loop(~, ~)
        input_loop.Visible = 'Off';
    end

    %% AND FOR THE SUB-WINDOWS

    % Saving the epoch parameters
    function save_epoch_loop(~, ~)

        % Saving the results
        Loop.Titles{counter_loop} = epoch_title;
        Loop.Lengths{counter_loop} = epoch_length;
    
        % Displaying the register, and closing the epoch
        ep_loop_window.Visible = 'Off';
        update_validation
        input_loop.Visible = 'On';
        counter_loop = counter_loop + 1;

        
    end

    % Updating the validation window
    function update_validation(~, ~)
        % Generating the input for the register
        for ix = 1:counter_loop
            disp_str = strcat(Loop.Lengths{ix}, " seconds of ", ...
                Loop.Titles{ix});
            uicontrol("Parent", val_loop_window, "Style", 'text', ...
                "String", disp_str, "Position", ...
                [5, 300-ix*25, 150, 20])
        end
        val_loop_window.Visible = 'On';
        input_loop.Visible = 'On';

        ep_loop_window.Visible = 'Off';
    end

    function close_val(~, ~)
        val_loop_window.Visible = 'Off';
        ep_loop_window.Visible = 'On';
    end

    
    function erase_last(~, ~)
        if length(Loop.Titles) == 1
            Loop.Titles = {};
            Loop.Lengths = {};
        else
            Loop.Titles = Loop.Titles(:, 1:end-1);
            Loop.Lengths = Loop.Lengths(:, 1:end-1);
        end
        counter_loop = counter_loop-1;
        val_loop_window.Visible = 'Off';
    end

    % Cancelling
    function cancel_epoch_loop(~, ~)
        ep_loop_window.Visible = 'Off';
    end

    % Retunring to menu
    function to_menu(~, ~)

        % Clearing all windows
        input_epoch.Visible = 'Off';
        input_loop.Visible = 'Off';
        ep_loop_window.Visible = 'Off';

        % Clearing temporary files
        Loop.Titles = {};
        Loop.Lengths = {};
        Loop.Index = {};
        loop_name = [];
        loop_repeats = [];
        counter_loop = 1;
    end
    
    %% CREATING THE FUNCTIONS FOR THE SAVE WINDOW

    % Saving the filename for the output
    function save_taskname(src, ~)

        file_name = src.String;
    end
    
    % Selecting destination for the output
    function set_save_path(~, ~)
        save_loc = uigetdir([], ...
            "Select folder where you want to save the task");
        path_notif.String = strcat(save_loc(1:10), "(...)",...
            save_loc(end-20:end));
    end
    
    % Saving all
    function save_task(~, ~)

        % Generating the output path
        save_path = strcat(save_loc, "\", file_name, ".mat");
        
        % Generating the output
        Task.Frames = str2double(string(Task.Lengths))*Task.FPS;
        Task.Start = [1, cumsum(Task.Frames(1, 1:end-1))+1];
        Task.End(1, length(Task.End))
        base = zeros(length(Task.Frames), Task.End(1, length(Task.End)));
        for ix = 1:length(Task.Frames)
            base(ix, Task.Start(1, ix):Task.End(1, ix)) = 1;
        end

        Task.Pattern = base;
        
        % Saving
        save(save_path, "Task")
        fig.Visible = 'off';
        closereq()
        
              
    end

    % Cancelling save
    function cancel_save
        save_input_w.Visible = 'off';
    end

end

