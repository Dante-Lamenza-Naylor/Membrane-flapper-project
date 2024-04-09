classdef solver_basket < handle
    properties
        basket
        dt
        dt_trigger
        key
        size
        relevant_indices
        rewinds
        times
        time_size
        basket_size
        last_time

        STM
        STM_size
        STM_dt_trigger

        LTM
        LTM_size
        LTM_dt_trigger
        LTM_dt_tracker
        
        tmax

        parameters

        figures
    end

    methods

        function obj = solver_basket(trigger_input)
            obj.basket = {};
            obj.dt_trigger = trigger_input;
            obj.dt = 0;
            obj.size = 0;
            obj.rewinds = 0;
            obj.basket_size = 0;
            obj.key = 0;
        end

        function SetUpMemory(obj,STM_dt, LTM_dt, tmax, rows)
            obj.tmax = tmax;

            obj.STM_dt_trigger = STM_dt;
            obj.LTM_dt_trigger = LTM_dt;

            LTM_whole_size = ceil(tmax / LTM_dt);
            % STM_whole_size = ceil(LTM_dt / STM_dt) + 4;
            STM_whole_size = 100;

            obj.LTM = cell(rows,LTM_whole_size);
            obj.STM = cell(rows,STM_whole_size);

            obj.LTM_size = 1;
            obj.LTM_dt_tracker = 0;
            
%             for i = (1:rows)
%                 obj.LTM{i,1} = 0;
%             end
        end

        function update_time(obj,t_input)
            % Pushes time into basket. Rewinds timeline if necessary.
            if isempty(obj.basket)
                obj.basket{1,1} = t_input;
                obj.dt = obj.dt + t_input;
                obj.size = 1;
                obj.basket{2,1} = 0; %Left disturbance field
                obj.basket{3,1} = 0; %Right disturbance field
                obj.basket{4,1} = 0; %Left Potential
                obj.basket{5,1} = 0; %Right potential
                obj.relevant_indices(1) = 1;
           
            % If t_input is less than any previous t in the list....
            else
                N = obj.size;
                change_dt = true;
                while t_input <= obj.basket{1,N}
                    if change_dt
                        dterror = obj.basket{1,N} - t_input;
                        obj.dt = mod(obj.dt - dterror,obj.dt_trigger);%(obj.basket{1,N}-obj.basket{1,N-1});                        
                        change_dt = false;
                    end
                    obj.basket(:,N) = [];
                    N = N-1;
                    obj.size = obj.size - 1;
                    if obj.size < obj.relevant_indices(length(obj.relevant_indices))
                        obj.relevant_indices(length(obj.relevant_indices)) = [];
                         if obj.dt < obj.dt_trigger * 1e-2
                            obj.dt = obj.dt_trigger; % Prevent the solver from avoiding the trigger
                         end
                        obj.rewinds = obj.rewinds + 1;
                    end
                end
                obj.basket{1,N+1} = t_input;
                obj.dt = obj.dt + (obj.basket{1,N+1} - obj.basket{1,N});
                obj.size = obj.size+1;
                obj.check_dt
            end

        end

        function new_update_time(obj,t_input)
            % Pushes time into basket. Rewinds timeline if necessary.
            if obj.basket_size == 0
                obj.basket{1,1} = t_input;
%                 obj.times(1) = 0;
                obj.dt = obj.dt + t_input;
                obj.time_size = 1;
                obj.basket_size = 1;
%                 obj.basket{2,1} = 0; %Left disturbance field
%                 obj.basket{3,1} = 0; %Right disturbance field
%                 obj.basket{4,1} = 0; %Left Potential
%                 obj.basket{5,1} = 0; %Right potential
                obj.relevant_indices(1) = 1;
                obj.last_time = 0;
                
            else
                Ntime = obj.time_size;
                Nbasket = obj.basket_size;
                change_dt = true;
                % While input time is less than the latest remembered
                % time...
                while t_input < obj.last_time
                    % If dt is supposed to change...
                    if change_dt
                        % calculate a new dt, remembering to do it modulo
                        % the dt trigger, and signal that dt is no longer
                        % supposed to change
                        dterror = obj.last_time - t_input;
                        obj.dt = mod(obj.dt - dterror,obj.dt_trigger);%(obj.basket{1,N}-obj.basket{1,N-1});                        
                        change_dt = false;
                    end
                    % If the input time is less than the latest time where
                    % the basket was supposed to do something special...
                    while t_input <= obj.basket{1,Nbasket}
                        % Delete what the basket remembers after that point
                        % in time and set sizes accordingly
%                         obj.basket(:,Nbasket) = [];
                        obj.basket_size = obj.basket_size - 1;
                        Nbasket = Nbasket - 1;
                        % Delete the memory of 'special' indices of the times
                        % matrix with all times past that point
                        obj.relevant_indices(length(obj.relevant_indices)) = [];
                         if obj.dt < obj.dt_trigger * 1e-2
                            obj.dt = obj.dt_trigger; % Prevent the solver from avoiding the trigger
                         end
                         % Keep track of how many rewinds happened
                        obj.rewinds = obj.rewinds + 1;
                    end
                    % Set time sizes accordingly - indeed, in every case,
                    % if this tree is triggered, times will lose one element
                    Ntime = Ntime-1;
                    obj.time_size = obj.time_size - 1;
                    obj.last_time = t_input;
                end
                % Add new t_input to end of times matrix, incriment dt,
                % adjust sizes accordingly, and check to see if the basket
                % should signal the 'key' flag
%                 obj.times(1,Ntime+1) = t_input;
%                 (obj.times(Ntime+1) - obj.times(Ntime));
%                 obj.time_size = obj.time_size+1;
                obj.dt = obj.dt + (t_input - obj.last_time);
                obj.last_time = t_input;
                obj.new_check_dt(t_input);
            end

        end


        function check_dt(obj)
            % Returns True/False, checks if dt_trigger condition is met.
            if obj.dt >= obj.dt_trigger
                obj.key = true;
                obj.dt = 0;
                obj.relevant_indices(length(obj.relevant_indices) + 1) = obj.size;
            else
                obj.key = false;
            end
        end

        function new_check_dt(obj,t_input)
            % Returns True/False, checks if dt_trigger condition is met.
            if obj.dt >= obj.dt_trigger
                obj.key = true;
                obj.dt = mod(obj.dt,obj.dt_trigger);
                obj.relevant_indices(length(obj.relevant_indices) + 1) = obj.time_size;
                obj.basket_size = obj.basket_size + 1;
                obj.basket{1,obj.basket_size} = t_input;
            else
                obj.key = false;
            end
        end

        function M_update_time(obj,t_input)
            % Pushes time into basket. Rewinds timeline if necessary.
            if isempty(obj.STM_size)
                obj.STM{1,1} = t_input;
%                 obj.times(1) = 0;
                obj.dt = obj.dt + t_input;
                obj.STM_size = 1;
%                 obj.basket_size = 1;
%                 obj.basket{2,1} = 0; %Left disturbance field
%                 obj.basket{3,1} = 0; %Right disturbance field
%                 obj.basket{4,1} = 0; %Left Potential
%                 obj.basket{5,1} = 0; %Right potential
%                 obj.relevant_indices(1) = 1;
                obj.last_time = 0;
                
            else
%                 Ntime = obj.time_size;
%                 Nbasket = obj.basket_size;
                change_dt = true;
                % While input time is less than the latest remembered
                % time...
                while t_input < obj.last_time
                    % If dt is supposed to change...
                    if change_dt
                        % calculate a new dt, remembering to do it modulo
                        % the dt trigger, and signal that dt is no longer
                        % supposed to change
                        dterror = obj.last_time - t_input;
                        obj.dt = mod(obj.dt - dterror,obj.dt_trigger);
                        obj.LTM_dt_tracker = mod(obj.LTM_dt_tracker - dterror, obj.LTM_dt_trigger);         
                        change_dt = false;
                    end
                    % If the input time is less than the latest time where
                    % the basket was supposed to do something special...
                    while t_input <= obj.STM{1,obj.STM_size}
                        while obj.LTM_size > 1 && t_input <= obj.LTM{1,obj.LTM_size-1}
                            obj.LTM_size = obj.LTM_size - 1;
                        end
                        % Delete what the basket remembers after that point
                        obj.STM_size = obj.STM_size - 1;
                        % Prevent the solver from avoiding the trigger
                         if obj.dt < obj.dt_trigger * 1e-2
                            obj.dt = obj.dt_trigger; 
                         end
                         % Keep track of how many rewinds happened
                        obj.rewinds = obj.rewinds + 1;
                    end
%                     % Set time sizes accordingly - indeed, in every case,
%                     % if this tree is triggered, times will lose one element
%                     Ntime = Ntime-1;
%                     obj.time_size = obj.time_size - 1;
                    obj.last_time = t_input;
                end
                % Add new t_input to end of times matrix, incriment dt,
                % adjust sizes accordingly, and check to see if the basket
                % should signal the 'key' flag
%                 obj.times(1,Ntime+1) = t_input;
%                 (obj.times(Ntime+1) - obj.times(Ntime));
%                 obj.time_size = obj.time_size+1;
                obj.dt = obj.dt + (t_input - obj.last_time);
                obj.LTM_dt_tracker = obj.LTM_dt_tracker + (t_input - obj.last_time);
                obj.last_time = t_input;
                obj.STMCommit(t_input);
            end

        end

        function STMsave(obj,row,thing)
            if obj.key
                obj.STM{row,obj.STM_size} = thing;
            end
        end

        function STMCommit(obj,t_input)
            % Returns True/False, checks if dt_trigger condition is met.
            if obj.dt >= obj.dt_trigger
                obj.key = true;

                obj.dt = mod(obj.dt,obj.dt_trigger);
                
                if obj.STM_size == length(obj.STM(1,:))
                    obj.STM(:,1:end-1) = obj.STM(:,2:end);
                else
                    obj.STM_size = obj.STM_size + 1;
                end
                
                obj.STM{1,obj.STM_size} = t_input;
                
                if (obj.LTM_dt_tracker >= obj.LTM_dt_trigger)
                    obj.LTMCommit;
                    obj.LTM_dt_tracker = mod(t_input,obj.LTM_dt_trigger);
                end

            else
                obj.key = false;
            end
        end

        function LTMCommit(obj)
            obj.LTM(:,obj.LTM_size) = obj.STM(:,obj.STM_size-1);
            obj.LTM{1,obj.LTM_size} = obj.STM{1,obj.STM_size};              % Doing a bit of cheatsy doodling here. Chiefly to fix issues with time repeating once a rewind comes.
            obj.LTM_size = obj.LTM_size + 1;                                
%             obj.DumpArray;
        end

        function DumpArray(obj)
            obj.STM(:,1:3) = obj.STM(:,obj.STM_size-2:obj.STM_size);
            obj.STM_size = 3;
        end


    end
end