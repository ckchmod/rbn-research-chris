classdef boolCellGrid_v2 < matlab.mixin.Copyable
    %   Boolean cell grid - this function implements a grid of cells modeled
    % as Random Boolean Networks, with randomized but uniform intra- and
    % inter-cellular connections
    %   This class also supports several different topologies, as defined
    % by the 'neighbors' functions
    
    properties
        % Supplied by caller
        topology    % Options: line, symmetric, single (hex in progress)
        numCells    % Scalar, number of total cells
        numGenes    % How many nodes there are within a cell
        k           % How many intracellular connections each cell has
        p           % Probability of function value taking on value of 1
        bandwidth   % How many cells communicate intercellularly
        initState   % Initial State
        initTruth   % Initial Truth Table
        initVar     % Initial Wiring Table
        perturb     % Perturbation
        %pagerank    % Pagerank
        
        % Generated in the constructor
        allCells    % Array of all the cells
        neighbors   % Array whose rows are the neighbors of that cell
        timenow     % Current time
        initTtable  % Truth table
        initvarF    % Intracellular connectivity
        initStates  % Initial States
        initInCells % Which cells receive input
        initOutCells% Which cells produce output
        criticality % Criticality measured by 2Kp(1-p)=1
        crit_val    % LHS of critcality value
        %page_rank   % Speed up constant (default .95)
        E_disordered % energy
        E_magnetization %average magnetization
        J           % Hamiltonian
        temperature
        h_interac
        
        movie_end_steps
        
        % For comparison with the 'drosophila.m' code
        allStates   % All states        
    end
    
    methods        
        % Constructor
        function obj = boolCellGrid_v2(topology,numCells,numGenes,k,p,bandwidth, initState, initTruth, initVar, perturb, J, temperature, h_interac)
            
            rng('shuffle');            
            obj.topology  = topology;
            obj.numCells  = numCells;
            obj.numGenes  = numGenes;
            obj.k         = k;
            obj.p         = p;
            obj.bandwidth = bandwidth;
            obj.perturb = perturb;
            
            % Get the randomized initial states
            if isempty(initState)
                initialStates = genInit(obj, topology);
            else
                dim = size(initState);
                assert(isempty(find(...
                    ~( (initState==0)+(initState==1) ), 1 )),...
                    'Initial state matrix should be all 0''s and 1''s')
                assert((dim(1)==numCells && dim(2)==numGenes),...
                    'Initial state matrix should have size (numCells x numGenes)')
                
                initialStates = initState;
            end
            obj.initStates = initialStates;
            
            % Generate the intercellular connectivity
            if isempty(initTruth)
                [outCells, inCells] = genInOutCells(obj);
                Ttable = genTtable(obj, inCells);
            elseif strcmp(topology, 'ising')
                dim = size(initTruth);
                inCells = 1;
                outCells = 1;
                Ttable = initTruth;               
            else
                dim = size(initTruth);
                inCells = 1:obj.bandwidth;
                outCells = (obj.numGenes - obj.bandwidth +1):obj.numGenes;
                
                assert(isempty(find(...
                    ~( (initTruth==-1)+(initTruth==0)+(initTruth==1) ), 1 )),...
                    'Truth table matrix should be all -1''s, 0''s, and 1''s')
                assert((dim(1)==2^k && dim(2) == numGenes),...
                    'Truth table matrix should have size (2^k x numGenes)')
                Ttable = initTruth;
            end
            
            % Generate the intracellular connectivity
            if isempty(initVar)
                varF = genvarF(obj, inCells);
            else
                dim = size(initVar);
                assert(isempty(find(...
                    ~( (initVar>=-1).*(initVar<=numGenes) ), 1 )),...23wde
                    'Initial connectivity matrix (varF) should only connect to nodes within the cell')
                assert((dim(1)==k && dim(2) == numGenes),...
                    'Connectivity matrix should have size (k x numGenes)')                
                varF = initVar;
            end
            
            % Set the properties just generated, which might change at some
            %   point in the simulation (thus they are just the initial values
            obj.initTtable    = Ttable;
            obj.initvarF      = varF;
            obj.initOutCells  = outCells;
            obj.initInCells   = inCells;
            
            % Create the cell objects that will be placed in this object's
            %   list (this object is a grid)
            allCells = cell(numCells,1);
            for cellPos=1:numCells
                allCells{cellPos}=boolCell_v2(numGenes,k,p,bandwidth,Ttable,varF,outCells,inCells,perturb);             
                allCells{cellPos}.setPos(cellPos);
                allCells{cellPos}.setState(initialStates(cellPos,:),1);
            end
            obj.allCells = allCells;          
            
            %Get and set the list of all neighbors
            obj.setNeighbors(numCells, topology);
            
            %Set initial time to 1
            obj.timenow = 1;
            
            %Output Mean-Field Approx Critcality of individual Boolean Network
            obj.crit_val = 2*k*p*(1-p);
                
            if (obj.crit_val < 1)
                obj.criticality = 'Ordered';
            elseif (obj.crit_val > 1)
                obj.criticality = 'Chaotic';
            else
                obj.criticality = 'Critical';
            end
            
            obj.movie_end_steps = 2986;    
            %obj.page_rank_default = 1;
            %obj.page_rank_const = .95;
             
            obj.temperature = temperature; % needs to change
            obj.J = J;
            obj.h_interac = h_interac;
            obj.E_magnetization = 0.0;
            h_f_temp = 1;
       
            if (strcmp(topology,'ising'))
                obj.E_disordered = zeros(numCells,1);
                neigh = obj.neighbors;
                for i = 1:obj.numCells
                    N_S_W_E = [obj.allCells{neigh(i,1)}.states(obj.bandwidth,obj.timenow), obj.allCells{neigh(i,2)}.states(obj.bandwidth,obj.timenow), ...
                             obj.allCells{neigh(i,3)}.states(obj.bandwidth,obj.timenow), obj.allCells{neigh(i,4)}.states(obj.bandwidth,obj.timenow)];
                    energy_i = energy(obj.allCells{i}.states(1,1), N_S_W_E, J, h_interac, h_f_temp);
                    
                    obj.E_disordered(i) = energy_i;
                end
                
                intercellular_nodes = zeros(obj.numCells, 1);
                for jCell=1:numCells
                    intercellular_nodes(jCell) = obj.allCells{jCell}.states(obj.bandwidth, obj.timenow);
                end
                
                figure(1);
                subplot(1,2,1);
                imagesc(reshape( intercellular_nodes, [sqrt(obj.numCells), sqrt(obj.numCells)])), title('Intercellular Node State (Spin)'), colorbar, colormap winter;
                subplot(1,2,2);
                imagesc(reshape(obj.E_disordered, [sqrt(obj.numCells), sqrt(obj.numCells)])), title('Initial Energy'), colorbar, colormap winter;
            end     
            
        end   
        
        %---------------------------------------------
        % Full update function
        %---------------------------------------------
        function obj = update_all(obj,numSteps)
            %Steps the simulation forward a given number of steps      
            %page_rank_default = 1;
            page_rank_const = 1; %page_rank_const .95
            tstart = obj.timenow;          
            for jT = tstart+1:(tstart+numSteps)  

                % First update the intracellular dynamics
                for jCell=1:obj.numCells
                    %page_rank_p = .95;
                    thisCell = obj.allCells{jCell};
                    if ( rand <= page_rank_const )
                        thisCell.update_genes(jT);
                    else
                        thisCell.page_rank(jT);
                    end                       
                end                
                % Then update the intercellular communication dynamics
                obj.update_intercell(jT); 
                  % If pertrubation exists, add perturbation
                if (obj.perturb > 0)
                    % Add Perturbation
                    for jCell=1:obj.numCells
                        thisCell = obj.allCells{jCell};
                        thisCell.perturb_genes(jT);
                    end
                else
                    % Do nothing
                end

            end           
            obj.timenow = jT;          
            obj.get_states;
            
        end
        
        %---------------------------------------------
        % Plotting function
        %---------------------------------------------
        function plot_cells(obj, save, delta ,dt)
            % Plots the cell states
            %   Only 'lines' implemented so far
            % Change the function so you can select range of states
            
            if nargin == 1
                save = false;
                dt = 0.0;
            end
            
            if (save == true)
                % This needs to be fixed
                % End timestep has to be same as jT end timestep
                mov(1:5000) = struct('cdata', [], 'colormap', []);
            end      
            
            isLine = strcmp(obj.topology,'line');
            
            if isLine
                for jT = 1:delta:size(obj.allStates,3)
                    
                    % Plot
                    imagesc(obj.allStates(:,:,jT).');
                    colorbar
                    colormap(hot);
                    title(sprintf('State at Step %d',jT));
                    drawnow;
                    pause(dt);
                    if (save == true)
                        mov(jT) = getframe(gcf);
                    end
                end
                v=VideoWriter('RBN.avi');
                v.FrameRate = dt;
                open(v)
                writeVideo(v, mov)
                close(v)  
            else
                %cesaro instantiation uncomment to use it
                b = zeros(sqrt(obj.numCells),sqrt(obj.numCells));
                movie_end_steps = obj.movie_end_steps;
                for jT = 1:delta:movie_end_steps; %size(obj.allStates,3)
                    
                    % We want to have a matrix of all the states, because
                    %   our actual cells are on a grid
                    stateVec = obj.allStates(:,:,jT); %Get the vector of states
                    stateVecFlat = zeros(obj.numCells,1);
                    for jGene=1:obj.numGenes
                        % Translate the gene state into a big integer
                        stateVecFlat = stateVecFlat + stateVec(:,jGene)*2^(jGene-1);
                    end
                    uniqueVec = unique(stateVecFlat);          
                    % Now get the matrix form
                    gridSide = sqrt(obj.numCells);
                    stateMat = reshape(stateVecFlat,[gridSide,gridSide]);
                    
                    %temp tbd (Cesaro Summation - uncomment to use it) 
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     b = b + stateMat;
%                     c = b/(jT);
%                     stateMat = c;
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    % Plot
                    imagesc(stateMat,[1 2^obj.numGenes]); 
                    c= colorbar;
                    %colormap(winter(length(uniqueVec)));
                    colormap(winter(2^obj.numGenes));
                    title(sprintf('State at Step %d',jT))
                    c.Label.String = 'Gene States';
                    drawnow
                    pause(dt)
                    if (save == true)
                        mov(jT) = getframe(gcf);
                    end
                end
            end
            if (save == true)
                v=VideoWriter('RBN.avi');
                v.FrameRate = 1/dt;
                open(v)
                writeVideo(v, mov)
                close(v)       
            end
        end
        
        %---------------------------------------------
        % Hamming distance
        %---------------------------------------------
        function rhs = ham(A,B)
            % Hamming distance is the number of positions at which A and B differ. Since
            %   A and B are matrices of zero and 1 this is just the sum of the absolute value of the difference
            %   between each entry in the final time step
            
            % The caller might have passed the boolCellGrid, not just the
            %   matrices
            if isa( A,'boolCellGrid_v2' )
                A = A.allStates;
            end
            if isa( B,'boolCellGrid_v2' ) 
                B = B.allStates;
            end
            
            % Check to make sure we have the same size matrices
            szA = size(A);
            szB = size(B);
            assert( szA(1)==szB(1) && szA(2)==szB(2) , ...
                'The states need to be the same size');   
            rhs = sum(sum(abs(A(:,:,end)-B(:,:,end))));
        end
        
        %---------------------------------------------
        % Insert mutant cell(s)
        %---------------------------------------------
        function obj = insert_mutants(obj, mutCellGrid, mutPos)
            % This function takes a grid of (probably identical) cells and
            %   adds a cell of the type inside 'mutCellGrid' in position(s) 'mutPos'
            
            % Safety checks
            assert( isa(obj,'boolCellGrid_v2'),...
                'First argument should be an object of class boolCellGrid_v2');
            assert( isa(mutCellGrid,'boolCellGrid_v2'),...
                'Second argument should be an object of class boolCellGrid_v2');
            assert( isvector(mutPos), ...
                'Third argument should be a vector with length >= 1');
                       
            % Overwrite the cells at mutPos with COPIES of mutCell (i.e. we
            %   can't have just a 'handle' pass by reference class)
            for j=1:length(mutPos)
                thisPos = mutPos(j);
                obj.allCells{thisPos} = ...
                    copy(mutCellGrid.allCells{thisPos});
            end           
        end
        
        %---------------------------------------------
        % Find Steady State Distribution of Each Cell
        %---------------------------------------------
        function rhs = ssDist(A)
            % Check/pass in a matrix, allStates
            if isa( A,'boolCellGrid_v2' )
                A = A.allStates;
            end
            
            % Return a Steady State Distributions
            numberCells = length(A(:,1,1));
            numberGenes = length(A(1,:,1));
            numberT = length(A(1,1,:));
            decMatrix = zeros(numberCells, numberT);
            
            % Convert states to decimal representation
            for ii = 1:numberCells
                for tt = 1:numberT
                    tempDec = 0;
                    for jj = 1:numberGenes
                        tempDec = tempDec + A(ii, jj, tt)*2^(numberGenes-jj);
                    end
                    decMatrix(ii,tt) = tempDec;
                end
            end
            
            % Return the Steady State Probability Distribution
            counts = zeros(numberCells, 2^numberGenes);
            for i = 1:length(decMatrix(:,1))
                for j = 1:length(decMatrix(1,:))
                    counts(i, decMatrix(i,j)+1) = counts(i, decMatrix(i,j)+1) + 1;
                end
            end
            rhs = counts/(sum(counts(1,:)));
        end      
    end
        
    methods (Access=private)
        
        %---------------------------------------------
        % Set the neighbors
        %---------------------------------------------
        function obj = setNeighbors(obj, numCells,topology)
            %Returns a list where the ROWS are the linear indices of the
            %neighbors. This uses Periodic Boundary Conditions, so each point
            %has the same number of neighbors
            
            switch topology
                
                case 'ising'
                      gridSide = round(sqrt(numCells));                    
                    assert(gridSide-sqrt(numCells)<1e-4,...
                        'Only square grids are supported for now');                    
                    matSize(1:2) = [gridSide;gridSide];                    
                    neighList = zeros(matSize(1)*matSize(2),4);    
                    
                    for jY = 1:matSize(1)
                        for jX = 1:matSize(2)
                            
                            %Get the linear index of the current point
                            jLin = sub2ind(matSize,jY,jX);
                            
                            %x neighbors
                            if jY~=1
                                jYmin1 = sub2ind(matSize,jY-1,jX);
                            else
                                jYmin1 = sub2ind(matSize,jY-1+matSize(1),jX);
                            end
                            if jY~=matSize(1)
                                jYplu1 = sub2ind(matSize,jY+1,jX);
                            else
                                jYplu1 = sub2ind(matSize,jY+1-matSize(1),jX);
                            end
                            
                            %y neighbors
                            if jX~=1
                                jXmin1 = sub2ind(matSize,jY,jX-1);
                            else
                                jXmin1 = sub2ind(matSize,jY,jX-1+matSize(2));
                            end
                            if jX~=matSize(2)
                                jXplu1 = sub2ind(matSize,jY,jX+1);
                            else
                                jXplu1 = sub2ind(matSize,jY,jX+1-matSize(2));
                            end
                            %disp(jLin)
                            neighList(jLin,:) = [jYmin1, jYplu1, jXmin1, jXplu1];
                        end
                    end
                    
                    
                case 'symmetric'
                    
                    gridSide = round(sqrt(numCells));                    
                    assert(gridSide-sqrt(numCells)<1e-4,...
                        'Only square grids are supported for now');                    
                    matSize(1:2) = [gridSide;gridSide];                    
                    neighList = zeros(matSize(1)*matSize(2),4);                    
                    for jY = 1:matSize(1)
                        for jX = 1:matSize(2)
                            
                            %Get the linear index of the current point
                            jLin = sub2ind(matSize,jY,jX);
                            
                            %x neighbors
                            if jY~=1
                                jYmin1 = sub2ind(matSize,jY-1,jX);
                            else
                                jYmin1 = sub2ind(matSize,jY-1+matSize(1),jX);
                            end
                            if jY~=matSize(1)
                                jYplu1 = sub2ind(matSize,jY+1,jX);
                            else
                                jYplu1 = sub2ind(matSize,jY+1-matSize(1),jX);
                            end
                            
                            %y neighbors
                            if jX~=1
                                jXmin1 = sub2ind(matSize,jY,jX-1);
                            else
                                jXmin1 = sub2ind(matSize,jY,jX-1+matSize(2));
                            end
                            if jX~=matSize(2)
                                jXplu1 = sub2ind(matSize,jY,jX+1);
                            else
                                jXplu1 = sub2ind(matSize,jY,jX+1-matSize(2));
                            end
                            %disp(jLin)
                            neighList(jLin,:) = [jYmin1, jYplu1, jXmin1, jXplu1];
                        end
                    end
                   
                case 'line'
                    
                    neighList = zeros(numCells,2);                 
                    for jX = 1:numCells                       
                        %1-d neighbors
                        if jX~=1
                            jXmin1 = jX-1;
                        else
                            jXmin1 = numCells;
                        end
                        if jX~=numCells
                            jXplu1 = jX+1;
                        else
                            jXplu1 = 1;
                        end                     
                        neighList(jX,:) = [jXmin1, jXplu1];
                    end
                    
                case 'single'
                    
                    neighList = zeros(numCells,2);                 
                    for jX = 1:numCells                       
                        % This part is the same as Line, but still gets the
                        % right neighList
                        if jX~=1
                            jXmin1 = jX-1;
                        else
                            jXmin1 = numCells;
                        end
                        if jX~=numCells
                            jXplu1 = jX+1;
                        else
                            jXplu1 = 1;
                        end                     
                        neighList(jX,:) = [jXmin1, jXplu1];
                    end    
                    
                case 'orthogonal'
                    
                    gridSide = round(sqrt(numCells));
                    assert(gridSide-sqrt(numCells)<1e-4,...
                        'Only square grids are supported for now');                    
                    matSize(1:2) = [gridSide;gridSide];                   
                    neighList = zeros(matSize(1)*matSize(2),2);
                    for jY = 1:matSize(1)
                        for jX = 1:matSize(2)
                            
                            %Get the linear index of the current point
                            jLin = sub2ind(matSize,jY,jX);
                            
                            %x neighbors; only to the West
                            if jY~=1
                                jYmin1 = sub2ind(matSize,jY-1,jX);
                            else
                                jYmin1 = sub2ind(matSize,jY-1+matSize(1),jX);
                            end
                            
                            %y neighbors; only to the South
                            if jX~=1
                                jXmin1 = sub2ind(matSize,jY,jX-1);
                            else
                                jXmin1 = sub2ind(matSize,jY,jX-1+matSize(2));
                            end                         
                            neighList(jLin,:) = [jYmin1, jXmin1];
                        end
                    end
                otherwise
                    error('Your topology isn''t supported')
            end         
            %Set the object property
            obj.neighbors = neighList;           
        end
                
        %---------------------------------------------
        % Generate the truth table for the cells
        %---------------------------------------------
        function Ttable = genTtable(obj, inCells)
       
            %Create a random truth table, assuming that each cell has
            %the same number of connections coming in, k
            
            Ttable = zeros(2^obj.k,obj.numGenes);

            for i = 1:2^obj.k
                for j = 1:obj.numGenes
                    %x = rand;
                    if (obj.p > rand) 
                        Ttable(i,j) = 1;
                    end
                end
            end
            
            %Ttable = randi([0,1], 2^obj.k, obj.numGenes);

            %Get rid of the internal connectivity for the cells
            %receiving input from their neighbors
            Ttable(:,inCells) = -1*ones( 2^obj.k,length(inCells) );

            %Check to make sure there's at least one way to turn on any
            %given gene
            %   Note that the inCells are special and will always have
            %   full -1 Ttables
%             for jGene = length(inCells)+1:obj.numGenes
%                 if isempty(find(Ttable(:,jGene), 1))
%                     Ttable(randi((obj.k)^2),jGene) = 1;
%                 end
%             end
        end     
        
        %---------------------------------------------
        % Generate the output and input connections for the cells
        %---------------------------------------------
        function varF = genvarF(obj, inCells)

            varF = -1*ones(obj.k, obj.numGenes);
            for jGene = 1:obj.numGenes

                if ~ismember(jGene,inCells)
                    %Randomize the connectivity, giving each node k
                    %connections
                    connectList = randperm(obj.numGenes); %Randomize the genes
                    connectList = connectList(connectList~=jGene); %Get rid of recurrent connections

                    varF(:,jGene) = connectList(1:obj.k).';
                else
                    %Do nothing; leave it at all -1
                end
            end
        end      
        
        %---------------------------------------------
        % Generate the output and input connections for the cells
        %---------------------------------------------
        function initStates = genInit(obj, topology)    
            initStates = randi([0,1],obj.numCells,obj.numGenes);
%             if strcmp(topology, 'ising')
%                initStates(:,1) = reshape(checkerboard_2(sqrt(obj.numCells)), obj.numCells, 1);
%             end
        end

        %---------------------------------------------
        % Generate the output and input connections for the cells
        %---------------------------------------------
        function [outC, inC] = genInOutCells(obj)          
            %The connectivity is already randomized, so let's keep it
            %easy and have the first x nodes receive input and the last
            %x send output
            inC = 1:obj.bandwidth;
            outC = (obj.numGenes-obj.bandwidth+1):obj.numGenes;
        end

        %---------------------------------------------
        % Inter-cellular communication
        %---------------------------------------------
        function obj = update_intercell(obj, timestep)
            %This function updates just the intercellular communication
            %between the cells that are on a grid
            
            %Get the table of neighbors
            topol = obj.topology;
            neigh = obj.neighbors;
            
            if strcmp(topol,'single')
                % single cell
                thisCell = obj.allCells{1};
                thisIn = thisCell.inCells(1); 
                thisOut = thisCell.outCells(1);
                
                %% Single cell has fixed/persistent input
                %thisCell.states(thisIn, timestep) = ...
                %    0; %manually change this to 0 or 1
                %thisCell.states(thisOut, timestep-1);
                
                %% Single cell has sampled receing input
                thisCell.states(thisIn, timestep) =...
                    randsample(2,1,true, [0 1])-1;
            elseif strcmp(topol,'ising')
               
                % ising model test
                intercellular_nodes_no_override = zeros(obj.numCells, 1);  
                intercellular_nodes = zeros(obj.numCells, 1);  
                for jCell=1:obj.numCells
                    intercellular_nodes_no_override(jCell) = obj.allCells{jCell}.states(obj.bandwidth, timestep);               
                    obj.allCells{jCell}.states(obj.bandwidth, timestep) = obj.allCells{jCell}.states(obj.bandwidth, timestep-1);
                    intercellular_nodes(jCell) = obj.allCells{jCell}.states(obj.bandwidth, timestep);
                end
                
                temperature_init = obj.temperature;
                
%                 figure(2)
%                 imagesc(reshape(intercellular_nodes, [sqrt(obj.numCells), sqrt(obj.numCells)])), colorbar, colormap winter
%                 title_time = strcat('Before calculation should be 0, temperature: ', ' ', string( temperature_init ) );
%                 title(title_time)

                thisOut = obj.bandwidth; 
                thisIn = obj.bandwidth;

                mc = randi(obj.numCells);
                neigh = obj.neighbors(mc, :);
                
                h_const = obj.h_interac;
                h_f = intercellular_nodes_no_override(mc);
                
                N_S_W_E = [obj.allCells{neigh(1)}.states(thisOut,timestep-1), obj.allCells{neigh(2)}.states(thisOut,timestep-1), ...
                         obj.allCells{neigh(3)}.states(thisOut,timestep-1), obj.allCells{neigh(4)}.states(thisOut,timestep-1)];
                energy_old = energy(obj.allCells{mc}.states(thisOut,timestep-1), N_S_W_E, obj.J, h_const, h_f);

                flip =  obj.allCells{mc}.states(thisOut,timestep-1);
                
                if( flip == 0)
                    % flip from 0 to 1
                    flip = 1;
                else
                    % flip from 1 to 0
                    flip = -1;
                end
                
                energy_new = energy(flip, N_S_W_E, obj.J, h_const, h_f);
                dE = energy_new - energy_old;
                k_B = -1; % -1 signifies dE < 0 passed
       
                if (dE < 0)
                    %unflip 
                    if( flip == -1)
                        obj.allCells{mc}.states(thisOut,timestep) = 0;
                    else
                        obj.allCells{mc}.states(thisOut,timestep) = 1;
                    end
                else
                    k_B = exp(-(dE)/temperature_init);
                    if (rand() < k_B) 
                        %display('flipped')
                        if( flip == -1)
                            obj.allCells{mc}.states(thisOut,timestep) = 0;
                        else
                            obj.allCells{mc}.states(thisOut,timestep) = 1;
                        end
                    end
                end
                
                if (mod(timestep, 100) == 0)
                    intercellular_nodes = zeros(obj.numCells, 1);  
                    for jCell=1:obj.numCells
                        intercellular_nodes(jCell) = obj.allCells{jCell}.states(obj.bandwidth, timestep);
                    end

                    temperature_init = obj.temperature;
                    figure(3)
                    imagesc(reshape(intercellular_nodes, [sqrt(obj.numCells), sqrt(obj.numCells)])), colorbar, colormap winter
                    title_time = strcat('Ater calculation, temperature: ', ' ', string( temperature_init ) );
                    title(title_time) 
                    suptitle(strcat('t= ', string(timestep), ',  mc= ', string(mc), ',   dE= ', string(dE), ',   k_B= ', string(k_B)));
                end
                  
            else 
                for jCell = 1:obj.numCells
                    %Get the cell object we'll update
                    %   Note: pass by REFERENCE
                    thisCell = obj.allCells{jCell};

                    %Update the nodes that receive extracellular communication, as
                    %passed by the caller inside the table 'interCell'
                    for jIn = 1:length(thisCell.inCells)

                        thisIn = thisCell.inCells(jIn); %The node receiving input
                        thisOut = thisCell.outCells(jIn); %The node OF THE NEIGHBOR that ouputs

                        neighStates = zeros(size(neigh,2),1);
                        % rows = size(a,1); cols = size(a,2)
                        % for i=1: in range rows
                        % end
                        neigh;
                        for jNeigh = 1:size(neigh,2)
                            %Get the index of the neighbor we want to query
                            thisNeigh = neigh(jCell,jNeigh);
                            neighStates(jNeigh) = obj.allCells{thisNeigh}.states(thisOut,timestep);
                        end
                        %Apply an 'or' function to the states of all of those
                        %output nodes, and that is the state of the node we're
                        %updating

                        %thisCell.states(thisIn,timestep) = ...
                        %     double(logical(sum(neighStates)));
                         
                        % Apply Linear Threshold of 2 (TEMPORARY)
                        neighStates; timestep;
                        if ( strcmp(topol,'symmetric'))
                            thisCell.states(thisIn,timestep) = ...
                                double(logical(sum(neighStates) >= 2));
                            
                            % thisCell.states(thisIn,timestep) = 1;                           
                            % test_tbd = checkerboard_2(sqrt(obj.numCells));
                            % thisCell.states(thisIn,timestep) = test_tbd(thisCell.cellPos);
                            
                        elseif strcmp(topol,'line')
                            thisCell.states(thisIn,timestep) = ...
                                double(logical(sum(neighStates) >= 1));
                        else
                            error('Your topology isn''t supported')
                        end
                            
                    end
                end         
            end           
        end      
        
        %---------------------------------------------
        % Get the states out from the single cell objects
        %---------------------------------------------
        function get_states(obj)
            obj.allStates = zeros(obj.numCells,obj.numGenes,obj.timenow);
            for jCell = 1:obj.numCells
                obj.allStates(jCell,:,:) = obj.allCells{jCell}.states;
            end
        end
    end    
end

