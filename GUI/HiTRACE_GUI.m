% Copyright (c) Rhiju Das and Sungroh Yoon 
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
% * Redistributions of source code must retain the above copyright notice, 
%   this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright 
%   notice, this list of conditions and the following disclaimer in the 
%   documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clear all;

global loadingTimer;
global stepNo;
global SETTING;
global DATA;
global setter;
global res;
stepNo = -1;

if(~isdeployed)
    addpath(genpath('../Scripts'));
    addpath(genpath('etc/'));
end

while(1)
    if(stepNo < 0)
        %%%% Environment variable initialize

        res = findResource;
        setter = SettingWindow;
        uiwait(setter.output);

        
	if(~isempty(SETTING))
        %%% Viewer Figure
            uihandles = Viewer;
            v = uihandles.output;        
	    
            if(~isdeployed) 
	      % is this necessary?
	      %matlabpool('close', 'force');
            end

            numOfJob = SETTING.numOfJob;
            if matlabpool( 'size' ) == 0 && numOfJob ~= 0;
                if(res.ClusterSize < numOfJob)
                    matlabpool( res.ClusterSize );
                else
                    matlabpool( numOfJob );
                end
            end

            loadingTimer = 0;
            delete(uihandles.axes3);
            pause(0.01);
        end
    end 
    
    if(~isempty(SETTING))
        refcol = SETTING.refcol;
        ymin = SETTING.ymin;
        ymax = SETTING.ymax;
        sequence = SETTING.sequence;
        offset = SETTING.offset;
        filenames = SETTING.list;
	primer_distance_from_end = SETTING.primer_distance_from_end;
        
        eBaseline = SETTING.enableBaseline;
        eRefine = SETTING.enableRefine;

        mutposFile = SETTING.mutposFile;
        %marksFile = SETTING.marksFile;
        
        enableDP = SETTING.enableDP;
        slack = SETTING.slack;
        maxShift = SETTING.maxShift;
        windowSize = SETTING.windowSize;
        res = findResource;
        
        SETTING.enableBaseline;
        SETTING.enableRefine;

        SETTING.mutposFile;
%        SETTING.marksFile;

        if(~isempty(SETTING.targetBlock))
	  
	    %str = sprintf('targetBlock = [%s];', SETTING.targetBlock);
            %str = strrep(str, '-', ':');
            %eval(str);

	    % actually we need targetBlock to be a cell...	    
	    s = SETTING.targetBlock;
	    s_string = 'targetBlock = {';
	    r = s;
	    while( ~isempty(r) )
	      [t,r] = strtok( r, ';' );
	      s_string = [ s_string,  '[',strrep( t, '-', ':' ),']' ];
	      if ~isempty(r); s_string = [ s_string, ',' ]; end;
	    end
	    s_string = [s_string, '}' ];
	    eval( s_string );
	    
	end
    else
        break;
    end
    
    timeline = tic;

    PLOT_STUFF = 0;	  
    
    if(stepNo < 1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP A (PREPROCESSING) + STEP B (CORRELATION OPTIMIZED LINEAR ALIGNMENT)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %--------------------------------------------------------------------------
        % Steps A., B.1 & B.2: Preprocessing, Intra- and Inter- batch alignment
        %--------------------------------------------------------------------------
        % In Step A, signal windowing and basline adjustment is performed as 
        % preprocessing. In step B, HiTRACE aligns preprocessed profiles by a 
        % linear shifting and scaling of the time axis, maximizing the correlation 
        % between the reference and each profile. Alignment quality improves 
        % progressively. Step B.1 is for intra-batch alignment, and Step B.2 is for
        % inter-batch alignment. 
        %--------------------------------------------------------------------------
               
        [ d, da] = quick_look( filenames, ymin, ymax, [], [], refcol, PLOT_STUFF );
        DATA.d = d;
        DATA.da = da;
        if length( d ) == 0
	  errordlg('Could not find data?','Error!');
	  return; % Sungroh I need your help with this! Graceful exit!
	else	  
	  uihandles.step = 1;
	  uihandles.d = d;
	  uihandles.da = da;
	  guidata(v, uihandles);
	  
	  figure(v);
	  
	  h1 = axes('position', [0.1, 0.2, 0.35, 0.7]);
	  image( 50 * d);
	  axis( [ 0.5 size(da,2)+0.5 ymin ymax] );
	  title( 'Signal')
	  
	  h2 = axes('position', [0.55, 0.2, 0.35, 0.7]);
	  image( 50 * da);
	  axis( [ 0.5 size(da,2)+0.5 ymin ymax] );
	  title( 'Reference ladder')
	  
	  colormap( 1- gray(100));
	  
	  set(uihandles.stepLabel, 'String', 'Loading and alignment (step 1/5)');
	  
	  pause(0.01);
	end
    end
    
    if(stepNo < 2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % baseline subtract removes offset
        if(eBaseline)
	  d_bsub = baseline_subtract_v2( d, ymin, ymax, 2e6, 2e4, PLOT_STUFF );
        else
            d_bsub = d;
        end
        
        DATA.d_bsub = d_bsub;
        
        uihandles.step = 2;
        uihandles.d_bsub = d_bsub;
        guidata(v, uihandles);
        
        figure(v);
        
        if(stepNo ~= 1)
            delete(h1);
            delete(h2);
        end
        
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);
        
        if(eBaseline)
            str = 'Signal (baseline subtraction enabled)';
        else
            str = 'Signal (baseline subtraction skipped)';
        end
        
        image( 50 * d_bsub);
        axis( [ 0.5 size(d_bsub,2)+0.5 ymin ymax] );        
        title( str );
        
        colormap( 1- gray(100));
        
        set(uihandles.stepLabel, 'String', 'Baseline subtraction (step 2/5)');
        
        pause(0.01);
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP C (DP-BASED NONLINEAR ADJUSTMENTS) + STEP D (PEAK FITTING)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %--------------------------------------------------------------------------
        % Step C.1: Nonlinear alignment
        %--------------------------------------------------------------------------
        % With current CE equipment, we found that it was not feasible to get 
        % complete alignment for profiles with better than single-band resolution 
        % with just a linear scaling of the time axis. To correct for these issues, 
        % we perform another round of refining the alignment to handle this 
        % nonlinearity after completion of the linear alignment procedures. This 
        % step handles the nonlinearities in electrophoretic mobility that arise 
        % due to small differences in temperature, geometry, and loading between 
        % different capillaries. This step breaks the time axis of a profile into 
        % small windows and utilizes a dynamic programming algorithm to align these 
        % separate segments.
        %--------------------------------------------------------------------------

    if(stepNo < 3)
        if SETTING.enableDP

	  if(~exist('targetBlock')) targetBlock = []; end;
	  if( length( targetBlock ) == 0 ); targetBlock = [ 1: size( d_bsub, 2 ) ]; end;
	  if ~iscell( targetBlock ) 
	    targetBlock = { targetBlock };
	  end
	  
	  % need to put in a warning here ... if ddATP ladders, etc., are in align_by_DP step they might be messed up.
	  if(~isempty(mutposFile)) 
            data_type = textread(mutposFile,'%s'); % this probably should only occur once...
	    showed_error = 0;
	    for n = 1:length( targetBlock )
	      for j = 1:length( targetBlock{n} )
		m = targetBlock{n}(j);
		if ( m <= length(data_type) &  ~isempty( strfind( data_type{m}, 'dd' ) ) )
		  showed_error = 1;
		  errordlg('You might be applying dynamic programming fine alignment to a lane with a sequencing ladder, and could get weird results','OK'); % Sungroh, can we have an option "change align blocks" that returns to the setting window?
		  break;
		end
	      end
	      if showed_error; break; end;
	    end
	  end

	  	  
	  d_align= align_by_DP( d_bsub(ymin:ymax, : ), targetBlock, slack, maxShift, windowSize, PLOT_STUFF);
	  DATA.d_align = d_align;
      
      
	  str  = 'Signal (nonlinear alignment by DP)';
        
	else
	  d_align = d_bsub(ymin:ymax, : );
	  DATA.d_align = d_align;
            str = 'Signal (nonlinear alignment by DP skipped)';
        end
        
        uihandles.step = 3;
        uihandles.d_align = d_align;
        guidata(v, uihandles);
        
        figure(v);
        
        if(stepNo ~= 2)
            delete(h1);
        end
                
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);
        
        
        image( 50 * d_align);
        axis( [ 0.5 size(d_bsub,2)+0.5 ymin ymax] );        
        title( str );
        
        colormap( 1- gray(100));
        
        set(uihandles.stepLabel, 'String', 'Nonlinear alignment by DP (step 3/5)');
        pause(0.01);
    end
        %--------------------------------------------------------------------------
        % Steps C.2 & D: Automated transfer of band annotation and deconvolution
        %--------------------------------------------------------------------------
        % In the band annotation procedure (step C.2), HiTRACE automatically 
        % assigns the bands on each profile by transferring annotation of a single 
        % reference profile by using a second dynamic programming algorithm. This 
        % step not only increases user convenience but also helps improve alignment 
        % quality and facilitates the last deconvolution process (step D), which 
        % fits peaks to predefined peak models to quantify final band intensities.
        %--------------------------------------------------------------------------
        % FITTING = 1;
        % if FITTING
        %     [ Fitted_Area ] = do_fitting( d_SHAPE_shift_realign_norm, numDataset, numRep );
        % end


        % Locations of mutations in MedLoop
    
    if(stepNo < 4)
        set(uihandles.stepLabel, 'String', 'Manual annotation (step 4/5)');
        
        if(stepNo ~= 3)
            delete(h1);
        end
        
        figure(v);
        
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);

        if ~exist( 'data_type' ) data_type = {}; end;
        % these marks will show up to guide assignments, they are pairs (i,j).
        marks = [];
	if(~isempty(mutposFile)) 
            data_type = textread(mutposFile,'%s');
	end	    
	
        %if(~isempty(marksFile))
        %    marks = load(marksFile);
        %end
	
	
	seqpos = ( (length(sequence)-primer_distance_from_end) : -1 :1 ) + offset;
	if ~exist( 'structure' );  % later fix this!
	  structure = sequence;
	  for m = 1:length(sequence); structure(m) = '.';end;
	end
	[ marks, area_pred, mutpos] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence, offset, seqpos, data_type);	
	marks = [ marks; [1:length(sequence)]'+offset, [1:length(sequence)]'+offset ];
	
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % manual annotation

	%xsel = [];
	% this should let us refine xsel if we change our mind on assignments... but not working. Sungroh, help!
	if exist( 'uihandles.xsel' ) xsel = uihandles.xsel;
	else xsel = [];  end;
	
	if length(xsel) == 0; cla; end;

	period = 1; USE_GUI = 1;
        xsel = mark_sequence( d_align, xsel, sequence(1:end-primer_distance_from_end), 0, offset, period, marks, mutpos, USE_GUI);

        numpeaks = length(xsel);
                
        uihandles.step = 4;
        uihandles.xsel = xsel;
        uihandles.numpeaks = numpeaks;
        uihandles.mutpos = mutpos;
        uihandles.marks = marks;
        guidata(v, uihandles);
        pause(0.01);
    end
    
    if(stepNo < 5)
        % do the fit.
        
        if(eRefine)
            [area_peak, prof_fit] = do_the_fit_fast( d_align, xsel', 0.0, PLOT_STUFF );            
            str = 'Fitting without refinement (step 5/5)';
        else
            [area_peak, prof_fit] = do_the_fit_GUI( d_align, xsel' );
            str = 'Fitting with refinement (step 5/5)';
        end
        
        
        if(stepNo ~= 4)
            delete(h1);
        end

	xsel_start = [];
	if length( xsel_start ) > 0; xsel_start = xsel(1,:); end;
        
        uihandles.step = 5;
        uihandles.area_peak = area_peak;
        uihandles.prof_fit = prof_fit;
        uihandles.xsel_start = xsel_start;
        guidata(v, uihandles);
        
        figure(v);

	h1 = axes('Position', [0.08, 0.2, 0.27, 0.7]);
	scalefactor = 40/mean(mean(d_align));
	image( scalefactor * d_align );
	if length( xsel_start ) > 0; ylim( [ min(xsel_start)-100, max(xsel_start)+100 ] ); end;
	title( 'data' );
	  
	h2 = axes('Position', [0.38, 0.2, 0.27, 0.7]);
	image( scalefactor * prof_fit );
	if length( xsel_start ) > 0; ylim( [ min(xsel_start)-100, max(xsel_start)+100 ] ); end;
	title( 'fit' );
	  
	h3 = axes('Position', [0.68, 0.2, 0.27, 0.7]);
	image( scalefactor * ( d_align - prof_fit) );
	if length( xsel_start ) > 0; ylim( [ min(xsel_start)-100, max(xsel_start)+100 ] ); end;
	title( 'residuals' );
	
        colormap( 1 - gray(100) );
        
        set(uihandles.stepLabel, 'String', str);
        
        uihandles.axesHandles = [h1 h2 h3];
        guidata(v, uihandles);
        
        set(uihandles.prevBtn, 'Enable', 'On');
        set(uihandles.nextBtn, 'Enable', 'On');
        set(uihandles.restartBtn, 'Enable', 'On');
        pause(0.01);
    end
    
    fin=toc(timeline);
    str = sprintf('Finished (running time = %f seconds). In order to return to the initial setting window, please close the viewer window.', fin);
    set(setter.saveBtn, 'Enable', 'On');
    setter.area_peak = area_peak;
    guidata(setter.output, setter);

    if ( exist( 'area_peak' ) & length( area_peak ) > 0 )
      h = questdlg(str, 'Finish!', 'Done', 'Save data', 'Save data');
      if(strcmp(h,'Save data'))
        [name path] = uiputfile('Output.rdat', 'Save to RDAT file');
        if(name)
	  fullname = strcat(path, name);
	  seqpos = length(sequence)-primer_distance_from_end - [0:(size(area_peak,1)-1)] + offset;
	  rdat = fill_rdat(name, sequence, offset, seqpos,area_peak);
	  output_rdat_to_file(fullname,rdat);
        end
      end
      uiwait(v);
    end
    set(setter.saveBtn, 'Enable', 'Off');
end