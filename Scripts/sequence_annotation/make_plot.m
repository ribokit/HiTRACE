function  annotation_handles = make_plot( d_align, xsel, sequence, offset, area_pred, annotation_handles, font_size )

% an 'annotation handle' consists of the text, horizontal line, and circle markers
%  that go with each sequence position.
%
% each handle will have { xsel position, handle for line, handle for text, handle for circle1, ... };
%
% the annotation handles are indexed by their sequence position [without using offset], but
% going backward from the very first one [at length(Sequence) ]. This is to allow
% additional bands that are 5' to the nominal sequence beginning.
%

numlanes = size( d_align, 2 ) ;

seqpos = get_seqpos( sequence, offset, xsel );

% why not just set contrast_factor based on colormap?
% might be much faster.

xlim = get(gca,'xlim');
xscale = abs(xlim(1) - xlim(2));
hold on

checked_handle = zeros( length( annotation_handles ),1 );
for i = length( xsel ):-1:1
    
    seq_idx = seqpos(i) - offset;
    handle_idx = length( sequence ) - seq_idx + 1;
    checked_handle( handle_idx ) = 1;
    if handle_idx <= length( annotation_handles ) && ~isempty( annotation_handles{handle_idx} )
        
        % This will move up or down handles that have already been created.
        handles = annotation_handles{ handle_idx };
        
        if ( handles{1} == xsel(i) ) continue; end;
        %fprintf( 'Updating handle: %d\n', handle_idx );
        handles{1} = xsel(i);
        annotation_handles{handle_idx} = handles;
        
        for m = 2:length( handles )
            h = handles{m};
            if isprop( h, 'YData' )
                newposition = get(h,'Ydata');
                newposition = newposition * 0 + xsel(i);
                set( h, 'YData', newposition );
            elseif isprop( h, 'Position' )
                newposition =  get( h, 'Position' );
                newposition( 2 ) = xsel(i);
                set( h, 'Position', newposition );
            end
        end
        
    else % need to create the handle.
        
        %fprintf( 'New handle: %d\n', handle_idx );
        
        handles = {};
        handles{1} = xsel(i);
        
        mycolor = [ 1 0 1]; % magenta for unknown.
        
        hold on
        h = plot( [0.5 numlanes+0.5 ], [xsel(i) xsel(i)], 'color',mycolor );
        handles = [ handles, h ];
        
        if ( seq_idx >=1  &&  seq_idx <= length(sequence) )
            
            % eterna colors. :)
            seqchar = sequence( seq_idx   );
            switch seqchar
                case {'A','a'}
                    mycolor = [0 0 1];
                case {'C','c'}
                    mycolor = [0 0.5 0];
                case {'U','T','u','t'}
                    mycolor = [1 0.5 0];
                case {'G','g'}
                    mycolor = [1 0 0];
            end
            set(h(1),'color',mycolor);
            
            txt_to_show = seqchar;
            seqnum = seqpos(i);
            txt_to_show = [seqchar,num2str( seqnum)];
            
            h = text( 0.5, xsel(i), txt_to_show );
            set(h,'HorizontalAlignment','right');
            set(h,'fontweight','bold','fontsize', font_size);%,'clipping','on');
            handles = [handles, h ];
            
            mark_points = find( area_pred(seq_idx,:) > 0.5 );
            h = plot( mark_points, mark_points*0 + xsel(i) , 'ro' );
            handles = [handles, h ];
            
        end
        
        annotation_handles{ handle_idx } = handles;
    end
    
end

% delete any handles that are no longer tagged to an 'xsel'
delete_handle_idx = find( ~checked_handle );
for i = 1:length( delete_handle_idx )
    handle_idx = delete_handle_idx(i);
    handles = annotation_handles{ handle_idx };
    for m = 2:length( handles)
        %set( handles{m},'visible','off' );
        delete( handles{m} );
    end
    annotation_handles{ handle_idx } = {};
end

axis off
hold off;