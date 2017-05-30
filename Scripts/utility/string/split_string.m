function cols = split_string( l, delimiter );

cols = strsplit( l, delimiter );
return;

% if ~exist('delimiter', 'var'); delimiter = ' '; end;
% delimiter = sprintf(delimiter);
% 
% remain = l;
% cols = {};
% while length( remain ) > 0
%   [token, remain] = strtok(remain, delimiter);
%   cols = [cols, token];
% end
% 

%in_space = 1;
%count = 0;
%for i = [1:length(l)]
%  if ( strcmp( l(i), delimiter ) )
%    in_space = 1;
%  else 
%    if ( in_space )
%      count = count + 1;
%      cols{ count }  = '';
%      in_space = 0;
%    end
%    cols{count} = [cols{count}, l(i)];
%  end
%end

