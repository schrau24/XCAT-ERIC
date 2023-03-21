function [varargout] = bart_evalc(cmd, varargin)
% BART	Call BART command from Matlab.
%   [A, B] = bart('command', X, Y) call command with inputs X Y and outputs A B
%
% 2014-2016 Martin Uecker <uecker@med.uni-goettingen.de>


	if nargin==0 || all(cmd==0)
		disp('Usage: bart <command> <arguments...>');
		return
	end

	bart_path = getenv('TOOLBOX_PATH');

	if isempty(bart_path)
		if exist('/usr/local/bin/bart', 'file')
			bart_path = '/usr/local/bin';
		elseif exist('/usr/bin/bart', 'file')
			bart_path = '/usr/bin';
		else
			error('Environment variable TOOLBOX_PATH is not set.');
		end
	end

	% clear the LD_LIBRARY_PATH environment variable (to work around
	% a bug in Matlab).

	if ismac==1
		setenv('DYLD_LIBRARY_PATH', '');
	else
		setenv('LD_LIBRARY_PATH', '');
	end

	name = tempname;

	in = cell(1, nargin - 1);

	for i=1:nargin - 1,
		in{i} = strcat(name, 'in', num2str(i));
		writecfl(in{i}, varargin{i});
	end

	in_str = sprintf(' %s', in{:});

	out = cell(1, nargout);

	for i=1:nargout-1,
		out{i} = strcat(name, 'out', num2str(i));
	end

	out_str = sprintf(' %s', out{:});

	if ispc
		% For cygwin use bash and modify paths
        bartstr = ['bash.exe --login -c ', ...
			strrep(bart_path, filesep, '/'), ...
	                '"', '/bart ', strrep(cmd, filesep, '/'), ' ', ...
			strrep(in_str, filesep, '/'), ...
                	' ', strrep(out_str, filesep, '/'), '"'];
		T = evalc('ERR = system(bartstr);');
    else
        bartstr = ['nice -n19 ', bart_path, '/bart ', cmd, ' ', in_str, ' ', out_str];
		T = evalc('ERR = system(bartstr);');
	end

	for i=1:nargin - 1,
		if (exist(strcat(in{i}, '.cfl'),'file'))
			delete(strcat(in{i}, '.cfl'));
		end

		if (exist(strcat(in{i}, '.hdr'),'file'))
			delete(strcat(in{i}, '.hdr'));
		end
	end

	for i=1:nargout,
        if i == 1
            varargout{i} = T;
        else
            if ERR==0
                varargout{i} = readcfl(out{i-1});
            end
            if (exist(strcat(out{i-1}, '.cfl'),'file'))
                delete(strcat(out{i-1}, '.cfl'));
            end
            if (exist(strcat(out{i-1}, '.hdr'),'file'))
                delete(strcat(out{i-1}, '.hdr'));
            end
        end
	end

	if ERR~=0
		error('command exited with an error');
	end
end
