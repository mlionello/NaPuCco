classdef MPrint < handle
%        a = MPrint();
%        a.println("hello there")
%           Hello there
%        a.print("I am working on %d out of %d steps\r", 1, 5); 
%           Hello there
%           I am working on 1 out of 5 steps
%        a.print("I am working on %d out of %d steps\r", 2, 5); 
%           Hello there
%           I am working on 2 out of 5 steps
%        a.println("I had enough of all this")
%           Hello there
%           I had enough of all this
%        a.print("I am working on ")
%        a.print("%d out of %d steps\r", 1, 5); 
%           Hello there
%           I had enough of all this
%           I am working on 1 out of 5 steps
%        a.print("%d out of %d steps\r", 2, 5); 
%           Hello there
%           I had enough of all this
%           I am working on 2 out of 5 steps     
%       a.flush(); a.print("zero")
%           Hello there
%           I had enough of all this
%           I am working on 2 out of 5 steps
%           zero

    properties
        carriagereturn;
        strlen;
        newlinewhenflushed;
        logfile;
        tmp_cr;
    end

    methods
        function obj = MPrint(logfilename, opts)
            arguments
                logfilename = nan;
                opts.newlinewhenflushed = 1;
            end
            obj.carriagereturn = "";
            obj.newlinewhenflushed = opts.newlinewhenflushed;
            obj.logfile = logfilename;
            obj.tmp_cr=0;
        end

        function print(obj, varargin)
            tolog = 1;
            cr = 0;
            message = string(varargin{1});
            if strlength(message)>2 && endsWith(message, "\r")
                cr = 1;
                message = extractBefore(message, "\r");
            end
            if string(varargin{end})=="tolog"
                tolog = 1;
                varargin = varargin(1: end-1);
            end
            if length(varargin) > 1
                A = varargin(2: end);
                for i = 1:  length(A)
                    message = compose(message, A{i});
                end
            end
            obj.strlen = strlength(message);
            difflenght = max(strlength(obj.carriagereturn)/2-obj.strlen,0);
            if ~endsWith(message, "\n")
                fprintf(obj.carriagereturn + message + ...
                    repmat(' ', 1, difflenght));
                obj.carriagereturn = "";
                if cr
                    obj.carriagereturn = repmat('\b', 1, obj.strlen + difflenght);
                end
            else
                fprintf(obj.carriagereturn + message);
                obj.carriagereturn = "";
                if cr
                    obj.carriagereturn = repmat('\b', 1, obj.strlen);
                end
            end
            try
                if tolog
                    fid = fopen(obj.logfile, 'a');
                    fprintf(fid, message + "\n");
                    fclose(fid);
                end
            catch 
                if isnan(obj.logfile)
                    warning('no logfile specified')
                else
                    warning("not possible to record log to file")
                end
            end
        end

        function tolog(obj, varargin)
            message = varargin(1);
            message = string(message{1});
            if length(varargin) > 1
                A = varargin(2: end);
                for i = 1:  length(A)
                    message = compose(message, A{i});
                end
            end
            fid = fopen(obj.logfile, 'a');
            fprintf(fid, strjoin(message) + "\n");
            fclose(fid);
        end

        function flush(obj)
            obj.carriagereturn = "";
            if obj.newlinewhenflushed
                fprintf(obj.carriagereturn + "\n")
            else
                fprintf(obj.carriagereturn)
            end
        end

        function clear(obj)
            if obj.tmp_cr>0
                fprintf(repmat('\b', 1, obj.tmp_cr));
                obj.tmp_cr = 0;
            end
        end

        function printcr_wait(obj, varargin)
            fmt = varargin;
            message = string(fmt{1});
            if strlength(message)>2 && endsWith(message, "\r")
                message = extractBefore(message, "\r");
            end
            if string(fmt{end})=="tolog"
                fmt = fmt{1: end-1};
            end
            if length(fmt) > 1
                A = fmt(2: end);
                for i = 1:  length(A)
                    message = compose(message, A{i});
                end
            end
            obj.tmp_cr = strlength(message);
            obj.print(varargin{:});
        end


        function println(obj, varargin)
            varargin{1} = varargin{1} + "\n";
            obj.print(varargin{:});
        end
    end
end