function [data, headers]=read_data(dataset)
%http://lib.stat.cmu.edu/jasadata/
%http://www.stat.columbia.edu/~gelman/book/data/
%http://lib.stat.cmu.edu/datasets/

switch dataset
    case 'light'
        headerlines = 4;
        filename = 'light.asc';
        fid = fopen(filename, 'r');
        [s, pos]=textscan(fid,  '%d', 'HeaderLines', headerlines);
        fclose(fid);
        data = double(s{1}');
        headers = {'time'};
    otherwise
        error('unknown data set');
end
