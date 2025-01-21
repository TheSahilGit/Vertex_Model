



% Open the file for reading
filename = 'para1_in.dat'; % Replace with your actual file name
fid = fopen(filename, 'r');

% Initialize an array to store the first column
firstColumn = [];

% Read the file line by line
while ~feof(fid)
    line = fgetl(fid); % Read a single line
    
    % Use regular expressions to extract the first column
    tokens = regexp(line, '^\s*(\d+)', 'tokens');
    
    if ~isempty(tokens)
        firstColumn(end+1) = str2double(tokens{1}{1}); % Convert to a number and store
    end
end

% Close the file
fclose(fid);

% Display the first column
disp('First Column:');
disp(firstColumn);
