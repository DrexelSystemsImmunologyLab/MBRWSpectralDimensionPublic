function [basefilename,extension] = rmfileextension(filename)
%RMFILEEXTENSION Separate the file extension from the rest of the file name.

extensionstart = regexp(filename,'[.][^.]+$');
if isempty(extensionstart)
    basefilename = filename;
    extension = '';
else
    basefilename = filename( 1:(extensionstart-1) );
    extension = filename(extensionstart:end);
end

end

