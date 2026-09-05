function v3_write_json(file,value)
[folder,~,~]=fileparts(file);if ~isempty(folder)&&~isfolder(folder),mkdir(folder);end
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0);cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
