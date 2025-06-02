function log_epoch(outdir, outfile, step, inText)
if ~exist([outdir '/' outfile]) || step == 0
    fclose(fopen([outdir '/' outfile] ,'w+'));
end
fid = fopen([outdir '/' outfile],'a');

switch step
    case 0 % log CSC file directory 
        fprintf(fid, 'Directory used for CSC files:\n');
        fprintf(fid, '%s\n', inText{1});
        fprintf(fid, 'Behavioral files used:\n');
        fprintf(fid, '%s\n', inText{2}{:}); 
    case 1 % log jacksheet contacts to be analyzed
        fprintf(fid, ['Electrodes analyzed: ' num2str(length(inText)) '\n']);
        fprintf(fid, '%s\n', inText{:});
    case 2 % log CSC file used for TTL
        fprintf(fid, 'CSC File used for TTL:\n');    
        fprintf(fid, '%s\n', inText);
    case 3

end

fclose(fid);


end