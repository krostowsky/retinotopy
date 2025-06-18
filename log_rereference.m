function log_rereference(outdir, outfile, step, inText)
if ~exist([outdir '/' outfile]) || step == 0
    fclose(fopen([outdir '/' outfile] ,'w+'));
end
fid = fopen([outdir '/' outfile],'a');

switch step
    case 0 % log number of intracranial contacts & their CSC file
        fprintf(fid, 'Number of intracranial contacts:\n');
        fprintf(fid, '%s\n', inText{1});
        fprintf(fid, 'CSC files for intracranial contacts:\n');
        for i = 1:size(inText{2}, 1)
            fprintf(fid, '%s\n', sprintf('%s %s\n', inText{2}{i,1}, sprintf('%d ', inText{2}{i,2})));
        end
    case 1 % log sorting of CSC files (should be the same as before, this is just to ensure correct order)
        fprintf(fid, 'CSC files for intracranial contacts (sorted):\n');
        for i = 1:size(inText{1}, 1)
            fprintf(fid, '%s\n', sprintf('%s %s\n', inText{1}{i,1}, sprintf('%d ', inText{1}{i,2})));
        end
    case 2 % log if corrupt
        fprintf(fid, 'no corrupt channels found\n');    
    case 3
        fprintf(fid, 'corrupt channels found, removed: %s\n', inText);
end

fclose(fid);

end