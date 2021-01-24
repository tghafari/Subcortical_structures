% Finds potential error outputs of FIRST 

filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/Processed_Data/';
allSubs=35;
errorLogs=cell(allSubs,1);
for nSub=1:allSubs
    cd(filePath)
    if numel(num2str(nSub))==1; Subject=['S0' num2str(nSub)]; else; Subject=['S' num2str(nSub)]; end %Create subject name
    if isfolder([filePath Subject '.anat']) %check if fsl_anat has been run on the subject's MRI
        cd([Subject '.anat/first_results'])
        cmd='cat *.logs/*.e*';
        system(cmd);
        disp(Subject)
    else
        warning('no .anat folder')
    end
end
 
