function allpics(dirname,extname)

files=dir([dirname '/*.' extname]);

for i=1:size(files),
    disp(['Processing file ' dirname '/' files(i).name])
    makepic([dirname '/' files(i).name]);
end