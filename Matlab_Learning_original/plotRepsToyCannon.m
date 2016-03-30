clear all, close all

                           
filenames = {
'RepsToyCannon_NewSamp20_eps100.mat'                                     
'RepsToyCannon_NewSamp20_eps50.mat'
'RepsToyCannon_NewSamp20_eps75.mat'                                      
'RepsToyCannon_NewSamp40_eps100.mat'                                     
'RepsToyCannon_NewSamp40_eps50.mat'                                      
'RepsToyCannon_NewSamp40_eps75.mat'                                      
'RepsToyCannon_NewSamp60_eps50.mat'                                      
'RepsToyCannon_NewSamp60_eps75.mat'
'RepsToyCannon_NewSamp60_eps100.mat'};

for i = 1:length(filenames)

    load(filenames{i})
    figure, 
    plot(median(R.mean, 1))
    title(filenames(i))
end
    