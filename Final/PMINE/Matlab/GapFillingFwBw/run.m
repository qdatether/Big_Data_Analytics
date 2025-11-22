dir_load = load_directory;
dir_save = save_directory;
load(dir_load)
segmentsCor = zeros(size(segmentsArt));

noteM = [1 0.22]; % code = 1 for TMS, 0.22 s it appears
tbeforeTMS = 0.01;
tafterTMS = 0.03;
taus = 1/1450;
tbeforerec = 0.21; % 0.2 and 0.24 in old
tafterrec = 0.20;
memb = 100;
nnei = 2;
lowpasspar = 0.01; % was = 0.3
highpasspar = 100; % was 40
firorder = 65;  % was = 60
samplefreq = 1/taus;

for i = 1:5
    xM = segmentsArt(:,:,i)';
    xnewM = fillgapsnnfbAnn(xM,noteM,tbeforeTMS,tafterTMS,taus,tbeforerec,tafterrec,memb,nnei);
    segmentsCor(:,:,i) = xnewM';
    segmentsCorFilt100(:,:,i) = eegfilt(xnewM',samplefreq,lowpasspar,highpasspar,0,firorder);    
end

save(dir_save,'segmentsCor','segmentsCorFilt100');