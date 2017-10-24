function [data_pop] = getPopulationLevelData()

%from Shook et al. 2012
data_pop.age = [60 90 120 180 270 360 720].*24; %hours
data_pop.S = [6200 5805 NaN 4607 NaN 1617 607];
data_pop.S_sem = [200 519.5 NaN 369 NaN 170.1 109.6];

%from Daynac et al.2016:
%now QS and AS in our models refer to QS and AS from Daynac:
data_pop.QS = [347.5 NaN 296.8 322.2 267.9 248.0 NaN];
data_pop.QS_sem = [112.2 NaN 94.1 81.4 54.3 126.7 NaN];

data_pop.AS = [307 NaN 241.6 352.2 276.2 305.2 NaN];
data_pop.AS_sem = [47.2 NaN 119.9 43.6 45.4 81.8 NaN];   

data_pop.TAP = [784.1 NaN 480.4 353.4 224.5 171.2 NaN];
data_pop.NB = [19573.2 NaN 9611.9 7522.3 5575.2 5014.4 NaN];
data_pop.NB1 = [5850 NaN 2750 1800 1750 800 NaN]./2;
data_pop.NB2 = data_pop.NB-data_pop.NB1;

data_pop.TAP_sem = [47.9 NaN 99.4 29.4 38.7 27.6 NaN];
data_pop.NB_sem = [972.3 NaN 852.5 522.8 772.0 952.8 NaN];
data_pop.NB1_sem = [250 NaN 150 120 80 40 NaN]./2;
data_pop.NB2_sem = data_pop.NB_sem-data_pop.NB1_sem;

end