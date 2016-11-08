clear all, close all

javaaddpath('../poi_library/poi-3.8-20120326.jar');
javaaddpath('../poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('../poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('../poi_library/xmlbeans-2.3.0.jar');
javaaddpath('../poi_library/dom4j-1.6.1.jar');
javaaddpath('../poi_library/stax-api-1.0.1.jar');
addpath ../

[num, txt, raw] = xlsread('RobotMeasurements.xlsx');

jani = dlmread('Batch1ForceMeasurement.txt');

figure, subplot(3,1,1), plot(jani(:, 1)/1e3, jani(:, 2:4)), title('position'), legend('X', 'Y', 'Z')
subplot(3,1,2), plot(jani(:, 1)/1e3, (jani(:, 5:7))), title('Forces'), legend('X', 'Y', 'Z')
subplot(3,1,3), plot(jani(:, 1)/1e3, (jani(:, 8:10))), title('troque'), legend('X', 'Y', 'Z')


feri = dlmread('Batch1ParametersLog.txt');
context = [context; feri(:, 1)];
