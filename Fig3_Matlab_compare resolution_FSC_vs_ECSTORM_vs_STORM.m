[storm_table,A,B]=make_2_STORM_image();

EC_STORM_1=imread('UV8000frameA.tif');
EC_STORM_2=imread('UV8000frameB.tif');

[FSC_STORM_raw_EC,frequency_STORM]=FSC_STORM(EC_STORM_1,EC_STORM_2);
EC_STORM_s=smooth(frequency_STORM,FSC_STORM_raw_EC, 0.2,'loess');

[FSC_STORM_raw_UV,frequency_STORM]=FSC_STORM(EC_STORM_1,EC_STORM_2);
UV_STORM_s=smooth(frequency_STORM,FSC_STORM_raw_UV, 0.2,'loess');

threshold=ones(numel(frequency_STORM),1)*0.143;


plot(frequency_STORM,FSC_STORM_raw_EC,'kd');
hold on 
plot(frequency_STORM,EC_STORM_s,'b-');


plot(frequency_STORM,FSC_STORM_raw_UV,'ro');
 
plot(frequency_STORM,UV_STORM_s,'k-');
plot(frequency_STORM,threshold);
xlim([1e-5 0.4e-1]);
ylim([-0.1 1.1]);
xlabel('Spatial frequency (nm^-^1)');
ylabel('Fourier ring correlation');

legend('~','EC','~','UV');
hold off
