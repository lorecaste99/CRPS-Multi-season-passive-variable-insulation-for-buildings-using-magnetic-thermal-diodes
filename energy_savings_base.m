%%
clc
clear all
close all

%acquiring temperature data every 30 minutes
Data=readtable('C:\Users\caste\Box\laptop\Documents\PhD\Research\IBUILD\SGTD\US map heat flow\data files\los_angeles_CA_2022.csv',Range="A1:K17524");

%Saving the variables from the data sheets
T_hourly=Data{3:end,6}; %temperature
DNI_hourly=Data{3:end,8}; %Direct normal irradiance
dew_point_hourly=Data{3:end,9}; %dew point
cloud_type_hourly=Data{3:end,7}; %cloud type (integer value)
zenith_hourly=Data{3:end,10}; %zenith angle
zenith_hourly(zenith_hourly>90)=90; %the formulas do not allow for a zenith angle above 90 degrees
alpha_hourly=90-zenith_hourly; %solar altitude angle

%taking longitude and latitude from data sheet
longitude=-Data{1,7};
latitude=Data{1,6};

%calculating the declination and equation of time for every day of the year
N=1:365;
declination=23.45.*sind(360.*(284+N)./365);
tau=360*N/365;
ET=-7.3412*sind(tau)+0.4944*cosd(tau)-9.3795*sind(2*tau)-3.2568*cosd(2*tau)-0.3179*sind(3*tau)-0.0774*cosd(3*tau)-0.1739*sind(4*tau)-0.1283*cosd(4*tau);

%creating time vectors
x = (datetime(2022,1,1,00,00,00):seconds(1):datetime(2022,12,31,23,30,00))';
x.Format = 'MMM dd, HH:mm';
date_30_mins = (datetime(2022,1,1,00,00,00):minutes(30):datetime(2022,12,31,23,30,00))';

%repeating ET and declination elements to match the time vectors
ET=(repelem(ET,48))';
declination=(repelem(declination,48))';

%calculating the incidence angle
solar_time= (date_30_mins) + minutes(4)*(-longitude) + minutes(ET);
omega=(hours(timeofday(solar_time))).*15-180;
solar_azimuth=real(asind((cosd(declination)).*sind(omega)./cosd(alpha_hourly)));
cos_incidence=cosd(alpha_hourly).*cosd(-solar_azimuth);

%calculating the incoming radiaton hitting the south facing vertical wall
rad_hourly=DNI_hourly.*cos_incidence;

%calculating average monthly temperature for each month
T_january=mean(T_hourly(1:1488,1));
T_february=mean(T_hourly(1488:2832,1));
T_march=mean(T_hourly(2832:4320,1));
T_april=mean(T_hourly(4320:5760,1));
T_may=mean(T_hourly(5760:7248,1));
T_june=mean(T_hourly(7248:8688,1));
T_july=mean(T_hourly(8688:10176,1));
T_august=mean(T_hourly(10176:11664,1));
T_september=mean(T_hourly(11664:13104,1));
T_october=mean(T_hourly(13104:14592,1));
T_november=mean(T_hourly(14592:16032,1));
T_december=mean(T_hourly(16032:17520,1));

%repeating the elements to have the right vector size
T_january=repelem(T_january,2678400);
T_february=repelem(T_february,2419200);
T_march=repelem(T_march,2678400);
T_april=repelem(T_april,2592000);
T_may=repelem(T_may,2678400);
T_june=repelem(T_june,2592000);
T_july=repelem(T_july,2678400);
T_august=repelem(T_august,2678400);
T_september=repelem(T_september,2592000);
T_october=repelem(T_october,2678400);
T_november=repelem(T_november,2592000);
T_december=repelem(T_december,2678400-1799);

%storing the average temperatures for reference
T_year=[T_january T_february T_march T_april T_may T_june T_july T_august T_september T_october T_november T_december];

%interpolating variables from NREL data to get a time resolution of 1 s
T_ambient=interp1(date_30_mins,T_hourly,x);
dew_point=interp1(date_30_mins,dew_point_hourly,x)+273.15; %this has to be in kelvin for later formulas
cloud_type=interp1(date_30_mins,cloud_type_hourly,x);
absorbance=0.6;
q_rad=interp1(date_30_mins,rad_hourly,x)*absorbance; %I am already accounting for absorption

%calculating sky temperature using formulas provided in SI
e_correction=(1+0.024.*(cloud_type)-0.0035.*(cloud_type.^2)+0.00028.*(cloud_type.^3));
T_sky=(((0.787+0.764.*log(dew_point./273)).*e_correction.*((T_ambient+273.15).^4)).^0.25); %this is in kelvin

%PVI circuit parameters for base case
R_value=10; 
RSI_value_off=R_value/5.678; %convert to units of m^2*K/W
switch_ratio=2.5;
RSI_value_on=RSI_value_off/switch_ratio;

%if optimized scenario transition temperatures desired --> uncomment
% T_max_1=24;
% T_min_1=21;
% T_max_2=22;
% T_min_2=19;

%if base case scenario transition temperatures desired --> uncomment
T_max_1=25.6;
T_min_1=21.2;
T_max_2=22.1;
T_min_2=17.7;

%building setpoint temperatures for base case scenario
R_rev=RSI_value_off;
T_summer=24;
T_winter=20;

%heat transfer parameters
h_ext=32;
h_int=8.29;
R_int=1/h_int;
emittance=0.8;
F_sky=sqrt(0.5)*0.5;
F_out=1-F_sky;
T_sky_c=T_sky-273;

%initial conditions for code
T_wall_out(1)=T_winter;
T_wall_in(1)=T_winter;
T_b(1)=T_winter;
h_rad_out(1)=3;
h_rad_sky(1)=1.5;

T_wall_out_no_circ(1)=T_winter;
T_wall_in_no_circ(1)=T_winter;
h_rad_out_no_circ(1)=3;
h_rad_sky_no_circ(1)=1.5;
%%
for i=2:1:size(T_ambient)

    %circuit
    if T_year(i)>=20;
        T_building=T_summer; %set this to summer scenario is T_month > 20
    else
        T_building=T_winter; %set this to summer scenario is T_month < 20
    end

    T_b(i)=T_building; %storing building temperature for plots later

    %calculating the h values for radiation
    h_rad_out(i)=F_out*emittance*(5.67e-8)*(((T_ambient(i)+273.15)^2)+(T_wall_out(i-1)+273.15)^2)*(T_ambient(i)+273.15+T_wall_out(i-1)+273.15);
    h_rad_sky(i)=F_sky*emittance*(5.67e-8)*((T_sky(i)^2)+(T_wall_out(i-1)+273.15)^2)*(T_sky(i)+T_wall_out(i-1)+273.15);

    %setting the argument for the function "circuit_R" at the bottom of the code
    diode1_arg=[T_wall_out(i-1),T_wall_in(i-1),T_max_1,T_min_1,T_max_2,T_min_2,RSI_value_off,switch_ratio];
    R_wall=circuit_R(diode1_arg);

    %calculating the inner and outer wall temperatures
    T_wall_out(i)=(T_ambient(i)*(h_ext+h_rad_out(i))+q_rad(i)+h_rad_sky(i)*T_sky_c(i)+T_building/R_wall/(1+R_wall/R_int))/(1/R_wall-1/(1+R_wall/R_int)+h_ext+h_rad_sky(i)+h_rad_out(i));
    T_wall_in(i)=(T_wall_out(i)/(1+R_wall/R_int))+T_building/(1+R_int/R_wall);

    %calculating the q entering or leaving building
    q_gain_loss(i)=(T_wall_out(i)-T_wall_in(i))/R_wall; %2 ways to compute this but they should be equal to each other

    %separating this in winter or summer months
    if T_b(i)==T_winter
        q_winter(i)=q_gain_loss(i);
        q_summer(i)=0;
    else
        q_winter(i)=0;
        q_summer(i)=q_gain_loss(i);
    end

   %everything is repeated for the NO PVI scenario
   h_rad_out_no_circ(i)=F_out*emittance*(5.67e-8)*(((T_ambient(i)+273.15)^2)+(T_wall_out_no_circ(i-1)+273.15)^2)*(T_ambient(i)+273.15+T_wall_out_no_circ(i-1)+273.15);
   h_rad_sky_no_circ(i)=F_sky*emittance*(5.67e-8)*((T_sky(i)^2)+(T_wall_out_no_circ(i-1)+273.15)^2)*(T_sky(i)+T_wall_out_no_circ(i-1)+273.15);

   T_wall_out_no_circ(i)=(T_ambient(i)*(h_ext+h_rad_out_no_circ(i))+q_rad(i)+h_rad_sky_no_circ(i)*T_sky_c(i)+T_building/R_rev/(1+R_rev/R_int))/(1/R_rev-1/(1+R_rev/R_int)+h_ext+h_rad_sky_no_circ(i)+h_rad_out_no_circ(i));
   T_wall_in_no_circ(i)=(T_wall_out(i)/(1+R_rev/R_int))+T_building/(1+R_int/R_rev);
   q_gain_loss_no_circ(i)=(T_wall_out_no_circ(i)-T_wall_in_no_circ(i))/R_rev;

    if T_b(i)==T_winter
        q_winter_no_circ(i)=q_gain_loss_no_circ(i);
        q_summer_no_circ(i)=0;
    else
        q_winter_no_circ(i)=0;
        q_summer_no_circ(i)=q_gain_loss_no_circ(i);
    end
end
%%

%calculating the cumulative energy savings
summer_q_integral=sum((q_summer))
winter_q_integral=sum((q_winter))
summer_q_integral_static=sum((q_summer_no_circ))
winter_q_integral_static=sum((q_winter_no_circ))

dQ=abs(q_gain_loss-q_gain_loss_no_circ);
tot_dQ=sum(dQ)


%calculating the HVAC load reduction in percentage
sum_summer=sum((q_summer));
sum_winter=sum((q_winter));

if sum_summer>=0
    sum_summer=sum_summer;
else 
    sum_summer=0;
end

if sum_winter<=0
    sum_winter=sum_winter;
else 
    sum_winter=0;
end

sum_summer_no_circ=sum(q_summer_no_circ);
sum_winter_no_circ=sum((q_winter_no_circ));

if sum_summer_no_circ>=0
    sum_summer_no_circ=sum_summer_no_circ;
else 
    sum_summer_no_circ=0;
end

if sum_winter_no_circ<=0
    sum_winter_no_circ=sum_winter_no_circ;
else 
    sum_winter_no_circ=0;
end

q_tot=abs(abs(sum_summer)+abs(sum_winter));
% q_tot=abs(abs(sum(q_winter)+abs(sum(q_summer));
q_tot_no_circ=abs(sum_summer_no_circ)+abs(sum_winter_no_circ);
%q_tot_no_circ=abs(sum((q_summer_no_circ)))+abs(sum((q_winter_no_circ)));
dQ_new=q_tot-q_tot_no_circ;
perc_inc_dQ=dQ_new/q_tot_no_circ*100

%%
%plots
close all

%change time interval to interval of interest
t1=1;
t2=2678400;

%creating the sum vectors
static_sum = zeros(1, t2-t1);
pvi_sum = zeros(1, t2-t1);

static_sum(1) = q_gain_loss_no_circ(t1);
pvi_sum(1) = q_gain_loss(t1);

%recurring sum of the cumulative energy entering or leaving the building
for i = 2:(t2-t1)
    static_sum(i) = static_sum(i-1) + q_gain_loss_no_circ(t1+i);
    pvi_sum(i) = pvi_sum(i-1) + q_gain_loss(t1+i);
end

%heat flow versus time graph
figure(1)
plot(x(t1:t2),q_gain_loss(t1:t2),'r','LineWidth',2)
hold on
plot(x(t1:t2),q_gain_loss_no_circ(t1:t2),'b','LineWidth',2)
set(gca, 'FontName','cambria math','FontSize',20)
xlim([datetime(2022,1,1,00,00,00) datetime(2022,1,31,00,00,00)])

%temperature versus time graph
figure(2)
plot(x(t1:t2),T_ambient(t1:t2),'Color','b','LineWidth',2)
hold on
plot(x(t1:t2),T_b(t1:t2),'Color','g','LineWidth',2)
hold on
plot(x(t1:t2),T_wall_out(t1:t2),'Color','black','LineWidth',2)
set(gca, 'FontName','cambria math','FontSize',20)
xlim([datetime(2022,1,1,00,00,00) datetime(2022,1,31,00,00,00)])

%cumulative energy versus time
figure(3)
plot(x(t1:t2-1),pvi_sum,'r','LineWidth',2)
hold on
plot(x(t1:t2-1),static_sum,'b','LineWidth',2)
set(gca, 'FontName','cambria math','FontSize',20)
%xlim([datetime(2022,7,1,00,00,00) datetime(2022,7,31,00,00,00)])
xlim([datetime(2022,1,1,00,00,00) datetime(2022,1,31,00,00,00)])

%this is the PVI argument function which calculates the PVI resistance as a
%function of transition temperatures and temperature conditions
function R_circuit=circuit_R(diode_arg)
%diode_arg=T_1,T_2,T_min_top,T_max_bottom,G_on_function

T_ambient=diode_arg(1);
T_building=diode_arg(2);
T_max_1=diode_arg(3);
T_min_1=diode_arg(4);
T_max_2=diode_arg(5);
T_min_2=diode_arg(6);
RSI_value_off=diode_arg(7);
switch_ratio=diode_arg(8);

RSI_value_on=RSI_value_off/switch_ratio;

if T_ambient>=T_max_1 && T_building<=T_min_1
    R_circuit=RSI_value_on;
elseif T_building>=T_max_2 && T_ambient<=T_min_2
    R_circuit=RSI_value_on;
else
    R_circuit=RSI_value_off;
end
end

