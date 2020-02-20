function code_for_figure_s1()

%%%%%%%%%%%%%%%%%%%%%
%FIGURES S1E AND S1F%
%%%%%%%%%%%%%%%%%%%%%

load('airFlowSensor_versus_cannula_test.mat')
figure

%Histograms of instantaneous breathing rates (Figure S1E)
subplot(2,2,1)
hold on
histogram(BR_exh_airflow)
histogram(BR_exh_cannula)
ranksum(BR_exh_airflow,BR_exh_cannula)
subplot(2,2,3)
hold on
histogram(BR_inh_airflow)
histogram(BR_inh_cannula)
ranksum(BR_inh_airflow,BR_inh_cannula)

%Shift between cannula and air flow sensor (Figure S1F)
subplot(2,2,[2 4])
hold on
histogram(Shift_exh,[-10:2:30]/1000)
histogram(Shift_inh,[-10:2:30]/1000)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
%FIGURE S1B%
%%%%%%%%%%%%

%PID test of the facial mask, example (Figure S1B)
load('pid_test_facialMask.mat')
nbOdors = 16;
o = 1;
figure
for s = 1:2
    f = find(odorlist(:,1)==o & odorlist(:,2)==s);
    subplot(1,2,s)
    hold on
    plot(PIDtracesR(:,f),'r');
    ylim([-0.1 0.5])
    
    f = find(odorlist(:,1)==o & odorlist(:,2)~=s);
    subplot(1,2,s)
    hold on
    plot(PIDtracesL(:,f),'b');
    ylim([-0.1 0.5])
end