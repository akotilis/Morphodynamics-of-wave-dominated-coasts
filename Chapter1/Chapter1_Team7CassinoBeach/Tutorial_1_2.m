%% 1.2.1 Prelimanary visualisations
clear all;
close all;

lowtide = load('lowTide.txt');
hightide = load('highTide.txt');

f_s =2;
duration = length(lowtide)/f_s;
t = linspace(1/f_s,duration,length(lowtide));

P1_lt = lowtide([1:end], 1);
P3_lt = lowtide([1:end], 2);
P6_lt = lowtide([1:end], 5);

P1_ht = hightide([1:end], 1);
P3_ht = hightide([1:end], 2);
P6_ht = hightide([1:end], 5);

%% Plots

figure(1); 

subplot(3,2,1); plot(t, P1_lt); xlim([0,duration]);ylim([-2 2]); xlabel('time(s)'); ylabel('h(m)'); title('P1 lowtide');
subplot(3,2,3); plot(t, P3_lt); xlim([0,duration]);ylim([-2 2]); xlabel('time(s)'); ylabel('h(m)'); title('P3 lowtide');
subplot(3,2,5); plot(t, P6_lt); xlim([0,duration]);ylim([-2 2]); xlabel('time(s)'); ylabel('h(m)'); title('P6 lowtide');

subplot(3,2,2); plot(t, P1_ht); xlim([0,duration]);ylim([-2 2]); xlabel('time(s)'); ylabel('h(m)'); title('P1 hightide');
subplot(3,2,4); plot(t, P3_ht); xlim([0,duration]);ylim([-2 2]); xlabel('time(s)'); ylabel('h(m)'); title('P3 hightide');
subplot(3,2,6); plot(t, P6_ht); xlim([0,duration]);ylim([-2 2]); xlabel('time(s)'); ylabel('h(m)'); title('P6 hightide');

%% 1.2.2 Computation of wave statistics

hi_test = [31 43 45 23 3 34 66 24 13];

significant_height(hi_test)


