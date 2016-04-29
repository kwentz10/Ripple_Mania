%created by Katherine Wentz 3/7/16
%Tipple Generation

%% INITIALIZE

clearvars 
clear figure(1)

%-----Model Constants

mu=0.35;
alpha_min=10;
alpha_range=3;
D=0.25; %(mm) grain diameter
B=1; %fraction of impacting grains with high enough energy to play a significant role in the ejection of other grains
S=10; %impact rate number of impacts per mm2 s2
Nej=12; %number of ejecta per impact

%-----Initial Arrays

%Space Array
x_min=0; %(mm) 
x_max=1000; %(mm) 
dx=20*D; %(mm) bin spacing
spacearray=x_min:dx:x_max;%(mm)
j_max=length(spacearray); % number of bins

%Impact Array
i_max=10000; %number of impacts
impactarray=1:1:i_max;

%Elevation Array
zr_mean=50; %(mm) inital mean elevation of ripple surface 
zr_amp=1;
elevationarray=zr_amp*rand(1,j_max)+zr_mean; %(mm)

%Variable Elevation Array
elevationarray_int=elevationarray;

%Number of Grains in Each Bin Array
Narray=elevationarray./(pi*(D^2)).*(4*(1-mu)*dx);

%Array of Random Impact Elevations Above Ripple Elevation Surface
max_impactelev=x_max*tand(alpha_min);
impactelevarray=randi(floor(max_impactelev),1,i_max);

%Array of Random Impact Angles 
impactalpha=randi(alpha_range,1,i_max)+alpha_min;

%Time Array
tmin=1;
tmax=i_max/(x_max*D*B*S); 

%-----Initial Plot

figure(1)
plot(spacearray, elevationarray_int,'g','linewidth',2)
title('Ripple Evolution')
xlabel('Distance (mm)','fontname','arial','fontsize',21)
ylabel('Elevation (mm)','fontname','arial','fontsize',21)
ylim([0.9*zr_mean zr_mean+100])
set(gca,'fontsize',18,'fontname','arial')
legend('Ripples')
pause(0.05)



%% RUN 

nplots = 100; %number of plots
iplot = i_max/nplots; %defines which iterations I want to plot

for i=1:i_max
    
    %determine maximum and minimum heights of impact elevation
    min_impact_elev=impactelevarray(1); %minimum impact elevation set as first elevation in space array
    max_impact_elev=max(elevationarray+(spacearray*tand(impactalpha(i)))); %apply alpha to entire space array to find the greatest change in height along space (also add to elevation array to find total height)
    
    %determine random impact elevation using min max impact elevations
    impact_elev=min_impact_elev+((max_impact_elev-min_impact_elev)*rand); %choose a random impact elevation from minimum and maximum impact elevations based on impact angle (could hit anywhere along space array at that angle), but what elevation is below or above the spacearray(1) elevation? Might screw this up--need residuals?
   
    %offset horizontal impact_elev by trajectory line  
    trajectory=impact_elev-tand(impactalpha(i))*spacearray; %find the elevation where each space array value will hit at spacearray(1)
   
    
    below=find(trajectory<elevationarray); %determine where trajectory is lower than elevations in space
    above=1:below(1)-1; %indices 1 through the first point before trajectory is lower than elevation array
    impact_bin=below(1)-1; %find point where trajectory is lower than elevation array
    
    if impact_bin>j_max    %wrap around
        impact_bin=impact_bin-j_max;
    end
    
    %splash function
    Narray(impact_bin)=Narray(impact_bin)-Nej;
    bins=[];
    bins=impact_bin+3:impact_bin+6; %the bins in which grains land
    
    for ii=1:length(bins)  %wrap around
        if bins(ii)>j_max
            bins(ii)=bins(ii)-j_max;
        end
    end
   
    Narray(bins(1)) = Narray(bins(1)) + ceil(0.25*Nej);
    Narray(bins(2)) = Narray(bins(2)) + ceil(0.5*Nej);    
    Narray(bins(3)) = Narray(bins(3)) + ceil(0.125*Nej);    
    Narray(bins(4)) = Narray(bins(4)) + ceil(0.125*Nej);
    
    %new elevation of ripples
    elevationarray=(pi.*Narray*D^2)/(4*(1-mu)*dx); 
    
    %plot
    if rem(impactarray(i),iplot)==0
        
        clear figure(1)
        %impact trajectory
        plot(spacearray(above), trajectory(above),'--b','linewidth',2)
        hold on
        %ripple elevation
        plot(spacearray, elevationarray,'g','linewidth',2)
        title('Ripple Evolution')
        xlabel('Distance (mm)','fontname','arial','fontsize',21)
        ylabel('Elevation (mm)','fontname','arial','fontsize',21)
        ylim([0.9*zr_mean zr_mean+100])
        set(gca,'fontsize',18,'fontname','arial')
        legend('Impact Trajectory','Ripples')
        pause(0.05)
        hold off
      
    end
    
end


%% FINALIZE

figure(1)
plot(spacearray, elevationarray_int,'g','linewidth',2)
hold on
plot(spacearray,elevationarray,'k','linewidth',2)
title('Ripple Evolution')
xlabel('Distance (mm)','fontname','arial','fontsize',21)
ylabel('Elevation (mm)','fontname','arial','fontsize',21)
ylim([0.9*zr_mean zr_mean+100])
set(gca,'fontsize',18,'fontname','arial')
legend('Initial Elevation','Ripple Generation')
text(200,500,[tmax,'s'],'fontsize',14);






