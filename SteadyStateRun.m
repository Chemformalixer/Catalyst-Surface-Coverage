function [Coverage_fraction_O, Coverage_fraction_CO, CO2_Production_rate]=SteadyStateRun(N,max_trials,yCO,CO2_time_step)

totalO2=zeros(1,max_trials);
totalCO=zeros(1,max_trials);
totalCO2=zeros(1,max_trials);
%CO2_Prod_rate=zeros(1,max_trials);

%S=struct('CO',0,'CO2',0,'O2',0);
%G=N*S;
CO_Grid_Occupancy=zeros(N+5,N+5);
O2_Grid_Occupancy=zeros(N+5,N+5);
%CO2_Grid_Occupancy=zeros(N,N);
%Grid_Occupancy.CO2=zero(N,N);

%initialization
yO2=1-yCO;
CO2_COUNT=0;
j=0;
for i_trial=1:max_trials
    % Picking a molecule
    j=j+1;
    if j==CO2_time_step
        CO2_COUNT=0;
        j=0;
    end
    Gas_Molecule_Pick=rand();
    if Gas_Molecule_Pick < yCO
        Cord_x=floor(rand()*N)+3;
        Cord_y=floor(rand()*N)+3;

        if CO_Grid_Occupancy(Cord_x,Cord_y)==0 && O2_Grid_Occupancy(Cord_x,Cord_y)==0
            CO_Grid_Occupancy(Cord_x,Cord_y)=1;
            rxnflag=0;
            %chkd_once1=0; chkd_once2=0; chkd_once3=0; chkd_once4=0;
            while rxnflag==0 && (O2_Grid_Occupancy(Cord_x-1,Cord_y)+O2_Grid_Occupancy(Cord_x+1,Cord_y)+O2_Grid_Occupancy(Cord_x,Cord_y-1)+O2_Grid_Occupancy(Cord_x,Cord_y+1))>0
                order=floor(rand()*4)+1;
                
                if order==1 && O2_Grid_Occupancy(Cord_x-1,Cord_y)==1 %seperate chkd_once1 from the other conditions
                    O2_Grid_Occupancy(Cord_x-1,Cord_y)=0;
                    CO_Grid_Occupancy(Cord_x,Cord_y)=0;
                    %chkd_once1=1;
                    rxnflag=1;
                elseif order==2 && O2_Grid_Occupancy(Cord_x+1,Cord_y)==1
                    O2_Grid_Occupancy(Cord_x+1,Cord_y)=0;
                    CO_Grid_Occupancy(Cord_x,Cord_y)=0;
                    %chkd_once2=1;
                    rxnflag=1;
                elseif order==3 && O2_Grid_Occupancy(Cord_x,Cord_y-1)==1
                    O2_Grid_Occupancy(Cord_x,Cord_y-1)=0;
                    CO_Grid_Occupancy(Cord_x,Cord_y)=0;
                    %chkd_once3=1;
                    rxnflag=1;
                elseif order==4 && O2_Grid_Occupancy(Cord_x,Cord_y+1)==1
                    O2_Grid_Occupancy(Cord_x,Cord_y+1)=0;
                    CO_Grid_Occupancy(Cord_x,Cord_y)=0;
                    %chkd_once4=1;
                    rxnflag=1;
                end
            end
            CO2_COUNT=rxnflag+CO2_COUNT;
        end
    elseif Gas_Molecule_Pick > yCO
        horizontal=0;
        vertical=0;
        x_y_flip=rand();
        if x_y_flip>0.5
            ad_x=floor(rand()*2);
            ad_y=0;
            horizontal=1;
            if ad_x==0
                ad_x=-1;
            end
        else
            ad_y=floor(rand()*2);
            ad_x=0;
            vertical=1;
            if ad_y==0
                ad_y=-1;
            end
            
        end
        
        
        x1=floor(rand()*N)+3;
        y1=floor(rand()*N)+3;

        x2=ad_x+x1;
        y2=ad_y+y1;
        
        Cord_x1=min(x1,x2);
        Cord_y1=min(y1,y2);
        Cord_x2=max(x1,x2);
        Cord_y2=max(y1,y2);
        
        if CO_Grid_Occupancy(Cord_x1,Cord_y1)==0 && O2_Grid_Occupancy(Cord_x1,Cord_y1)==0 && CO_Grid_Occupancy(Cord_x2,Cord_y2)==0 && O2_Grid_Occupancy(Cord_x2,Cord_y2)==0
            O2_Grid_Occupancy(Cord_x1,Cord_y1)=1;
            O2_Grid_Occupancy(Cord_x2,Cord_y2)=1;
            rxnflag1=0;
            rxnflag2=0;
            if horizontal==1
                while ((rxnflag1==0) && (CO_Grid_Occupancy(Cord_x1-1,Cord_y1)+CO_Grid_Occupancy(Cord_x1,Cord_y1-1)+CO_Grid_Occupancy(Cord_x1,Cord_y1+1)>0)) || ((rxnflag2==0) && (CO_Grid_Occupancy(Cord_x2,Cord_y2-1)+CO_Grid_Occupancy(Cord_x2+1,Cord_y2)+CO_Grid_Occupancy(Cord_x2,Cord_y2+1))>0)
                    order=floor(rand()*6)+1;
                    
                    
                    if order==1 && CO_Grid_Occupancy(Cord_x1-1,Cord_y1)==1 %seperate chkd_once1 from the other conditions
                        CO_Grid_Occupancy(Cord_x1-1,Cord_y1)=0;
                        O2_Grid_Occupancy(Cord_x1,Cord_y1)=0;
                        rxnflag1=1;
                    elseif order==2 && CO_Grid_Occupancy(Cord_x2+1,Cord_y2)==1
                        CO_Grid_Occupancy(Cord_x2+1,Cord_y1)=0;
                        O2_Grid_Occupancy(Cord_x2,Cord_y2)=0;
                        rxnflag2=1;
                    elseif order==3 && CO_Grid_Occupancy(Cord_x1,Cord_y1-1)==1
                        CO_Grid_Occupancy(Cord_x1,Cord_y1-1)=0;
                        O2_Grid_Occupancy(Cord_x1,Cord_y1)=0;
                        rxnflag1=1;
                    elseif order==4 && CO_Grid_Occupancy(Cord_x2,Cord_y2-1)==1
                        CO_Grid_Occupancy(Cord_x2,Cord_y1-1)=0;
                        O2_Grid_Occupancy(Cord_x2,Cord_y2)=0;
                        rxnflag2=1;
                    elseif order==5 && CO_Grid_Occupancy(Cord_x1,Cord_y1+1)==1
                        CO_Grid_Occupancy(Cord_x1,Cord_y1+1)=0;
                        O2_Grid_Occupancy(Cord_x1,Cord_y1)=0;
                        rxnflag=1;
                    elseif order==6 && CO_Grid_Occupancy(Cord_x2,Cord_y2+1)==1
                        CO_Grid_Occupancy(Cord_x2,Cord_y1+1)=0;
                        O2_Grid_Occupancy(Cord_x2,Cord_y2)=0;
                        rxnflag2=1;
                        
                    end
                end
                CO2_COUNT=CO2_COUNT+rxnflag2+rxnflag1;
                
            elseif vertical==1
                
                rxnflag1=0;
                rxnflag2=0;
                
                
                while ((rxnflag1==0) && (CO_Grid_Occupancy(Cord_x1,Cord_y1-1)+CO_Grid_Occupancy(Cord_x1-1,Cord_y1)+CO_Grid_Occupancy(Cord_x1+1,Cord_y1)>0)) || ((rxnflag2==0) && (CO_Grid_Occupancy(Cord_x2-1,Cord_y2)+CO_Grid_Occupancy(Cord_x2,Cord_y2+1)+CO_Grid_Occupancy(Cord_x2+1,Cord_y2))>0)
                    order=floor(rand()*6)+1;
                    
                    if order==1 && CO_Grid_Occupancy(Cord_x1,Cord_y1-1)==1
                        CO_Grid_Occupancy(Cord_x1,Cord_y1-1)=0;
                        O2_Grid_Occupancy(Cord_x1,Cord_y1)=0;
                        rxnflag1=1;
                    elseif order==2 && CO_Grid_Occupancy(Cord_x2,Cord_y2+1)==1
                        CO_Grid_Occupancy(Cord_x2,Cord_y2+1)=0;
                        O2_Grid_Occupancy(Cord_x2,Cord_y2)=0;
                        rxnflag2=1;
                    elseif order==3 && CO_Grid_Occupancy(Cord_x1-1,Cord_y1)==1
                        CO_Grid_Occupancy(Cord_x1-1,Cord_y1)=0;
                        O2_Grid_Occupancy(Cord_x1,Cord_y1)=0;
                        rxnflag1=1;
                    elseif order==4 && CO_Grid_Occupancy(Cord_x2-1,Cord_y2)==1
                        CO_Grid_Occupancy(Cord_x2-1,Cord_y2)=0;
                        O2_Grid_Occupancy(Cord_x2,Cord_y2)=0;
                        rxnflag2=1;
                    elseif order==5 && CO_Grid_Occupancy(Cord_x1+1,Cord_y1)==1
                        CO_Grid_Occupancy(Cord_x1+1,Cord_y1)=0;
                        O2_Grid_Occupancy(Cord_x1,Cord_y1)=0;
                        rxnflag1=1;
                    elseif order==6 && CO_Grid_Occupancy(Cord_x2+1,Cord_y2)==1
                        CO_Grid_Occupancy(Cord_x2+1,Cord_y2)=0;
                        O2_Grid_Occupancy(Cord_x2,Cord_y2)=0;
                        rxnflag2=1;
                        
                    end
                end
                CO2_COUNT=CO2_COUNT+rxnflag2+rxnflag1;
                
            end
        end
    end
    totalO2(i_trial)=sum(O2_Grid_Occupancy(:))/(N*N);
    totalCO(i_trial)=sum(CO_Grid_Occupancy(:))/(N*N);
    totalCO2(i_trial)=CO2_COUNT/(N*N);
end

Coverage_fraction_O=totalO2(end);
Coverage_fraction_CO=totalCO(end);
CO2_Production_rate=totalCO2(end);
