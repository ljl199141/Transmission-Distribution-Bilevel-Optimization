    clc;
    clear all;
    close all;
%% Transmission System Parameters
    nt=24; % time horizon
    mpc = case30;  
    [ng ~] = size(mpc.gen); 
    [bus,~] = size(mpc.bus);
    cg1 = mpc.gencost(:,6); % 1st order generator cost
    cg2 = mpc.gencost(:,5); % 2nd order generator cost
    crg=[6 6.75 7 5.25 5 5]';
    sampleL=[120.6 115.8 114.8 112.6 114.0 113.4 117.1 126.3 130.7 132.5 135.6 134.8 136.5 137.7 137.1 138.0 136.3 133.3 131.7 129.3 128.2 127.4 125.6 124.2];
    sf=sampleL(:,1:nt)./mean(sampleL(:,1:nt));
    sf=repmat(sf,bus,1);
    loads=repmat(mpc.bus(:,3),1,nt); 
    loads=loads.*sf;
    loads(2,:) = 0;
    loads(8,:) = 0;
    wbus1 = 5; % bus for wind1
    wbus2 = 4; % bus for wind1
    mbus1 = 9; % note
    mbus2 = 15; % note
    mbus3 = 8; % note
    wc = 1000; % wind spill cost
    lc = 1000; % load shedding cost
    GSF = makePTDF(mpc.baseMVA, mpc.bus, mpc.branch, bus);
    lineC = repmat(mpc.branch(:,6),1,nt);   
    Conoff=2*ones(1,ng);
    
%% DISCO1 System Parameters
   dbusA = 6;
   pdA = [0, 1.5, 2.5, 0, 1.5, 0.75];
   pdA = repmat(pdA',1,nt);
   cd2A = 0.002;
   cd1A = 25;
   cim = 30;
   
   PA = sdpvar(dbusA,nt,'full');
   PAup = 30;
   PAdn = 0;
   dgA = sdpvar(1,nt,'full');
   dgAup = 10;
   dgAdn = 0;
   CDA = [dgAdn<=dgA<=dgAup, PAdn<=PA<=PAup];
   for i = 1:dbusA-1
       if i ~= dbusA-1
           CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:)];
       else
           CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:) + dgA];
       end
   end
   
   OD1 = sum(dgA')*cd1A + sum(dgA'.*dgA')*cd2A + sum(PA(1,:))*cim;   
   
   
   bigM = [];
   DC = [];
   ST = [];
   l1A = sdpvar(1,nt,'full'); 
   b1A = binvar(1,nt,'full'); 
   m1A = 10000;
   DC = [DC,l1A.*(dgA-dgAdn) == 0];
   bigM = [bigM, l1A<=m1A*b1A, dgA-dgAdn<=m1A*(1-b1A)];
   
   l2A = sdpvar(1,nt,'full');  
   b2A = binvar(1,nt,'full');
   DC = [DC,l2A.*(dgA-dgAup) == 0];
   bigM = [bigM, l2A<=m1A*b2A, -dgA+dgAup<=m1A*(1-b2A)];
   
   l3A = sdpvar(dbusA,nt,'full'); 
   b3A = binvar(dbusA,nt,'full'); 
   DC = [DC,l3A.*(PA-PAdn) == 0];
   bigM = [bigM, l3A<=m1A*b3A, PA-PAdn<=m1A*(1-b3A)];
   
   l4A = sdpvar(dbusA,nt,'full');  
   b4A = binvar(dbusA,nt,'full');
   DC = [DC,l4A.*(PA-PAup) == 0];
   bigM = [bigM, l4A<=m1A*b4A, -PA+PAup<=m1A*(1-b4A)];
   
   mu1A = sdpvar(5,nt,'full'); 
 
   DF = [l1A>=0,l2A>=0,l3A>=0,l4A>=0];
   ST = [ST,2*cd2A*dgA+cd1A*ones(1,nt)+l2A-l1A-mu1A(5,:) == 0];
   ST = [ST, cim*ones(1,nt)-mu1A(1,:)+l4A(1,:)-l3A(1,:)==0, mu1A(1,:)-mu1A(2,:)+l4A(2,:)-l3A(2,:) == 0,...
       mu1A(2,:)-mu1A(3,:)+l4A(3,:)-l3A(3,:) == 0, mu1A(3,:)-mu1A(4,:)+l4A(4,:)-l3A(4,:) == 0, ...
       mu1A(4,:)-mu1A(5,:)+l4A(5,:)-l3A(5,:) == 0, mu1A(5,:)+l4A(6,:)-l3A(6,:) == 0];
%% Wind Forecast and Scenarios
    %wind forecast
    wl = 1; % wind penetration scaling factor
    load('ObsDays.mat');
    k = 54;
    wf1=wl*ObsDays(k,:,11);
    corr1=corrcoef(ObsDays(:,:,11));
    std_dev = [0.12,0.15,0.18,0.5,0.6,0.67,0.72,0.76,0.79,0.82,0.83,0.8315,0.833,0.835,0.836,0.838,0.839,0.841,0.842,0.844,0.845,0.847,0.848,0.85];
    for j=1:nt
        for k=1:nt
        covarr1(j,k)=corr1(j,k)*std_dev(k)*std_dev(j);
        end
    end
    beta=0.01;
    eulernum=exp(1);
    epsilon=0.01;
    Ndelta=1*nt;%num of uncertainties,one wind each period
    Nneed=ceil((1/epsilon)*(eulernum/(eulernum-1))*(log(1/beta)+4*(Ndelta+1)-1));%Number of scenarios
    
    %create wind error set
    avg=zeros(nt,1)';
    for i=1:Nneed
        err1(i,:)=mvnrnd(avg,covarr1);
    end
    for j=1:Nneed
        wn1(j,:)=wf1+wf1.*err1(j,:);
        pos=wn1(j,:)>=0;
        wst(j,:)=wn1(j,:).*pos;
    end
    winderror=wst-repmat(wf1,Nneed,1);   
    %create wind error set
    pr1=95;%upper percentile(throughout)
    pr2=5;%lower percentile(throughout)
    windup1=prctile(winderror,pr1);
    winddown1=prctile(winderror,pr2);

    % wind2
    wl=0; % wind penetration scaling factor
    wf2=wl*ObsDays(54,:,11);
    for j=1:Nneed
        wn1(j,:)=wf2+wf2.*err1(j,:);
        pos=wn1(j,:)>=0;
        wst(j,:)=wn1(j,:).*pos;
    end
    winderror=wst-repmat(wf2,Nneed,1);   
    
    %create wind error set
    windup2=prctile(winderror,pr1);
    winddown2=prctile(winderror,pr2);
    
    % wind1
    scale = 1;
    wf1 = scale*(wf1);
    windup1 = scale*(windup1);
    winddown1 = scale*(winddown1);

    % wind2
    scale2 = 0;
    wf2 = scale2*(wf2);
    windup2 = scale2*(windup2);
    winddown2 = scale2*(winddown2);
    
    % total wind
    wf = wf1+wf2;
    windup = windup1+windup2;
    winddown = winddown1+winddown2;
%% Transmission Optimization Variables and Their Bounds
    % Generators
    Gmax = [];
    for i = 1:ng
        Gmax = [Gmax;mpc.gen(i,9)*ones(1,nt)];
    end
    Gmin=[zeros(ng,nt)];
    pg=sdpvar(ng,nt,'full');
%     onoff=ones(ng,nt);
    onoff=binvar(ng,nt,'full');

    Rgmax=0.1*Gmax;
    rgup=sdpvar(ng,nt,'full');
    rgdn=sdpvar(ng,nt,'full');
    Rdn=0.3*Gmax;
    Rup=0.3*Gmax;
    
    % import & export prices
    Pim = sdpvar(1,nt,'full');
    Pimmax = 10*ones(1,nt);
    Pimmin = 0*ones(1,nt);
    
    % DR price
    Pdr=sdpvar(1,nt,'full');
    Pdrmax = 10;
    
    % LB & UB 
    CO = [Gmin<=pg<=Gmax,0<=Pdr<=Pdrmax,0<=rgup<=Rgmax,0<=rgdn<=Rgmax,Pimmin<=Pim<=Pimmax];      
    %% Transmission Constraints
    % power flow with forecast wind
    genbus = mpc.gen(:,1);
    genvec = zeros(1,bus);
    for i = 1:length(genbus)
        genvec(genbus(i))=1;
    end
    
    Pinj=sdpvar(bus,nt,'full'); % bus nodal matrix with forecast wind        
    count = 1;
    for i=1:bus
      if genvec(i) ~= 1
          if genvec(i) == wbus1
              CO = [CO,Pinj(i,:)==-loads(i,:)+wf1];
          elseif genvec(i) == wbus2
              CO = [CO,Pinj(i,:)==-loads(i,:)+wf2];
%           elseif genvec(i) == mbus1
%               CO = [CO,Pinj(i,:)==-loads(i,:)];
%           elseif genvec(i) == mbus2
%               CO = [CO,Pinj(i,:)==-loads(i,:)];
%           elseif genvec(i) == mbus3
%               CO = [CO,Pinj(i,:)==-loads(i,:)];
          else
              CO=[CO,Pinj(i,:)==-loads(i,:)];
          end
      else
          CO=[CO,Pinj(i,:)==-loads(i,:)+pg(count,:)];
          count = count+1;
      end
    end
    CO = [CO,-lineC<=GSF*Pinj<=lineC];
    
    Pinj1=sdpvar(bus,nt,'full'); % bus nodal matrix with forecast wind        
    count = 1;
    for i=1:bus
      if genvec(i) ~= 1
          if genvec(i) == wbus1
              CO = [CO,Pinj1(i,:)==-loads(i,:)+wf1+winddown1];
          elseif genvec(i) == wbus2
              CO = [CO,Pinj1(i,:)==-loads(i,:)+wf2+winddown2];
%           elseif genvec(i) == mbus1
%               CO = [CO,Pinj1(i,:)==-loads(i,:)+net+DRdn];
%           elseif genvec(i) == mbus2
%               CO = [CO,Pinj1(i,:)==-loads(i,:)+netC+DRdnC];
%           elseif genvec(i) == mbus3
%               CO = [CO,Pinj1(i,:)==-loads(i,:)+netB+DRdnB];
          else
              CO=[CO,Pinj1(i,:)==-loads(i,:)];
          end
      else
          CO=[CO,Pinj1(i,:)==-loads(i,:)+pg(count,:)+rgup(count,:)];
          count = count+1;
      end
    end
    CO = [CO,-lineC<=GSF*Pinj1<=lineC];

    Pinj2=sdpvar(bus,nt,'full'); % bus nodal matrix with forecast wind        
    count = 1;
    for i=1:bus
      if genvec(i) ~= 1
          if genvec(i) == wbus1
              CO = [CO,Pinj2(i,:)==-loads(i,:)+wf1+windup1];
          elseif genvec(i) == wbus2
              CO = [CO,Pinj2(i,:)==-loads(i,:)+wf2+windup2];
%           elseif genvec(i) == mbus1
%               CO = [CO,Pinj2(i,:)==-loads(i,:)+net-DRup];
%           elseif genvec(i) == mbus2
%               CO = [CO,Pinj2(i,:)==-loads(i,:)+netC-DRupC];
%           elseif genvec(i) == mbus3
%               CO = [CO,Pinj2(i,:)==-loads(i,:)+netB-DRupB];
          else
              CO=[CO,Pinj2(i,:)==-loads(i,:)];
          end
      else
          CO=[CO,Pinj2(i,:)==-loads(i,:)+pg(count,:)-rgdn(count,:)];
          count = count+1;
      end
    end
    CO = [CO,-lineC<=GSF*Pinj2<=lineC];
    
    % Generator Reserve Constraints
    CO = [CO,Gmin.*onoff<=pg<=Gmax.*onoff]; % Generator Bounds with Rgs
    CO = [CO,pg+rgup<=Gmax.*onoff,Gmin.*onoff<=pg-rgdn]; 
    CO = [CO,(pg(:,2:nt)+rgup(:,2:nt))-(pg(:,1:nt-1)-rgdn(:,1:nt-1))<=Rup(:,2:nt)]; % Up ramping Constraints with Rgs
    CO = [CO,-Rdn(:,2:nt)<=(pg(:,2:nt)-rgdn(:,2:nt))-(pg(:,1:nt-1)+rgup(:,1:nt-1))]; % Dn Ramping Constraints with Rgs
    CO = [CO,windup-sum(rgdn)==0,-winddown-sum(rgup)==0]; % note
  
    % Generator Constraints
    CO=[CO,-Rdn(:,2:nt)<=pg(:,2:nt)-pg(:,1:nt-1)<=Rup(:,2:nt)]; % ramping CO

    CO=[CO,sum(pg)-sum(loads)+wf==0];   %note
    %% Transmission Objective
    OO=sum(onoff')*Conoff'+sum(pg')*cg1+sum(rgup' + rgdn')*crg+sum(windup-sum(rgdn))*wc+sum(-winddown-sum(rgup))*lc;%note
    O1 = sum(pg'.*pg')*cg2;
%% Solve Problem
%    optimize([CDA],OD1)
%    optimize([CDA,ST,bigM,DF],0)
%    OD1
%    dual = sum(-sum(dgA.*dgA)*cd2A+sum(pdA(2:6,:).*mu1A)+sum(l3A*PAdn-l4A*PAup)-sum(l2A*dgAup));
%    OO1 = dual-(sum(dgA')*cd1A + sum(dgA'.*dgA')*cd2A)
   
   optimize([CO],OO+O1) % note
    % TransAvg = mean(costs)
    % pene1 = max(value(wf))/149
    % pene2 = (max(value(windup))+max(value(wf)))/149
    % MGAvg = mean(MGcost)
%% Plots
    % MG Gen VS Load:
%     t = 1:nt;
%     figure
%     subplot(6,1,1)
%     plot(t,value(wf),t,value(sum(pg)),t,value(ex))
%     legend('Wind Forecast','Total Generation','MG Export') 
%     title('Generation in the Transmission System')
% 
%     subplot(6,1,2)
%     plot(t,value(im),t,sum(loads))
%     legend('MG Import','Total Load')
%     title('Load in the Transmission System')
% 
%     subplot(6,1,3)
%     plot(t,value(mpg),t,value(im),t,max(0,-value(pb)))
%     legend('MG Generation','MG Import','Battery Discharging') 
%     title('Generation in the MG')
% 
%     subplot(6,1,4)
%     plot(t,value(ex),t,NCL,t,value(pd),t,max(0,value(pb)))
%     legend('MG Export','MG Non-controllable Load','MG Controllable Load','Battery Charging')
%     title('Load in the MG')
%     
%     subplot(6,1,5)
%     plot(t,value(windup),t,value(sum(rgdn)),t,value(DRup))
%     legend('Windup','Rgdn','DRup')
%     title('Windup and Rgdn in the Transmission System')
% 
%     subplot(6,1,6)
%     plot(t,value(-winddown),t,value(sum(rgup)),t,value(DRdn),t,value(Pdr))
%     legend('Winddown','Rgup','DRdn','DR Price') 
%     title('Winddn and Rgup in the MG')

% TransAvg = mean(costs)
% Pene1 = max(value(wf))/149
% Pene2 = (max(value(windup))+max(value(wf)))/149

% disp('lim1')
% lim1 = lineC-GSF*Pinj
% disp('lim2')
% lim2 = lineC-GSF*Pinj1
% disp('lim3')
% lim3 = lineC-GSF*Pinj2
% disp('lim11')
% lim11 = GSF*Pinj+lineC
% disp('lim22')
% lim22 = GSF*Pinj1+lineC
% disp('lim33')
% lim33 = GSF*Pinj2+lineC
% disp('Gmax-(pg+rgup)')
% Gmax-(pg+rgup)
% disp('pg-rgdn-Gmin')
% pg-rgdn-Gmin
% disp('Ramping1')
% Rup(:,2:nt)-(pg(:,2:nt)+rgup(:,2:nt))+(pg(:,1:nt-1)-rgdn(:,1:nt-1))
% disp('Ramping2')
% (pg(:,2:nt)-rgdn(:,2:nt))-(pg(:,1:nt-1)+rgup(:,1:nt-1))+Rdn(:,2:nt)
% disp('Reserve uplim')
% Rgmax-rgup
% disp('Reserve dnlim')
% Rgmax-rgdn