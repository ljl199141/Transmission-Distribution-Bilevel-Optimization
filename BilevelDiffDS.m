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
    wbus1 = 6; % bus for wind1
    wbus2 = 4; % bus for wind1
    dbus1 = 5; % note
    mbus2 = 15; % note
    mbus3 = 8; % note
    wc = 1000; % wind spill cost
    lc = 1000; % load shedding cost
    GSF = makePTDF(mpc.baseMVA, mpc.bus, mpc.branch, bus);
    lineC = repmat(mpc.branch(:,6),1,nt);   
    Conoff=2*ones(1,ng);    
%% DISCO1 System Parameters
   dbusA = 6;
   total = [0, 4, 5, 3, 4, 2];
   total = repmat(total',1,nt);
   dssize = 1;
   cd2A = 0.02;
   cd1A = 5;  
   PA = sdpvar(dbusA-1,nt,'full');
   pdA1 = sdpvar(dbusA,nt,'full');
   pdA1up = total*dssize;
   pdA1dn = 0.00001*total*dssize;
   pdA = total-pdA1up;
   CpdA1 = 3;
   PAup = repmat(0.8*[14.0000   8.2500    4.5000    6.5000    3.7500]',1,24); %[14.0000   10.2500    5.5000    5.5000    1.7500    0.0000]
   PAdn = 0;
   dgA = sdpvar(1,nt,'full');
   dgAup = 0;
   dgAdn = 0;
   drupA = sdpvar(dbusA,nt,'full');
   drdnA = sdpvar(dbusA,nt,'full');
   drA1 = 1;
   drA2 = 1;
   drscale = 0.5;
%% Wind Forecast and Scenarios
    %wind forecast
    wl = 1; % wind penetration scaling factor
    load('ObsDays.mat');
    k = 54;
    wf1=wl*ObsDays(k,:,11)+2;
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
    pr1=100;%upper percentile(throughout)
    pr2=0;%lower percentile(throughout)
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
    scale = 0.5;
    wf1 = scale*(wf1);
    windup1 = scale*(windup1);
    winddown1 = scale*(winddown1);
    
    windup1 = [2.14850277468970,2.38176933526874,1.74677853896453,2.60511542849926,2.35355880939503,2.41625640922597,9.26333481098347,4.72603747192732,3.64443832323432,3.54370589728489,3.77247031861537,5.32507128646877,3.89024919339414,3.43366439452552,3.79975232769689,5.91542122062703,18.0071565379186,20.3629232486021,27.8562658759252,16.9300600310706,14.6166013530343,12.8505562939878,9.60886741793157,16.5816874966563];
    winddown1 = [-2.29653290221721,-2.29980619610469,-1.86455639333191,-1.35000000000000,-1,-1,-3.40000000000000,-1.45000000000000,-1.05000000000000,-1,-1.10000000000000,-1.55000000000000,-1.15000000000000,-1,-1.10000000000000,-1.95000000000000,-5.70000000000000,-6.30000000000000,-9.10000000000000,-5.65000000000000,-4.55000000000000,-3.90000000000000,-3.15000000000000,-5.15000000000000];

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
    onoff=ones(ng,nt);
%     onoff=binvar(ng,nt,'full');

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
    
    cimA = sdpvar(1,nt,'full'); %note
    drpA = sdpvar(1,nt,'full'); %note
    % LB & UB 
    CO = [Gmin<=pg<=Gmax,0<=Pdr<=Pdrmax,0<=rgup<=Rgmax,0<=rgdn<=Rgmax,Pimmin<=Pim<=Pimmax,0<=cimA<=30,0<=drpA<=10];      
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
          if i == wbus1
              CO = [CO,Pinj(i,:)==-loads(i,:)+wf1];
          elseif i == wbus2
              CO = [CO,Pinj(i,:)==-loads(i,:)+wf2];
          elseif i == dbus1
              CO = [CO,Pinj(i,:)==-loads(i,:)-PA(1,:)-PA(2,:)];%note
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
          if i == wbus1
              CO = [CO,Pinj1(i,:)==-loads(i,:)+wf1+winddown1];
          elseif i == wbus2
              CO = [CO,Pinj1(i,:)==-loads(i,:)+wf2+winddown2];
          elseif i == dbus1
              CO = [CO,Pinj1(i,:)==-loads(i,:)-PA(1,:)-PA(2,:)+sum(drdnA)];%note
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
          if i == wbus1
              CO = [CO,Pinj2(i,:)==-loads(i,:)+wf1+windup1];
          elseif i == wbus2
              CO = [CO,Pinj2(i,:)==-loads(i,:)+wf2+windup2];
          elseif i == dbus1
              CO = [CO,Pinj2(i,:)==-loads(i,:)-PA(1,:)-PA(2,:)-sum(drupA)];%note
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
    CO = [CO,windup-sum(rgdn)-sum(drupA)==0,-winddown-sum(rgup)-sum(drdnA)==0]; % note
  
    % Generator Constraints
    CO=[CO,-Rdn(:,2:nt)<=pg(:,2:nt)-pg(:,1:nt-1)<=Rup(:,2:nt)]; % ramping CO

    CO=[CO,sum(pg)-sum(loads)+wf-PA(1,:)-PA(2,:)==0];   %note
    %% DISCO1 Constraints
  CDA = [dgAdn<=dgA<=dgAup, PAdn<=PA<=PAup, pdA1dn<=pdA1<=pdA1up, 0<=drupA<=drscale*pdA1, 0<=drdnA<=drscale*pdA1];
%    for i = 1:dbusA-1
%        if i == 1
%            CDA = [CDA, PA(4,:) == PA(1,:) - pdA(2,:) - pdA1(2,:)]; %mu1A(1,:)
%        elseif i == 2
%            CDA = [CDA, PA(3,:) == PA(2,:) - pdA(5,:) - pdA1(5,:)];
%        elseif i == 3
%            CDA = [CDA, 0 == PA(3,:) - pdA(6,:) - pdA1(6,:)];
%        elseif i == 4
%            CDA = [CDA, PA(5,:) == PA(4,:) - pdA(3,:) - pdA1(3,:)];
%        elseif i == 5
%            CDA = [CDA, 0 == PA(5,:) - pdA(4,:) - pdA1(4,:)]; %mu1A(5,:)
%        end
%    end  

%    mu1A = sdpvar(5,nt,'full'); 
%    bigM = [];
%    DC = [];
%    ST = [];
%    l1A = sdpvar(1,nt,'full'); 
%    b1A = binvar(1,nt,'full'); 
%    m1A = 10000000;
% %    DC = [DC,l1A.*(dgA-dgAdn) == 0];
%    bigM = [bigM, l1A<=m1A*b1A, dgA-dgAdn<=m1A*(1-b1A)];
%    
%    l2A = sdpvar(1,nt,'full');  
%    b2A = binvar(1,nt,'full');
% %    DC = [DC,l2A.*(dgA-dgAup) == 0];
%    bigM = [bigM, l2A<=m1A*b2A, -dgA+dgAup<=m1A*(1-b2A)];
%    
%    l3A = sdpvar(dbusA,nt,'full'); 
%    b3A = binvar(dbusA,nt,'full'); 
% %    DC = [DC,l3A.*(PA-PAdn) == 0];
%    bigM = [bigM, l3A<=m1A*b3A, PA-PAdn<=m1A*(1-b3A)];
%    
%    l4A = sdpvar(dbusA,nt,'full');  
%    b4A = binvar(dbusA,nt,'full');
% %    DC = [DC,l4A.*(PA-PAup) == 0];
%    bigM = [bigM, l4A<=m1A*b4A, -PA+PAup<=m1A*(1-b4A)];
%    
%    l5A = sdpvar(dbusA,nt,'full'); 
%    b5A = binvar(dbusA,nt,'full'); 
% %    DC = [DC,l5A.*(pdA1-pdA1dn) == 0];
%    bigM = [bigM, l5A<=m1A*b5A, pdA1-pdA1dn<=m1A*(1-b5A)];
%    
%    l6A = sdpvar(dbusA,nt,'full');  
%    b6A = binvar(dbusA,nt,'full');
% %    DC = [DC,l6A.*(pdA1-pdA1up) == 0];
%    bigM = [bigM, l6A<=m1A*b6A, -pdA1+pdA1up<=m1A*(1-b6A)];
%    
%    l7A = sdpvar(dbusA,nt,'full'); 
%    b7A = binvar(dbusA,nt,'full'); 
% %    DC = [DC,l7A.*(drupA) == 0];
%    bigM = [bigM, l7A<=m1A*b7A, drupA<=m1A*(1-b7A)];
%    
%    l8A = sdpvar(dbusA,nt,'full');  
%    b8A = binvar(dbusA,nt,'full');
% %    DC = [DC,l8A.*(drupA-drscale*pdA1) == 0];
%    bigM = [bigM, l8A<=m1A*b8A, -drupA+drscale*pdA1<=m1A*(1-b8A)];
%    
%    l9A = sdpvar(dbusA,nt,'full'); 
%    b9A = binvar(dbusA,nt,'full'); 
% %    DC = [DC,l9A.*(drdnA) == 0];
%    bigM = [bigM, l9A<=m1A*b9A, drdnA<=m1A*(1-b9A)];
%    
%    l10A = sdpvar(dbusA,nt,'full');  
%    b10A = binvar(dbusA,nt,'full');
% %    DC = [DC,l10A.*(drdnA-drscale*pdA1) == 0];
%    bigM = [bigM, l10A<=m1A*b10A, -drdnA+drscale*pdA1<=m1A*(1-b10A)];
%        
%    DF = [l1A>=0,l2A>=0,l3A>=0,l4A>=0,l5A>=0,l6A>=0,l7A>=0,l8A>=0,l9A>=0,l10A>=0];
%    ST = [ST,cimA-mu1A(1,:)+l4A(1,:)-l3A(1,:)==0, cimA-mu1A(2,:)+l4A(2,:)-l3A(2,:) == 0, mu1A(2,:)-mu1A(3,:)+l4A(3,:)-l3A(3,:) == 0, ...
%        mu1A(1,:)-mu1A(4,:)+l4A(4,:)-l3A(4,:) == 0, mu1A(4,:)-mu1A(5,:)+l4A(5,:)-l3A(5,:) == 0, mu1A(5,:)+l4A(6,:)-l3A(6,:) == 0];%PA
%    ST = [ST,2*CpdA1*pdA1(1,:)-2*CpdA1*pdA1up(1,:)+l6A(1,:)-l5A(1,:)-drscale*l8A(1,:)-drscale*l10A(1,:) == 0,...
%        2*CpdA1*pdA1(2,:)-2*CpdA1*pdA1up(2,:)+mu1A(1,:)+l6A(2,:)-l5A(2,:)-drscale*l8A(2,:)-drscale*l10A(2,:) == 0,...
%        2*CpdA1*pdA1(3,:)-2*CpdA1*pdA1up(3,:)+mu1A(4,:)+l6A(3,:)-l5A(3,:)-drscale*l8A(3,:)-drscale*l10A(3,:) == 0,...
%        2*CpdA1*pdA1(4,:)-2*CpdA1*pdA1up(4,:)+mu1A(5,:)+l6A(4,:)-l5A(4,:)-drscale*l8A(4,:)-drscale*l10A(4,:) == 0,...
%        2*CpdA1*pdA1(5,:)-2*CpdA1*pdA1up(5,:)+mu1A(2,:)+l6A(5,:)-l5A(5,:)-drscale*l8A(5,:)-drscale*l10A(5,:) == 0,...
%        2*CpdA1*pdA1(6,:)-2*CpdA1*pdA1up(6,:)+mu1A(3,:)+l6A(6,:)-l5A(6,:)-drscale*l8A(6,:)-drscale*l10A(6,:) == 0];%pdA1      
%    ST = [ST,2*drA2*drupA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l8A-l7A == 0,2*drA2*drdnA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l10A-l9A == 0];%drup,drdn
%    ST = [ST,cimA == drpA];%dgA
   
  base1 = 100;
  busNum = 6;
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
  bus1 = [
	1	3	0	0	0	0	1	1	0	135	1	1.05	0.95;
	2	1	4	0	0	0	1	1	0	135	1	1.05	0.95;
	3	1	5	0	0	0	1	1	0	135	1	1.05	0.95;
	4	1	3	0	0	0	1	1	0	135	1	1.05	0.95;
	5	1	4	0	0	0	1	1	0	135	1	1.05	0.95;
	6	1	2	0	0	0	1	1	0	135	1	1.05	0.95;];

%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
  branch1 = [
	1	2	0.02	0.06	0.03	11.2	130	130	0	0	1	-360	360;
	1	5	0.05	0.19	0.02	6.6	130	130	0	0	1	-360	360;
	5	6	0.06	0.17	0.02	3.6	65	65	0	0	1	-360	360;
	2	3	0.01	0.04	0	5.20	130	130	0	0	1	-360	360;
	3	4	0.05	0.2	0.02	3	130	130	0	0	1	-360	360;];

  GSF1 = makePTDF(base1, bus1, branch1, busNum);

  PinjDS=sdpvar(busNum,nt,'full'); % bus nodal matrix with forecast wind        
    for i=1:busNum
      if i==1
          CDA=[CDA,PinjDS(i,:)==-pdA1(i,:)+PA(1,:)+PA(2,:)];
      else
          CDA=[CDA,PinjDS(i,:)==-pdA1(i,:)];
      end
    end
  CDA = [CDA,GSF1*PinjDS==PA];
  CDA = [CDA,PA(1,:)+PA(2,:) == sum(pdA + pdA1)];
   
   mu1A = sdpvar(6,nt,'full'); 
   mu2A = sdpvar(5,nt,'full'); 
   mu3A = sdpvar(1,nt,'full'); 
   
   ST = [];
   bigM = [];
   DC = [];

   l1A = sdpvar(1,nt,'full'); 
   b1A = binvar(1,nt,'full'); 
   m1A = 10000000;
%    DC = [DC,l1A.*(dgA-dgAdn) == 0];
   bigM = [bigM, l1A<=m1A*b1A, dgA-dgAdn<=m1A*(1-b1A)];
   
   l2A = sdpvar(1,nt,'full');  
   b2A = binvar(1,nt,'full');
%    DC = [DC,l2A.*(dgA-dgAup) == 0];
   bigM = [bigM, l2A<=m1A*b2A, -dgA+dgAup<=m1A*(1-b2A)];
   
   l3A = sdpvar(dbusA-1,nt,'full'); 
   b3A = binvar(dbusA-1,nt,'full'); 
%    DC = [DC,l3A.*(PA-PAdn) == 0];
   bigM = [bigM, l3A<=m1A*b3A, PA-PAdn<=m1A*(1-b3A)];
   
   l4A = sdpvar(dbusA-1,nt,'full');  
   b4A = binvar(dbusA-1,nt,'full');
%    DC = [DC,l4A.*(PA-PAup) == 0];
   bigM = [bigM, l4A<=m1A*b4A, -PA+PAup<=m1A*(1-b4A)];
   
   l5A = sdpvar(dbusA,nt,'full'); 
   b5A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l5A.*(pdA1-pdA1dn) == 0];
   bigM = [bigM, l5A<=m1A*b5A, pdA1-pdA1dn<=m1A*(1-b5A)];
   
   l6A = sdpvar(dbusA,nt,'full');  
   b6A = binvar(dbusA,nt,'full');
%    DC = [DC,l6A.*(pdA1-pdA1up) == 0];
   bigM = [bigM, l6A<=m1A*b6A, -pdA1+pdA1up<=m1A*(1-b6A)];
   
   l7A = sdpvar(dbusA,nt,'full'); 
   b7A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l7A.*(drupA) == 0];
   bigM = [bigM, l7A<=m1A*b7A, drupA<=m1A*(1-b7A)];
   
   l8A = sdpvar(dbusA,nt,'full');  
   b8A = binvar(dbusA,nt,'full');
%    DC = [DC,l8A.*(drupA-drscale*pdA1) == 0];
   bigM = [bigM, l8A<=m1A*b8A, -drupA+drscale*pdA1<=m1A*(1-b8A)];
   
   l9A = sdpvar(dbusA,nt,'full'); 
   b9A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l9A.*(drdnA) == 0];
   bigM = [bigM, l9A<=m1A*b9A, drdnA<=m1A*(1-b9A)];
   
   l10A = sdpvar(dbusA,nt,'full');  
   b10A = binvar(dbusA,nt,'full');
%    DC = [DC,l10A.*(drdnA-drscale*pdA1) == 0];
   bigM = [bigM, l10A<=m1A*b10A, -drdnA+drscale*pdA1<=m1A*(1-b10A)];
       
   DF = [l1A>=0,l2A>=0,l3A>=0,l4A>=0,l5A>=0,l6A>=0,l7A>=0,l8A>=0,l9A>=0,l10A>=0];
   for i = 1:dbusA
       if i == 1
          ST = [ST,cimA-mu1A(i,:)-mu2A(i,:)+mu3A+l4A(i,:)-l3A(i,:) == 0];
          ST = [ST,2*CpdA1*pdA1(i,:)-2*CpdA1*pdA1up(i,:)+mu1A(i,:)-mu3A+l6A(i,:)-l5A(i,:)-drscale*l8A(i,:)-drscale*l10A(i,:) == 0];
       elseif i == 2
          ST = [ST,cimA-mu1A(1,:)-mu2A(i,:)+mu3A+l4A(i,:)-l3A(i,:) == 0];
          ST = [ST,2*CpdA1*pdA1(i,:)-2*CpdA1*pdA1up(i,:)+mu1A(i,:)-mu3A+l6A(i,:)-l5A(i,:)-drscale*l8A(i,:)-drscale*l10A(i,:) == 0];
       elseif i == 6
          ST = [ST,2*CpdA1*pdA1(i,:)-2*CpdA1*pdA1up(i,:)+mu1A(i,:)-mu3A+l6A(i,:)-l5A(i,:)-drscale*l8A(i,:)-drscale*l10A(i,:) == 0];
       else
          ST = [ST,-mu2A(i,:)+l4A(i,:)-l3A(i,:) == 0];
          ST = [ST,2*CpdA1*pdA1(i,:)-2*CpdA1*pdA1up(i,:)+mu1A(i,:)-mu3A+l6A(i,:)-l5A(i,:)-drscale*l8A(i,:)-drscale*l10A(i,:) == 0]; 
       end
    end
   ST = [ST,mu1A+GSF1'*mu2A == 0];%PinjDS  
   ST = [ST,2*drA2*drupA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l8A-l7A == 0,2*drA2*drdnA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l10A-l9A == 0];%drup,drdn
   ST = [ST,cimA == drpA];%dgA
%% Transmission Objective
    OO = sum(onoff')*Conoff'+sum(pg')*cg1+sum(rgup' + rgdn')*crg+sum(windup-sum(rgdn)-sum(drupA))*wc+sum(-winddown-sum(rgup)-sum(drdnA))*lc;%note
    O1 = sum(pg'.*pg')*cg2;
%% Solve Problem
%    optimize([CDA],OD1)
%    optimize([CDA,ST,bigM,DF],0)
%    CpdA1 = 0;
%    dual = -sum(dgA.*dgA)*cd2A- sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2-sum(sum(pdA1.*pdA1)*CpdA1)+...
%        sum(sum(pdA(2:6,:).*mu1A)+sum(l3A*PAdn-l4A.*PAup)+sum(l5A.*pdA1dn-l6A.*pdA1up)-l2A*dgAup)+sum(sum(pdA1up.*pdA1up))*CpdA1
   dual = -sum(dgA.*dgA)*cd2A- sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2-sum(sum(pdA1.*pdA1)*CpdA1)+...
       sum(sum(pdA.*mu1A)+sum(l3A*PAdn-l4A.*PAup)+sum(l5A.*pdA1dn-l6A.*pdA1up)-l2A*dgAup)+sum(sum(pdA1up.*pdA1up))*CpdA1
   OD1 = sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 +sum(dgA')*cd1A +...
   sum(dgA'.*dgA')*cd2A  + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 -sum(drpA.*sum(drupA+drdnA))+ sum((PA(1,:)+PA(2,:)).*cimA) %- sum(sum(pdA1up.*pdA1up))*CpdA1
   
   O2 = dual-(sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 +sum(dgA')*cd1A...
       + sum(dgA'.*dgA')*cd2A + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1) %- sum(sum(pdA1up.*pdA1up))*CpdA1
%    O2 = O2-(sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 - sum(sum(pdA1up.*pdA1up))*CpdA1 + sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2)
   
   optimize([CDA,ST,bigM,DF,CO],OO+O1-O2) % note -10000*scale
%    value(cimA)
%    value(drpA)
%    value(dgA)
%    value(drupA)
%    value(drdnA)
   TransAvg = value(OO+O1-O2)+OD1
    % pene1 = max(value(wf))/149
%    pene2 = (max(value(windup))+max(value(wf)))/149
   OD1
   sum(PA(1,:)+PA(2,:))
   sum(sum(drupA+drdnA))
%    disp('total windup')
%    sum(windup)
%    disp('total winddown')
%    sum(winddown)
%    disp('energy price')
%    EnergyPrice = value(cimA)
%    disp('drup')
%    value(sum(drupA))
%    disp('drdn')
%    value(sum(drdnA))
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