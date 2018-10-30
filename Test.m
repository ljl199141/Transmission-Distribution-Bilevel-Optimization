    clc;
    clear all;
    close all;
%% DISCO1 System Parameters
    nt = 12;
    mpc1 = case4gs;  
    [dbusA,~] = size(mpc1.bus);
    mpc1.bus(1,3) = 0;
    loadsA=repmat(mpc1.bus(:,3),1,nt); 
    GSF1 = makePTDF(mpc1.baseMVA, mpc1.bus, mpc1.branch, dbusA);
    PAup = repmat(mpc1.branch(:,6),1,nt); 
    [lineNum,~] = size(GSF1);
    total = loadsA;
    dssize = 1;
    cd2A = 0.02;
    cd1A = 5;
    PA = sdpvar(lineNum,nt,'full');
    pdA1 = sdpvar(dbusA,nt,'full');
    pdA1up = 0.01*total*dssize;
    pdA1dn = 0.00001*total*dssize;
    pdA = 0*(total-pdA1up);
    CpdA1 = 3;
    PAdn = 0;
    dgA = sdpvar(1,nt,'full');
    dgAup = 0;
    dgAdn = 0;
    drupA = sdpvar(dbusA,nt,'full');
    drdnA = sdpvar(dbusA,nt,'full');
    drA1 = 1;
    drA2 = 1;
    drscale = 0.3;
    cimA=14*ones(nt,1)';
    drpA=cimA;
    PinjDS=sdpvar(dbusA,nt,'full'); % bus nodal matrix with forecast wind   
%% DISCO1 Constraint
    CDA = [ PAdn<=PA,PA<=PAup, pdA1dn<=pdA1, pdA1<=pdA1up, 0<=drupA,drupA<=drscale*pdA1, 0<=drdnA,drdnA<=drscale*pdA1];

    mu1A = sdpvar(dbusA,nt,'full'); 
    mu2A = sdpvar(lineNum,nt,'full'); 
    mu3A = sdpvar(1,nt,'full'); 

    ST = [];
    bigM = [];
    DC = [];

    m1A = 10000000;

    l3A = sdpvar(lineNum,nt,'full'); 
    b3A = binvar(lineNum,nt,'full'); 
    %    DC = [DC,l3A.*(PA-PAdn) == 0];
    bigM = [bigM, l3A<=m1A*b3A, PA-PAdn<=m1A*(1-b3A)];

    l4A = sdpvar(lineNum,nt,'full');  
    b4A = binvar(lineNum,nt,'full');
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

    DF = [l3A>=0,l4A>=0,l5A>=0,l6A>=0];
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
% 
    DF = [l3A>=0,l4A>=0,l5A>=0,l6A>=0,l7A>=0,l8A>=0,l9A>=0,l10A>=0]; %
    
    for i=1:dbusA
      if i==1
          CDA=[CDA,PinjDS(i,:)==-pdA1(i,:)+PA(1,:)+PA(2,:)];
      else
          CDA=[CDA,PinjDS(i,:)==-pdA1(i,:)];
      end
    end
    CDA = [CDA,GSF1*PinjDS==PA];
    CDA = [CDA,PA(1,:)+PA(2,:) == sum(pdA + pdA1)];
    
    for i = 1:dbusA
       ST = [ST,2*CpdA1*pdA1(i,:)-2*CpdA1*pdA1up(i,:)+mu1A(i,:)-mu3A+l6A(i,:)-l5A(i,:)-drscale*l8A(i,:)-drscale*l10A(i,:) == 0]; 
    end %pdA1

    ST = [ST,2*drA2*drupA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l8A-l7A == 0,2*drA2*drdnA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l10A-l9A == 0];%drup,drdn

    for i = 1:lineNum
       if i == 1
          ST = [ST,cimA-mu1A(i,:)-mu2A(i,:)+mu3A+l4A(i,:)-l3A(i,:) == 0];
       elseif i == 2
          ST = [ST,cimA-mu1A(1,:)-mu2A(i,:)+mu3A+l4A(i,:)-l3A(i,:) == 0];
       else
          ST = [ST,-mu2A(i,:)+l4A(i,:)-l3A(i,:) == 0];
       end
    end %PA

    ST = [ST,mu1A+GSF1'*mu2A == 0];%PinjDS  
    
    dual1 = - sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2-sum(sum(pdA1.*pdA1)*CpdA1)+...
    sum(sum(pdA.*mu1A)+sum(l3A*PAdn-l4A.*PAup)+sum(l5A.*pdA1dn-l6A.*pdA1up))+sum(sum(pdA1up.*pdA1up))*CpdA1
    OD1 = sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 +...
    sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 -sum(drpA.*sum(drupA+drdnA))+ sum((PA(1,:)+PA(2,:)).*cimA)

    optimize([CDA],OD1)
    value(OD1)
    optimize([CDA,ST,bigM,DF],OD1)
    value(OD1)
    value(dual1)

    