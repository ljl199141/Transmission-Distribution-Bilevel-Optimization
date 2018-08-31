%    clc
%    clear all
%    close all
%    nt = 24;
%    dbusA = 6;
%    pdA = [0, 1.5, 2.5, 0, 1.5, 0.75];
%    pdA = repmat(pdA',1,nt);
%    cd2A = 0.002;
%    cd1A = 25;
%    cim = 30;
%    
%    PA = sdpvar(dbusA,nt,'full');
%    pdA1 = sdpvar(dbusA,nt,'full');
%    pdA1up = pdA;
%    pdA1dn = 0.5*pdA;
%    PAup = 30;
%    PAdn = 0;
%    dgA = sdpvar(1,nt,'full');
%    dgAup = 10;
%    dgAdn = 0;
%    CDA = [dgAdn<=dgA<=dgAup, PAdn<=PA<=PAup, pdA1dn<=pdA1<=pdA1up];
%    for i = 1:dbusA-1
%        if i ~= dbusA-1
%            CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:) - pdA1(i+1,:)];
%        else
%            CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:) - pdA1(i+1,:) + dgA];
%        end
%    end
%    
%    OD1 = sum(dgA')*cd1A + sum(dgA'.*dgA')*cd2A + sum(PA(1,:))*cim;   
%    
%    
%    bigM = [];
%    DC = [];
%    ST = [];
%    l1A = sdpvar(1,nt,'full'); 
%    b1A = binvar(1,nt,'full'); 
%    m1A = 10000;
%    DC = [DC,l1A.*(dgA-dgAdn) == 0];
%    bigM = [bigM, l1A<=m1A*b1A, dgA-dgAdn<=m1A*(1-b1A)];
%    
%    l2A = sdpvar(1,nt,'full');  
%    b2A = binvar(1,nt,'full');
%    DC = [DC,l2A.*(dgA-dgAup) == 0];
%    bigM = [bigM, l2A<=m1A*b2A, -dgA+dgAup<=m1A*(1-b2A)];
%    
%    l3A = sdpvar(dbusA,nt,'full'); 
%    b3A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l3A.*(PA-PAdn) == 0];
%    bigM = [bigM, l3A<=m1A*b3A, PA-PAdn<=m1A*(1-b3A)];
%    
%    l4A = sdpvar(dbusA,nt,'full');  
%    b4A = binvar(dbusA,nt,'full');
%    DC = [DC,l4A.*(PA-PAup) == 0];
%    bigM = [bigM, l4A<=m1A*b4A, -PA+PAup<=m1A*(1-b4A)];
%    
%    l5A = sdpvar(dbusA,nt,'full'); 
%    b5A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l5A.*(pdA1-pdA1dn) == 0];
%    bigM = [bigM, l5A<=m1A*b5A, pdA1-pdA1dn<=m1A*(1-b5A)];
%    
%    l6A = sdpvar(dbusA,nt,'full');  
%    b6A = binvar(dbusA,nt,'full');
%    DC = [DC,l6A.*(pdA1-pdA1up) == 0];
%    bigM = [bigM, l6A<=m1A*b6A, -pdA1+pdA1up<=m1A*(1-b6A)];
%    
%    mu1A = sdpvar(5,nt,'full'); 
%  
%    DF = [l1A>=0,l2A>=0,l3A>=0,l4A>=0,l5A>=0,l6A>=0];
%    ST = [ST,2*cd2A*dgA+cd1A*ones(1,nt)+l2A-l1A-mu1A(5,:) == 0];
%    ST = [ST,cim*ones(1,nt)-mu1A(1,:)+l4A(1,:)-l3A(1,:)==0, mu1A(1,:)-mu1A(2,:)+l4A(2,:)-l3A(2,:) == 0, mu1A(2,:)-mu1A(3,:)+l4A(3,:)-l3A(3,:) == 0, mu1A(3,:)-mu1A(4,:)+l4A(4,:)-l3A(4,:) == 0, mu1A(4,:)-mu1A(5,:)+l4A(5,:)-l3A(5,:) == 0, mu1A(5,:)+l4A(6,:)-l3A(6,:) == 0];
%    ST = [ST,l6A(1,:)-l5A(1,:) == 0, mu1A(1,:)+l6A(2,:)-l5A(2,:) == 0, mu1A(2,:)+l6A(3,:)-l5A(3,:) == 0, mu1A(3,:)+l6A(4,:)-l5A(4,:) == 0, mu1A(4,:)+l6A(5,:)-l5A(5,:) == 0, mu1A(5,:)+l6A(6,:)-l5A(6,:) == 0];
%    optimize([CDA],OD1)
% 
%    optimize([CDA,ST,bigM,DF],0)
%    OD1
%    dual = sum(-sum(dgA.*dgA)*cd2A+sum(pdA(2:6,:).*mu1A)+sum(l3A*PAdn-l4A*PAup)+sum(l5A.*pdA1dn-l6A.*pdA1up)-sum(l2A*dgAup))
%    OO1 = dual-(sum(dgA')*cd1A + sum(dgA'.*dgA')*cd2A)

%    clc
%    clear all
%    close all
%    nt = 24;
%    dbusA = 6;
%    pdA = [0, 1.5, 2.5, 0, 1.5, 0.75];
%    pdA = repmat(pdA',1,nt);
%    cd2A = 0.002;
%    cd1A = 25;
%    cim = 30;
%    
%    PA = sdpvar(dbusA,nt,'full');
%    pdA1 = sdpvar(dbusA,nt,'full');
%    pdA1up = pdA;
%    pdA1dn = 0.5*pdA;
%    CpdA1 = 1;
%    PAup = 30;
%    PAdn = 0;
%    dgA = sdpvar(1,nt,'full');
%    dgAup = 10;
%    dgAdn = 0;
%    CDA = [dgAdn<=dgA<=dgAup, PAdn<=PA<=PAup, pdA1dn<=pdA1<=pdA1up];
%    for i = 1:dbusA-1
%        if i ~= dbusA-1
%            CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:) - pdA1(i+1,:)];
%        else
%            CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:) - pdA1(i+1,:) + dgA];
%        end
%    end
%    
%    OD1 = sum(dgA')*cd1A + sum(dgA'.*dgA')*cd2A + sum(PA(1,:))*cim + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 - sum(sum((pdA1up).*(pdA1up)))*CpdA1;   
%    
%    
%    bigM = [];
%    DC = [];
%    ST = [];
%    l1A = sdpvar(1,nt,'full'); 
%    b1A = binvar(1,nt,'full'); 
%    m1A = 10000;
%    DC = [DC,l1A.*(dgA-dgAdn) == 0];
%    bigM = [bigM, l1A<=m1A*b1A, dgA-dgAdn<=m1A*(1-b1A)];
%    
%    l2A = sdpvar(1,nt,'full');  
%    b2A = binvar(1,nt,'full');
%    DC = [DC,l2A.*(dgA-dgAup) == 0];
%    bigM = [bigM, l2A<=m1A*b2A, -dgA+dgAup<=m1A*(1-b2A)];
%    
%    l3A = sdpvar(dbusA,nt,'full'); 
%    b3A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l3A.*(PA-PAdn) == 0];
%    bigM = [bigM, l3A<=m1A*b3A, PA-PAdn<=m1A*(1-b3A)];
%    
%    l4A = sdpvar(dbusA,nt,'full');  
%    b4A = binvar(dbusA,nt,'full');
%    DC = [DC,l4A.*(PA-PAup) == 0];
%    bigM = [bigM, l4A<=m1A*b4A, -PA+PAup<=m1A*(1-b4A)];
%    
%    l5A = sdpvar(dbusA,nt,'full'); 
%    b5A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l5A.*(pdA1-pdA1dn) == 0];
%    bigM = [bigM, l5A<=m1A*b5A, pdA1-pdA1dn<=m1A*(1-b5A)];
%    
%    l6A = sdpvar(dbusA,nt,'full');  
%    b6A = binvar(dbusA,nt,'full');
%    DC = [DC,l6A.*(pdA1-pdA1up) == 0];
%    bigM = [bigM, l6A<=m1A*b6A, -pdA1+pdA1up<=m1A*(1-b6A)];
%    
%    mu1A = sdpvar(5,nt,'full'); 
%  
%    DF = [l1A>=0,l2A>=0,l3A>=0,l4A>=0,l5A>=0,l6A>=0];
%    ST = [ST,2*cd2A*dgA+cd1A*ones(1,nt)+l2A-l1A-mu1A(5,:) == 0];
%    ST = [ST,cim*ones(1,nt)-mu1A(1,:)+l4A(1,:)-l3A(1,:)==0, mu1A(1,:)-mu1A(2,:)+l4A(2,:)-l3A(2,:) == 0, mu1A(2,:)-mu1A(3,:)+l4A(3,:)-l3A(3,:) == 0, mu1A(3,:)-mu1A(4,:)+l4A(4,:)-l3A(4,:) == 0, mu1A(4,:)-mu1A(5,:)+l4A(5,:)-l3A(5,:) == 0, mu1A(5,:)+l4A(6,:)-l3A(6,:) == 0];
%    ST = [ST,2*CpdA1*pdA1(1,:)-2*CpdA1*pdA1up(1,:)+l6A(1,:)-l5A(1,:) == 0, 2*CpdA1*pdA1(2,:)-2*CpdA1*pdA1up(2,:)+mu1A(1,:)+l6A(2,:)-l5A(2,:) == 0, 2*CpdA1*pdA1(3,:)-2*CpdA1*pdA1up(3,:)+mu1A(2,:)+l6A(3,:)-l5A(3,:) == 0, 2*CpdA1*pdA1(4,:)-2*CpdA1*pdA1up(4,:)+mu1A(3,:)+l6A(4,:)-l5A(4,:) == 0, 2*CpdA1*pdA1(5,:)-2*CpdA1*pdA1up(5,:)+mu1A(4,:)+l6A(5,:)-l5A(5,:) == 0, 2*CpdA1*pdA1(6,:)-2*CpdA1*pdA1up(6,:)+mu1A(5,:)+l6A(6,:)-l5A(6,:) == 0];
%    optimize([CDA],OD1)
% 
%    optimize([CDA,ST,bigM,DF],0)
%    OD1
%    dual = sum(-sum(dgA.*dgA)*cd2A-sum(pdA1.*pdA1)*CpdA1+sum(pdA(2:6,:).*mu1A)+sum(l3A*PAdn-l4A*PAup)+sum(l5A.*pdA1dn-l6A.*pdA1up)-sum(l2A*dgAup))
%    OO1 = dual-(sum(dgA')*cd1A + sum(dgA'.*dgA')*cd2A)

   clc
   clear all
   close all
   nt = 24;
   dbusA = 6;
   pdA = [0, 1.5, 2.5, 0, 1.5, 0.75];
   pdA = repmat(pdA',1,nt);
   cd2A = 0.002;
   cd1A = 25;
   cimA = 10*ones(1,24);
   
   PA = sdpvar(dbusA,nt,'full');
   pdA1 = sdpvar(dbusA,nt,'full');
   pdA1up = pdA;
   pdA1dn = 0.5*pdA;
   CpdA1 = 1;
   PAup = 30;
   PAdn = 0;
   dgA = sdpvar(1,nt,'full');
   dgAup = 10;
   dgAdn = 0;
   drupA = sdpvar(dbusA,nt,'full');
   drdnA = sdpvar(dbusA,nt,'full');
   drA1 = 1;
   drA2 = 1;
   drpA = [1.30000000000142,1.30000000000144,1.30000000000000,1.30000000000002,1,1,1.30000000000143,1.30000000000002,1.04999999999984,1,1.10000000000000,1.30000000000142,1.15000000000000,1,1.10000000000000,1.29999999999931,1.50000000000000,1.50000000000000,1.50000000000000,1.50000000000000,1.30000000000145,1.30000000000145,1.30000000000000,1.50000000000000];
   
   CDA = [dgAdn<=dgA<=dgAup, PAdn<=PA<=PAup, pdA1dn<=pdA1<=pdA1up, 0<=drupA<=0.2*pdA1, 0<=drdnA<=0.2*pdA1];
   for i = 1:dbusA-1
       if i ~= dbusA-1
           CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:) - pdA1(i+1,:)];
       else
           CDA = [CDA, PA(i+1,:) == PA(i,:) - pdA(i+1,:) - pdA1(i+1,:) + dgA];
       end
   end
   
   OD1 = sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 +sum(dgA')*cd1A +...
       sum(dgA'.*dgA')*cd2A + sum(PA(1,:).*cimA) + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 - sum(sum((pdA1up).*(pdA1up)))*CpdA1-sum(drpA.*sum(drupA+drdnA));   
   
   mu1A = sdpvar(5,nt,'full'); 
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
   
   l5A = sdpvar(dbusA,nt,'full'); 
   b5A = binvar(dbusA,nt,'full'); 
   DC = [DC,l5A.*(pdA1-pdA1dn) == 0];
   bigM = [bigM, l5A<=m1A*b5A, pdA1-pdA1dn<=m1A*(1-b5A)];
   
   l6A = sdpvar(dbusA,nt,'full');  
   b6A = binvar(dbusA,nt,'full');
   DC = [DC,l6A.*(pdA1-pdA1up) == 0];
   bigM = [bigM, l6A<=m1A*b6A, -pdA1+pdA1up<=m1A*(1-b6A)];
   
   l7A = sdpvar(dbusA,nt,'full'); 
   b7A = binvar(dbusA,nt,'full'); 
   DC = [DC,l7A.*(drupA) == 0];
   bigM = [bigM, l7A<=m1A*b7A, drupA<=m1A*(1-b7A)];
   
   l8A = sdpvar(dbusA,nt,'full');  
   b8A = binvar(dbusA,nt,'full');
   DC = [DC,l8A.*(drupA-0.2*pdA1) == 0];
   bigM = [bigM, l8A<=m1A*b8A, -drupA+0.2*pdA1<=m1A*(1-b8A)];
   
   l9A = sdpvar(dbusA,nt,'full'); 
   b9A = binvar(dbusA,nt,'full'); 
   DC = [DC,l9A.*(drdnA) == 0];
   bigM = [bigM, l9A<=m1A*b9A, drdnA<=m1A*(1-b9A)];
   
   l10A = sdpvar(dbusA,nt,'full');  
   b10A = binvar(dbusA,nt,'full');
   DC = [DC,l10A.*(drdnA-0.2*pdA1) == 0];
   bigM = [bigM, l10A<=m1A*b10A, -drdnA+0.2*pdA1<=m1A*(1-b10A)];
 
   DF = [l1A>=0,l2A>=0,l3A>=0,l4A>=0,l5A>=0,l6A>=0,l7A>=0,l8A>=0,l9A>=0,l10A>=0];
   ST = [ST,2*cd2A*dgA+cd1A*ones(1,nt)+l2A-l1A-mu1A(5,:) == 0];%dgA
   ST = [ST,cimA-mu1A(1,:)+l4A(1,:)-l3A(1,:)==0, mu1A(1,:)-mu1A(2,:)+l4A(2,:)-l3A(2,:) == 0, mu1A(2,:)-mu1A(3,:)+l4A(3,:)-l3A(3,:) == 0, mu1A(3,:)-mu1A(4,:)+l4A(4,:)-l3A(4,:) == 0, mu1A(4,:)-mu1A(5,:)+l4A(5,:)-l3A(5,:) == 0, mu1A(5,:)+l4A(6,:)-l3A(6,:) == 0];%PA
   ST = [ST,2*CpdA1*pdA1(1,:)-2*CpdA1*pdA1up(1,:)+l6A(1,:)-l5A(1,:)-l8A(1,:)-l10A(1,:) == 0, 2*CpdA1*pdA1(2,:)-2*CpdA1*pdA1up(2,:)+mu1A(1,:)+l6A(2,:)-l5A(2,:)-l8A(2,:)-l10A(2,:) == 0,...
       2*CpdA1*pdA1(3,:)-2*CpdA1*pdA1up(3,:)+mu1A(2,:)+l6A(3,:)-l5A(3,:)-l8A(3,:)-l10A(3,:) == 0, 2*CpdA1*pdA1(4,:)-2*CpdA1*pdA1up(4,:)+mu1A(3,:)+l6A(4,:)-l5A(4,:)-l8A(4,:)-l10A(4,:) == 0,...
       2*CpdA1*pdA1(5,:)-2*CpdA1*pdA1up(5,:)+mu1A(4,:)+l6A(5,:)-l5A(5,:)-l8A(5,:)-l10A(5,:) == 0, 2*CpdA1*pdA1(6,:)-2*CpdA1*pdA1up(6,:)+mu1A(5,:)+l6A(6,:)-l5A(6,:)-l8A(6,:)-l10A(6,:) == 0];%pdA1
   ST = [ST,2*drA2*drupA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l8A-l7A == 0,2*drA2*drdnA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l10A-l9A == 0];%drup,drdn
   optimize([CDA],OD1)

   optimize([CDA,ST,bigM,DF],0)
   OD1
   dual = -sum(dgA.*dgA)*cd2A+sum(-sum(pdA1.*pdA1)*CpdA1+sum(pdA(2:6,:).*mu1A)+sum(l3A*PAdn-l4A*PAup)+sum(l5A.*pdA1dn-l6A.*pdA1up)-l2A*dgAup)
   OO1 = dual-(sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 +sum(dgA')*cd1A...
       + sum(dgA'.*dgA')*cd2A + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 - sum(sum((pdA1up).*(pdA1up)))*CpdA1)+drpA*sum(sum(drupA'+drdnA'))
