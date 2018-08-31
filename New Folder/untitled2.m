 %% DISCO2 System Parameters
   dbusB = 6;
   pdB = [0, 1.5, 2.5, 0, 1.5, 0.75];
   pdB = repmat(pdB',1,nt);
   cd2B = 0.002;
   cd1B = 25;
   cimB = 30;
   
   PB = sdpvar(dbusB,nt,'full');
   pdB1 = sdpvar(dbusB,nt,'full');
   pdB1up = pdB;
   pdB1dn = 0.5*pdB;
   CpdB1 = 1;
   PBup = 30;
   PBdn = 0;
   dgB = sdpvar(1,nt,'full');
   dgBup = 10;
   dgBdn = 0;
   drupB = sdpvar(dbusB,nt,'full');
   drdnB = sdpvar(dbusB,nt,'full');
   drB1 = 1;
   drB2 = 1;
   drpB = 0;
%% DISCO1 Constraints
   CDB = [dgBdn<=dgB<=dgBup, PBdn<=PB<=PBup, pdB1dn<=pdB1<=pdB1up, 0<=drupB<=0.2*pdB1, 0<=drdnB<=0.2*pdB1];
   for i = 1:dbusB-1
       if i ~= dbusB-1
           CDB = [CDB, PB(i+1,:) == PB(i,:) - pdB(i+1,:) - pdB1(i+1,:)];
       else
           CDB = [CDB, PB(i+1,:) == PB(i,:) - pdB(i+1,:) - pdB1(i+1,:) + dgB];
       end
   end
   
   ODB = sum(sum(drupB'+drdnB'))*drB1 + sum(sum(drupB'.*drupB'+drdnB'.*drdnB'))*drB2 +sum(dgB')*cd1B +...
       sum(dgB'.*dgB')*cd2B + sum(sum((pdB1up-pdB1).*(pdB1up-pdB1)))*CpdB1 -...
       sum(sum((pdB1up).*(pdB1up)))*CpdB1 - sum(drpB.*sum(drupB+drdnB)) + sum(PB(1,:).*cimB);   
   
   mu1B = sdpvar(5,nt,'full'); 
   bigM = [];
   DC = [];
   ST = [];
   l1B = sdpvar(1,nt,'full'); 
   b1B = binvar(1,nt,'full'); 
   m1B = 10000;
   DC = [DC,l1B.*(dgB-dgBdn) == 0];
   bigM = [bigM, l1B<=m1B*b1B, dgB-dgBdn<=m1B*(1-b1B)];
   
   l2B = sdpvar(1,nt,'full');  
   b2B = binvar(1,nt,'full');
   DC = [DC,l2B.*(dgB-dgBup) == 0];
   bigM = [bigM, l2B<=m1B*b2B, -dgB+dgBup<=m1B*(1-b2B)];
   
   l3B = sdpvar(dbusB,nt,'full'); 
   b3B = binvar(dbusB,nt,'full'); 
   DC = [DC,l3B.*(PB-PBdn) == 0];
   bigM = [bigM, l3B<=m1B*b3B, PB-PBdn<=m1B*(1-b3B)];
   
   l4B = sdpvar(dbusB,nt,'full');  
   b4B = binvar(dbusB,nt,'full');
   DC = [DC,l4B.*(PB-PBup) == 0];
   bigM = [bigM, l4B<=m1B*b4B, -PB+PBup<=m1B*(1-b4B)];
   
   l5B = sdpvar(dbusB,nt,'full'); 
   b5B = binvar(dbusB,nt,'full'); 
   DC = [DC,l5B.*(pdB1-pdB1dn) == 0];
   bigM = [bigM, l5B<=m1B*b5B, pdB1-pdB1dn<=m1B*(1-b5B)];
   
   l6B = sdpvar(dbusB,nt,'full');  
   b6B = binvar(dbusB,nt,'full');
   DC = [DC,l6B.*(pdB1-pdB1up) == 0];
   bigM = [bigM, l6B<=m1B*b6B, -pdB1+pdB1up<=m1B*(1-b6B)];
   
   l7B = sdpvar(dbusB,nt,'full'); 
   b7B = binvar(dbusB,nt,'full'); 
   DC = [DC,l7B.*(drupB) == 0];
   bigM = [bigM, l7B<=m1B*b7B, drupB<=m1B*(1-b7B)];
   
   l8B = sdpvar(dbusB,nt,'full');  
   b8B = binvar(dbusB,nt,'full');
   DC = [DC,l8B.*(drupB-0.2*pdB1) == 0];
   bigM = [bigM, l8B<=m1B*b8B, -drupB+0.2*pdB1<=m1B*(1-b8B)];
   
   l9B = sdpvar(dbusB,nt,'full'); 
   b9B = binvar(dbusB,nt,'full'); 
   DC = [DC,l9B.*(drdnB) == 0];
   bigM = [bigM, l9B<=m1B*b9B, drdnB<=m1B*(1-b9B)];
   
   l10B = sdpvar(dbusB,nt,'full');  
   b10B = binvar(dbusB,nt,'full');
   DC = [DC,l10B.*(drdnB-0.2*pdB1) == 0];
   bigM = [bigM, l10B<=m1B*b10B, -drdnB+0.2*pdB1<=m1B*(1-b10B)];
 
   DF = [l1B>=0,l2B>=0,l3B>=0,l4B>=0,l5B>=0,l6B>=0,l7B>=0,l8B>=0,l9B>=0,l10B>=0];
   ST = [ST,2*cd2B*dgB+cd1B*ones(1,nt)+l2B-l1B-mu1B(5,:) == 0];%dgB
   ST = [ST,cimB-mu1B(1,:)+l4B(1,:)-l3B(1,:)==0, mu1B(1,:)-mu1B(2,:)+l4B(2,:)-l3B(2,:) == 0, mu1B(2,:)-mu1B(3,:)+l4B(3,:)-l3B(3,:) == 0, mu1B(3,:)-mu1B(4,:)+l4B(4,:)-l3B(4,:) == 0, mu1B(4,:)-mu1B(5,:)+l4B(5,:)-l3B(5,:) == 0, mu1B(5,:)+l4B(6,:)-l3B(6,:) == 0];%PB
   ST = [ST,2*CpdB1*pdB1(1,:)-2*CpdB1*pdB1up(1,:)+l6B(1,:)-l5B(1,:)-l8B(1,:)-l10B(1,:) == 0, 2*CpdB1*pdB1(2,:)-2*CpdB1*pdB1up(2,:)+mu1B(1,:)+l6B(2,:)-l5B(2,:)-l8B(2,:)-l10B(2,:) == 0,...
       2*CpdB1*pdB1(3,:)-2*CpdB1*pdB1up(3,:)+mu1B(2,:)+l6B(3,:)-l5B(3,:)-l8B(3,:)-l10B(3,:) == 0, 2*CpdB1*pdB1(4,:)-2*CpdB1*pdB1up(4,:)+mu1B(3,:)+l6B(4,:)-l5B(4,:)-l8B(4,:)-l10B(4,:) == 0,...
       2*CpdB1*pdB1(5,:)-2*CpdB1*pdB1up(5,:)+mu1B(4,:)+l6B(5,:)-l5B(5,:)-l8B(5,:)-l10B(5,:) == 0, 2*CpdB1*pdB1(6,:)-2*CpdB1*pdB1up(6,:)+mu1B(5,:)+l6B(6,:)-l5B(6,:)-l8B(6,:)-l10B(6,:) == 0];%pdB1
   ST = [ST,2*drB2*drupB+(drB1)*ones(dbusB,nt)-repmat(drpB,dbusB,1)+l8B-l7B == 0,2*drB2*drdnB+(drB1)*ones(dbusB,nt)-repmat(drpB,dbusB,1)+l10B-l9B == 0];%drup,drdn
   
   ODB
   dualB = -sum(dgB.*dgB)*cd2B+sum(-sum(pdB1.*pdB1)*CpdB1+sum(pdB(2:6,:).*mu1B)+sum(l3B*PBdn-l4B*PBup)+sum(l5B.*pdB1dn-l6B.*pdB1up)-l2B*dgBup)
   O2B = dual-(sum(sum(drupB'+drdnB'))*drB1 + sum(sum(drupB'.*drupB'+drdnB'.*drdnB'))*drB2 +sum(dgB')*cd1B...
       + sum(dgB'.*dgB')*cd2B + sum(sum((pdB1up-pdB1).*(pdB1up-pdB1)))*CpdB1 - sum(sum((pdB1up).*(pdB1up)))*CpdB1)