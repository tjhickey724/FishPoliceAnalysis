e1 = importdata('gamelogs20151130.txt',' ',0);
e2 = importdata('gamelogs20160229.txt',' ',0);
e2pro = importdata('gamelogs20160229_pro.txt',' ',0);
e3 = importdata('gamelogs20160404.txt',' ',0);
e3pro = importdata('gamelogs20160404_pro.txt',' ',0);





% renumber so that each subject in each experiment has an id from 1-numsubjects
% exp1           1:16
% exp2 novice    1:12
% exp2 pro       1:2
% exp3 novice    1:8
% exp3 pro       1:7
%
e1(e1(:,4)==2,4)=4;  % change id of user 2 to 4, as no data for users 1,3,4
e1(:,4) = e1(:,4)-3;   % e1 subjects go from 1-15
e1 = e1(e1(:,4)>=1,:) ; % throw out the experimenters demo...

e2(:,4) = e2(:,4)-2 ; % e2 subjects go from 1-12;
e2pro(:,4) = e2pro(:,4)-4; %e2pro subjects go from 1-2
e3(:,4) = e3(:,4)-0; %e3 subjects go from 1-8
e3pro(:,4) = e3pro(:,4) - 1; %e3pro go from 1-6;
e3pro = e3pro(e3pro(:,4)>=1,:);
% add the experiment number as the 20th column
e1(:,20)=1;
e2(:,20)=2;
e2pro(:,20)=3;
e3(:,20)=4;
e3pro(:,20)=5;

% now create a table where each user has a unique id
% first we copy all of the tables
e1a=e1(:,:);
e2a=e2(:,:);
e2pa = e2pro(:,:);
e3a=e3(:,:);
e3pa= e3pro(:,:);
e2a(:,4) = e2a(:,4)+15;
e2pa(:,4) = e2pa(:,4)+27;
e3a(:,4) = e3a(:,4)+29;
e3pa(:,4) = e3pa(:,4)+37;


e = [e1a; e2a; e2pa; e3a; e3pa];

% plot the distribution of users
subplot(6,1,1); histogram(e1(:,4),-0.5:20);
subplot(6,1,2); histogram(e2(:,4),-0.5:20);
subplot(6,1,3); histogram(e2pro(:,4),-0.5:20);
subplot(6,1,4); histogram(e3(:,4),-0.5:20);
subplot(6,1,5); histogram(e3pro(:,4),-0.5:20);
subplot(6,1,6);
histogram(e(:,4),-0.5:45);
 

