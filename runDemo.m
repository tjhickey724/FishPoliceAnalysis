function [ out ] = runDemo(varargin )
%RunDemo This runs various FishPolice Data Analyses
%   Detailed explanation goes here later ...
% load the database with a command like this ..
% z = importdata('ExpData20160229Matlab.txt',' ',1);
% then call it with runDemo(z.data,'rt1');   
    fpdata = varargin{1};
    cmd = varargin{2};

    timeInSecs = fpdata(:,1);
    eventType=fpdata(:,2);
    objectId=fpdata(:,3);
    gameId = objectId;  % alternate name for column 3
    
    userId=fpdata(:,4);
    level=fpdata(:,5);
    age=fpdata(:,6);
    mode=fpdata(:,7);
    gameTime=fpdata(:,8);
    eventTime=fpdata(:,9);
    eventId=fpdata(:,10);
    visual=fpdata(:,11);
    audio=fpdata(:,12);
    side=fpdata(:,13);
    action=fpdata(:,14);
    actionTime=fpdata(:,15);
    reaction=fpdata(:,16);
    key=fpdata(:,17);
    correct=fpdata(:,18);
    congruent=fpdata(:,19);

    rowTime = actionTime; %gameTime+eventTime;
    rowTimeN = (rowTime(:)-rowTime(1))./60000; % normalized row time
    levelTime = level*2000+reaction;
    levelTimeN = levelTime./2000;
    oddball = visual==0 | audio==0;
    
    [a,b] = histcounts(userId,1:100);
    users = b(a>0);
    [a,b] = histcounts(gameId,1:2000);
    games = b(a>0);
   % display(users);
   % display(games);
    
   % next we calculate the level in visual(1) and audio(2) each user rose
   % to:
   vlevel=[]; alevel=[];
   for u=users 
       vlevel(u) = max(level(userId==u&mode==1));
       alevel(u) = max(level(userId==u&mode==2));
   end
   % next we find the users who reached at least level 7 on visual and
   % audio
   users7 = [];
   userId7 = userId==(userId+1);
   u7V = userId==(userId+1);
   u7A = userId==(userId+1);
   
   for u=users 
       if vlevel(u)>=7 && alevel(u)>=7
           userId7 = userId7 | userId==u;
           users7 = [u,users7];
       end
       if vlevel(u)>=7
           u7V= u7V | userId==u;
       end
       if alevel(u)>=7
           u7A = u7A | userId==u;
       end
   end
   %display(users7);
   u7 = (userId==0);
   ua7= (userId==0);
   uv7= (userId==0);
   for u=users7
       %urt(u) = mean(reaction(selection&userId==u));
       u7 = u7 | userId==u;
   end

   % next we calculate the average reaction time by user
   % and use that to calculate relative reaction times in our graphs
   userRT = userId*0;
   for u=users
       userRT = userRT + (userId==u)*sum(reaction(userId==u&eventType==2))/sum(userId==u&eventType==2);
   end
   relRT = reaction -userRT;

   
   ifiPrev=zeros(length(eventId),1)+2000;
lastBirth=1;
lastAction=1;
lastGame= gameId(1);
for i=[2:length(fpdata)]
    %display([i,gameId(i),lastGame,lastAction,eventType(i),eventId(i)])
    if (gameId(i) ~= lastGame)
        if eventType(i)==1
            lastGame= gameId(i);
            lastBirth=i;
        end
        
    elseif eventType(i)==1 
        %display([i,gameId(i),ifi,lastBirth,gameId(lastBirth)]);
        %display([actionTime(i)-actionTime(lastAction)]);
        ifi = eventTime(i)-eventTime(lastBirth)-reaction(lastAction);
        for j=[lastBirth+1:i-1]
            ifiPrev(j)=ifiPrev(lastBirth);
            ifiPrev(i) = ifi;
        end
        if eventType(i+1)==2 || eventType(i+1)==4
            ifiPrev(i+1) = ifi;
        end
        lastBirth = i;
    elseif eventType(i)==2 || eventType(i)==4
        lastAction=i;
    end
end
    
    switch cmd
        case 'test'
            histogram(relRT(eventType==2),-1000:25:1000);
        case 'eventID'
            constraints = varargin{3};
            eventRT=zeros(60,1);
            eventCor=zeros(60,1);
            for i=[1:60]
                eventCor(i)=mean(correct(constraints&eventId==i));
                eventRT(i)=mean(reaction(constraints&eventId==i)/1000.0);
            end
            plot(eventRT,'-+r'); hold on;
            plot(eventCor,'-xb');
        case 'ifi'
            constraints = varargin{3};
            color = varargin{4};
            scatter(ifiPrev(constraints&eventType==2), reaction(constraints&eventType==2)/1000.0,strcat('.',color));
            hold on;
            p = polyfit(ifiPrev(constraints&eventType==2),correct(constraints&eventType==2),10);
            bounds = [1000-75*max(level(constraints&eventType==2)):10:2000-150*min(constraints&eventType==2)];
            plot(bounds,polyval(p,bounds),strcat('-',color));
            axis([0,2000,0,1]); hold off
            
        case 'ifireact'
            constraints = varargin{3};
            color = varargin{4};
            scatter(ifiPrev(constraints&eventType==2), reaction(constraints&eventType==2)/1000.0,strcat('.',color));
            hold on;
            p = polyfit(ifiPrev(constraints&eventType==2),0.001*reaction(constraints&eventType==2),10);
            bounds = [1000-75*max(level(constraints&eventType==2)):10:2000-150*min(level(constraints&eventType==2))];
            plot(bounds,polyval(p,bounds),strcat('+',color));
            axis([0,2000,0,1]); 
            
        case 'ifireactAll'
            constraints = varargin{3};
            hold off;
            runDemo(fpdata,'ifireact',constraints&level<=5,'m');
            hold on;
            runDemo(fpdata,'ifireact',constraints&level==6,'r');
            runDemo(fpdata,'ifireact',constraints&level==7,'b');
            runDemo(fpdata,'ifireact',constraints&level==8,'k');
            runDemo(fpdata,'ifireact',constraints&level==9,'c');
            legend('<=5','<=5','6','6','7','7','8','8','9','9');
            
           
            case 'ifiCorrect'
            constraints = varargin{3};
            color = varargin{4};
            scatter(ifiPrev(constraints&eventType==2), reaction(constraints&eventType==2)/1000.0,strcat('.',color));
            hold on;
            p = polyfit(ifiPrev(constraints&eventType==2),correct(constraints&eventType==2),4);
            bounds = [1000-75*max(level(constraints&eventType==2)):10:2000-150*min(level(constraints&eventType==2))];
            plot(bounds,polyval(p,bounds),strcat('+',color));
            axis([0,2000,0,1]); 
            
        case 'ifiCorrectAll'
            constraints = varargin{3};
            %hold off;
            runDemo(fpdata,'ifiCorrect',constraints&level<=5,'c');
            hold on;
            runDemo(fpdata,'ifiCorrect',constraints&level==6,'r');
            runDemo(fpdata,'ifiCorrectt',constraints&level==7,'b');
            runDemo(fpdata,'ifiCorrect',constraints&level==8,'k');
            runDemo(fpdata,'ifiCorrect',constraints&level==9,'m');
            legend('<=5','<=5','6','6','7','7','8','8','9','9');
            
        case 'ifilev'
            constraints = varargin{3};
            lev = varargin{4};
            color = varargin{5};
            scatter(ifiPrev(constraints&eventType==2&level==lev), reaction(constraints&eventType==2&level==lev)/1000.0,strcat('.',color));
            hold on;
            p = polyfit(ifiPrev(constraints&eventType==2&level==lev),correct(constraints&eventType==2&level==lev),10);
            bounds = [1000-75*lev:100:2000-150*lev];
            plot(bounds,polyval(p,bounds),strcat('-',color));
            axis([0,2000,0,1]); hold off
            
        case 'learning'
            constraints = varargin{3};
            pcorrect = [];
            levs = [1:10]
            for lev=levs
                selection =   level==lev & constraints;
                tries(lev) = sum(selection & (eventType==2 | eventType==3))/sum(selection & eventType==1);
                misses = constraints& level==lev & eventType==3;
                missrt = sum(actionTime(misses)-gameTime(misses));
                pcorrect(lev)=sum(correct(selection&eventType==2))/sum(selection&eventType==2)
                rt(lev) = sum(reaction(selection&correct==1))/sum(selection&correct==1);
                rt0(lev) = sum(reaction(selection&eventType==2)/sum(selection&eventType==2));
                rt2(lev) = (sum(reaction(selection&eventType==2&correct==1))+missrt)/(sum(selection&eventType==2&correct==1)+sum(misses));
                miss(lev) = sum(u7&level==lev&eventType==3)/sum(u7&level==lev&eventType==1);
            end
            plot(rt./1000);hold on; 
            plot(rt0./1000); plot(rt2./1000);plot(tries);  
            plot(pcorrect);%plot(miss);
            legend('rt','rt0','rt2','tries','pcorrect');
            display(tries);
            grid on;
            hold off;
            

            
        case 'accVrtFishType'
            constraints = varargin{3};
                     runDemo(fpdata,'accVrt',constraints&visual==8&audio==8,'+','r');
            hold on; runDemo(fpdata,'accVrt',constraints&visual==5&audio==5,'x','b');
            hold on; runDemo(fpdata,'accVrt',constraints&visual==8&audio==5,'d','m');
            hold on; runDemo(fpdata,'accVrt',constraints&visual==5&audio==8,'s','k');
            hold on; runDemo(fpdata,'accVrt',constraints&visual==0,'o','b');
            hold on; runDemo(fpdata,'accVrt',constraints&audio==0,'*','c');
            legend('FF fit','FF','FF err','SS fit','SS','SS err','FS fit','FS','FS err','SF fit','SF','SF err','Hidden fit','Hidden','Hidden err','Silent fit','Silent','Silent err');
            hold off;
            
        case 'accVrtFishTypeSmooth'
            constraints = varargin{3};
                     runDemo(fpdata,'accVrtSmooth',constraints&visual==8&audio==8,'+','r');
            hold on; runDemo(fpdata,'accVrtSmooth',constraints&visual==5&audio==5,'x','b');
            hold on; runDemo(fpdata,'accVrtSmooth',constraints&visual==8&audio==5,'d','m');
            hold on; runDemo(fpdata,'accVrtSmooth',constraints&visual==5&audio==8,'s','k');
            hold on; runDemo(fpdata,'accVrtSmooth',constraints&visual==0,'o','g');
            hold on; runDemo(fpdata,'accVrtSmooth',constraints&audio==0,'*','c');
            legend('FF fit','FF','FF err','SS fit','SS','SS err','FS fit','FS','FS err','SF fit','SF','SF err','Hidden fit','Hidden','Hidden err','Silent fit','Silent','Silent err');
            hold off;
        
        case 'accVrt'
            constraints = varargin{3};
            theMarker = varargin{4};
            theColor = varargin{5};
            
            selection = constraints & eventType==2;
            x = reaction(selection)/1000;
            y = correct(selection);
            x1=0.4:0.01:1.4;
            p = polyfit(x,y,6);y1=polyval(p,x1);plot(x1,y1,strcat('-',theColor));
            axis([0,2,0.0,1]); hold on;
            bins100=350:100:1450;
            mTF = histcounts(reaction(selection),bins100);
            mT = histcounts(reaction(selection&correct==1),bins100);
            display(strcat('number of keystrokes: ',num2str(sum(mTF))));
            %plot(0.6:0.1:1.2,mT./mTF,strcat('-',theMarker));
            p = mT./mTF;
            errs = 2*sqrt(p.*(1-p)./mTF);
            errorbar(0.4:0.1:1.4,p,errs,strcat(theMarker,theColor));
            plot(0.4:0.1:1.4,mTF./sum(mTF),strcat('-',theMarker,theColor));
            %plot(0.4:0.1:1.4,errs,strcat('*',theColor));
            grid on;
            
        case 'relaccVrt'
            constraints = varargin{3};
            theMarker = varargin{4};
            theColor = varargin{5};
            
            selection = constraints & eventType==2;
            x = relRT(selection)/1000;
            y = correct(selection);
            x1=-0.4:0.01:0.4;
            p = polyfit(x,y,6);y1=polyval(p,x1);
            plot(x1,y1,strcat('-',theColor));
            axis([-1.0,1.0,0.0,1]); hold on;
            bins100=-1050:100:1050;
            mTF = histcounts(relRT(selection),bins100);
            mT = histcounts(relRT(selection&correct==1),bins100);
            display(strcat('number of keystrokes: ',num2str(sum(mTF))));
            %plot(-1:0.1:1,mT./mTF,strcat('-',theMarker));
            probs = mT./mTF;
            errs = 2*sqrt(probs.*(1-probs)./mTF);
            %errorbar(-1:0.1:1,probs,errs,strcat(theMarker,theColor));
            plot(-1:0.1:1,probs,strcat(theMarker,theColor));
            %plot(-1:0.1:1,mTF./sum(mTF),strcat('-',theColor));
            plot(-1:0.1:1,errs,strcat('-*',theColor));
            grid on;
            
        case 'relaccVrtFishType'
            constraints = varargin{3};
                     runDemo(fpdata,'relaccVrt',constraints&visual==8&audio==8,'+','r');
            hold on; runDemo(fpdata,'relaccVrt',constraints&visual==5&audio==5,'x','b');
            hold on; runDemo(fpdata,'relaccVrt',constraints&visual==8&audio==5,'d','m');
            hold on; runDemo(fpdata,'relaccVrt',constraints&visual==5&audio==8,'s','k');
            hold on; runDemo(fpdata,'relaccVrt',constraints&visual==0,'o','b');
            hold on; runDemo(fpdata,'relaccVrt',constraints&audio==0,'*','c');
            legend('FF fit','FF','FF err','SS fit','SS','SS err','FS fit','FS','FS err','SF fit','SF','SF err','Hidden fit','Hidden','Hidden err','Silent fit','Silent','Silent err');
            hold off;
            
        case 'relaccVrtFishTypeSmooth'
            constraints = varargin{3};
                     runDemo(fpdata,'relaccVrtSmooth',constraints&visual==8&audio==8,'+','r');
            hold on; runDemo(fpdata,'relaccVrtSmooth',constraints&visual==5&audio==5,'x','b');
            hold on; runDemo(fpdata,'relaccVrtSmooth',constraints&visual==8&audio==5,'d','m');
            hold on; runDemo(fpdata,'relaccVrtSmooth',constraints&visual==5&audio==8,'s','k');
            hold on; runDemo(fpdata,'relaccVrtSmooth',constraints&visual==0,'o','g');
            hold on; runDemo(fpdata,'relaccVrtSmooth',constraints&audio==0,'*','c');
            legend('FF fit','FF','FF err','SS fit','SS','SS err','FS fit','FS','FS err','SF fit','SF','SF err','Hidden fit','Hidden','Hidden err','Silent fit','Silent','Silent err');
            hold off;
            
            
         case 'accVrtSmooth'
            constraints = varargin{3};
            theMarker = varargin{4};
            theColor = varargin{5};
            
            selection = constraints & eventType==2;
            x = reaction(selection)/1000;
            y = correct(selection);
            smoothT=[];
            smoothF=[];
            smooth=[];
            for i=[200:1400]
                smoothT(i)=sum(selection&reaction>i-50&reaction<i+50&correct==1);
                smoothF(i)=sum(selection&reaction>i-50&reaction<i+50&correct==0);
                smooth(i) = sum(selection&reaction>i-50&reaction<i+50);
            end
            N = sum(smooth)/100;
            smoothT = smoothT.*(smooth./N>0.01);
            
            vals = (1:1400)./1000;
            %lot(vals,smoothT/N,strcat('-',theColor));hold on;
            %plot(vals,smoothF/N,strcat('.',theColor));
            plot(vals,smooth/N,strcat('.',theColor));hold on
            plot(vals,smoothT./smooth,strcat('-',theColor));
            %plot(vals,smooth/N,strcat('.-',theColor));
            %return;
            x1=0.4:0.01:1.4;
            p = polyfit(x,y,6);y1=polyval(p,x1);
            %plot(x1,y1,strcat('-',theColor));
            axis([0,2,0.0,1]); 
            hold on;
            bins100=350:100:1450;
            mTF = histcounts(reaction(selection),bins100);
            mT = histcounts(reaction(selection&correct==1),bins100);
            display(strcat('number of keystrokes: ',num2str(sum(mTF))));
            %plot(0.6:0.1:1.2,mT./mTF,strcat('-',theMarker));
            p = mT./mTF;
            errs = 2*sqrt(p.*(1-p)./mTF);
            %errorbar(0.4:0.1:1.4,p,errs,strcat(theMarker,theColor));
            plot(0.4:0.1:1.4,p,strcat(theMarker,theColor));
            %plot(0.4:0.1:1.4,mTF./sum(mTF),strcat(theMarker,theColor));
            %plot(0.4:0.1:1.4,errs,strcat('*',theColor));
            grid on;
            
            
        case 'relaccVrtSmooth'
            constraints = varargin{3};
            theMarker = varargin{4};
            theColor = varargin{5};
            
            selection = constraints & eventType==2;
            x = relRT(selection); %/1000;
            y = correct(selection);
            
            smoothT=[];
            smoothF=[];
            smooth=[];
            for i=[-1000:1000]
                smoothT(i+1001)=sum(selection&relRT>i-50&relRT<i+50&correct==1);
                smoothF(i+1001)=sum(selection&relRT>i-50&relRT<i+50&correct==0);
                smooth(i+1001) = sum(selection&relRT>i-50&relRT<i+50);
            end
            N = sum(smooth)/100;
            smoothT = smoothT.*(smooth./N>0.01);
            
            vals = (-1000:1000)./1000;
            %lot(vals,smoothT/N,strcat('-',theColor));hold on;
            %plot(vals,smoothF/N,strcat('.',theColor));
            plot(vals,smooth/N,strcat('.',theColor));hold on
            plot(vals,smoothT./smooth,strcat('-',theColor));
            %return;
            %x1=-0.4:0.01:0.4;
            %p = polyfit(x,y,6);y1=polyval(p,x1);
            %plot(x1,y1,strcat('-',theColor));
            axis([-1.0,1.0,0.0,1]); hold on;
            bins100=-1050:100:1050;
            mTF = histcounts(relRT(selection),bins100);
            mT = histcounts(relRT(selection&correct==1),bins100);
            display(strcat('number of keystrokes: ',num2str(sum(mTF))));
            %plot(-1:0.1:1,mT./mTF,strcat('-',theMarker));
            probs = mT./mTF;
            errs = 2*sqrt(probs.*(1-probs)./mTF);
            %errorbar(-1:0.1:1,probs,errs,strcat(theMarker,theColor));
            plot(-1:0.1:1,probs,strcat(theMarker,theColor));
            %plot(-1:0.1:1,mTF./sum(mTF),strcat('-',theColor));
            %plot(-1:0.1:1,errs,strcat('-*',theColor));
            grid on;
            
            
        case 'showlevels'
            subplot(2,1,1);histogram(vlevel,0.5:1:10.5);title('highest level attained in visual mode');
            subplot(2,1,2); histogram(alevel,0.5:1:10.5);title('highest level attained in auditory mode');
            axis([0,12,0,14]);
            visualLevel = vlevel'; auditoryLevel=alevel'; userids=[1:57]';
            display(table(userids,visualLevel,auditoryLevel));
            display(table(users',vlevel(users)',alevel(users)'));
            
        case 'calcbias'
            urt =[];
            usd=[];
            u7 = [];
            for u=users7
                u7 = [u7; u];
                urt0 =  mean(reaction(eventType==2&userId==u&mode==1));
                usd0 = std(reaction(eventType==2&userId==u&mode==1));
                urt = [ urt; urt0];
                usd = [usd; usd0];
            end
            %display([u7,urt,usd])
            tab = table(u7,urt,usd);
            display(tab);
            scatter(urt',usd');
            
      case 'calcbias2'
          % calculate bestfit line for correctness vs reaction time for
          % visual mode
            subplot(1,1,1);
            urt =[];
            usd=[];
            u7 = [];
            counter=1;
            selection = eventType==2 & mode==1 & reaction>=400 & reaction <=1200 & visual==8;
            for u=users7
                subplot(5,5,counter);
                scatter(reaction(selection&userId==u),100.*correct(selection&userId==u));
                title(u);
                lsline;
                counter = counter + 1;
                axis([0,2000,0,200]);
                s0=sum(selection&userId==u&correct==1)
                hold on;
                histogram(reaction(selection&userId==u&correct==1))
                histogram(reaction(selection&userId==u&correct==0))
                %hold on;
                hold off;
            end
            hold off;
       case 'calcbias3'
          % calculate bestfit quadratic for correctness vs reaction time for
          % visual mode
            subplot(1,1,1);
            urt =[];
            usd=[];
            u7 = [];
            counter=1;
            selection = eventType==2 & level>3 & mode==1 & reaction>=400 & reaction <=1200 & visual==8;
            for u=users7
                subplot(5,5,counter);
                x=reaction(selection&userId==u);
                y = 100.*correct(selection&userId==u);
                scatter(x,y);
                title(u);
                p = polyfit(x,y,2);
                x1 = linspace(400,1200,17);
                y1 = polyval(p,x1);
                hold on
                plot(x1,y1);
                counter = counter + 1;
                axis([0,2000,0,200]);
                s0=sum(selection&userId==u&correct==1);
                %hold on;
                histogram(reaction(selection&userId==u&correct==1),0:50:2000)
                histogram(reaction(selection&userId==u&correct==0),0:50:2000)
                %hold on;
                hold off;
            end
            hold off;
            
      case 'calcbias4'
          % update reactions times by subtracting the mean reaction time
          % for each individual and then plotting a histogram of the
          % correct and incorrect responses wrt modified reaction time
          % 
            subplot(1,1,1);
            muid = max(userId);
            urt =zeros(1,muid);;
            usd=[];
            u7 = userId == 0;
            counter=1;
            
            selection = eventType==2 & mode==1 & reaction>=400 & reaction <=1200&visual==8;

            for u=users7
                urt(u) = mean(reaction(selection&userId==u));
                u7 = u7 | userId==u;
            end
            
            m = mean(reaction(selection&u7));
            %histogram(reaction(selection&u7)- m,-1000:50:1000);hold on;
            display(std(reaction(selection&u7)-m));
            %display(size(reaction(selection)));
            %display(size(userId(selection)));
            %display(size(urt(userId(selection))));
            modreact = reaction(selection&u7) - urt(userId(selection&u7))';
            display(std(modreact));
            modreactT = reaction(selection&correct==1&u7) - urt(userId(selection&correct==1&u7))';
            modreactF = reaction(selection&correct==0&u7) - urt(userId(selection&correct==0&u7))';
            histogram(modreact,-1000:10:1000); hold on;
            histogram(modreactT,-1000:10:1000);
            histogram(modreactF,-1000:10:1000);      
            hold off;
            
      case 'calcbias5'
          % find the best fit quadratics for accuracy wrt normalized reaction
          % times as above. Use
          % update reactions times by subtracting the mean reaction time
          % for each individual and then plotting a histogram of the
          % correct and incorrect responses wrt modified reaction time
          % 
            subplot(1,1,1);
            muid = max(userId);
            urt =zeros(1,muid);;
            usd=[];
            u7 = userId == 0;
            counter=1;
            
            selection = eventType==2 & mode==1 & visual==5 & reaction>=400 & level<=7 & level >=4 & reaction <=1200;

            for u=users7
                urt(u) = mean(reaction(selection&userId==u));
                u7 = u7 | userId==u;
            end
            
            m = mean(reaction(selection&u7));
            
            
            modreact  = reaction(selection&u7) - urt(userId(selection&u7))';
            modreactT = reaction(selection&correct==1&u7) - urt(userId(selection&correct==1&u7))';
            modreactF = reaction(selection&correct==0&u7) - urt(userId(selection&correct==0&u7))';
            histogram(modreact,-1000:10:1000); hold on;
            histogram(modreactT,-1000:10:1000); hold on
            histogram(modreactF,-1000:10:1000); 
            mT = histcounts(modreactT,-425:50:425);
            mF = histcounts(modreactF,-425:50:425);
            mTF = histcounts(modreact,-425:50:425);
            plot(-400:50:400,mT./mTF*100,'-*b');
     
            x=modreact;
            y = 100*correct(selection&u7);
            yneg = y(x<0);
            xneg = x(x<0);
            scatter(x,y);
            p = polyfit(xneg,yneg,2);
            x1 = -400:50:0;
            y1 = polyval(p,x1);
            hold on
            plot(x1,y1,'-*r');
            xpos = x(x>0);
            ypos = y(x>0);
            ppos = polyfit(xpos,ypos,2);
            x2 = 0:50:400;
            y2 = polyval(ppos,x2);
            plot(x2,y2,'-*r');
            %axis([-1000,1000,0,200]);   
            
            title(['mean=',num2str(m),'visual']);
            hold off;
            
      case 'gameplay1'  %runDemo(fpdata,'gameplay1')
            subplot(1,1,1);
            hold off;
            scatter(rowTime,levelTime,10,userId,'filled');
            grid on
      case 'gameplay2'
            subplot(1,1,1);
            hold off;
            uid = varargin{3}
            scatter(rowTime(userId==uid),levelTime(userId==uid),10,'b','filled');
        case 'gameplay3'
            uid = varargin{3}
            rad=15;
            scatter(rowTime(userId==uid&correct==1),levelTime(userId==uid&correct==1),rad,'b','filled');
            hold on ; scatter(rowTime(userId==uid&correct==0),levelTime(userId==uid&correct==0),rad,'r','filled');
            grid on ;set(gca,'YTick',0:500:18000); 
            set(gca,'XTick',min(rowTime):60*1000:max(rowTime))
            hold off
        case 'gameplay4'
            uid = varargin{3};
            rad=10;
            filter = userId==uid&correct==1;
            scatter(rowTime(filter&mode==1),levelTime(filter&mode==1),rad,'b','+');%'filled');
            hold on;
            scatter(rowTime(filter&mode==2),levelTime(filter&mode==2),rad,'g','+'); %'filled');
            hold on ; scatter(rowTime(userId==uid&correct==0),levelTime(userId==uid&correct==0),rad,'r','+'); %'filled');
            grid on ;set(gca,'YTick',0:500:18000); 
            set(gca,'XTick',min(rowTime):60*1000:max(rowTime));
            hold off;
        case 'gameplay4a'
            uid = varargin{3};
            rad=20;
            filter = userId==uid&correct==1;
            scatter([20],[10],10,'r','o');
            hold on;
            scatter(rowTimeN(filter&mode==1),levelTimeN(filter&mode==1),rad,'b','+');%'filled');
            hold on;
            scatter(rowTimeN(filter&mode==2),levelTimeN(filter&mode==2),rad,'g','+'); %'filled');
            hold on ; scatter(rowTimeN(userId==uid&correct==0),levelTimeN(userId==uid&correct==0),rad,'r','+'); %'filled');
            grid on ;
            %set(gca,'YTick',0:2000:20000); 
            set(gca,'YTick',0:1:10);
            %set(gca,'XTick',min(rowTimeN):5:max(rowTimeN));
            set(gca,'XTick',0:5:45);
            %set(gca,'XTick',0:5:50);
            axis([0,45,0,10]);
            grid minor;
            hold off;
        case 'gameplay5'
            [a,b]=histcounts(userId,1:1000);
            userids = 1:999;
            userids = userids(a>0);
            m = ceil(sqrt(max(userids)));
            for u=userids
                subplot(m,m,u);
                runDemo(fpdata,'gameplay4a',u);
            end
        case 'accuracy1'
            m = varargin{3}; % this is the mode 1=visual 2=auditory
            acong = correct(eventType==1& mode==m & action==1 & congruent==1&level>0);
            aincong = correct(eventType==1& mode==m & action==1 & congruent==2 & ~oddball&level>0);
            aodd = correct(eventType==1& mode==m & action==1 & oddball&level>0);
            rcong = reaction(eventType==1& mode==m & correct==1 & action==1 & congruent==1&level>0);
            rincong = reaction(eventType==1& mode==m & correct==1 & action==1 & congruent==2 & ~oddball&level>0);
            rodd = reaction(eventType==1& mode==m & correct==1 & action==1 & oddball&level>0);

            condition = {'congruent';'incongruent';'oddball'};
            ameans = [mean(acong);mean(aincong);mean(aodd)];
            aN = [length(acong);length(aincong);length(aodd)];
            
            rmeans = [mean(rcong);mean(rincong);mean(rodd)];
            rstdevs = [std(rcong);std(rincong);std(rodd)];
            rN = [length(rcong);length(rincong);length(rodd)];
            T = table(condition,ameans,aN, rmeans,rstdevs,rN);
            disp(T)





        case 'rt1'
            cong = reaction(action==1 & congruent==1);
            incong = reaction(action==1 & congruent==2 & ~oddball);
            odd = reaction(action==1 & oddball);
            disp({'congruent',mean(cong),std(cong),length(cong)});
            disp({'incongruent',mean(incong),std(incong),length(incong)});
            disp({'oddball',mean(odd),std(odd),length(odd)});
            condition = {'congruent';'incongruent';'oddball'}
            means = [mean(cong);mean(incong);mean(odd)];
            stdevs = [std(cong);std(incong);std(odd)];
            N = [length(cong);length(incong);length(odd)];
            T = table(condition,means,stdevs,N);
            disp(T)
        case 'rt2'
            subplot(1,1,1);
            hold off;
            for lev=0:9
                subplot(10,1,lev+1);
                [x1,y1] = histcounts(reaction(level==lev&correct==1&mode==1),0:50:2000);
                plot(y1,[0,x1]./sum(level==lev&correct==1&mode==1),'-b');
                hold on;
                [x1,y1] = histcounts(reaction(level==lev&correct==1&mode==2),0:50:2000);
                plot(y1,[0,x1]./sum(level==lev&correct==1&mode==2),'-r');
                grid on;
                %hold on
            end
        case 'rt3'  %plot reaction times by level
            subplot(1,1,1);
            hold off;
            for lev = 0:9
                subplot(10,1,lev+1);
                h = histogram(reaction(level==lev&correct==1),0:25:2000);
                h.Normalization = 'probability'
            end
        case 'rt4'  %plot reaction times by user
            subplot(1,1,1);
            hold off;
            for uid = 3:14
                subplot(12,1,uid-2);
                histogram(reaction(userId==uid&correct==1),0:25:2000);
            end
        case 'rt5'  % plots reaction times by user and level
            subplot(1,1,1);
            m = varargin{3}
            hold off;
            for lev = 0:9
                for uid = 0:11
                    hold off;
                    subplot(10,12,12*lev+uid+1);
                    histogram(reaction(userId==uid+3&mode==m&level==lev&correct==1),0:50:2000);
                    grid on
                end
            end
        case 'rt6' % reaction time for correct vs incorrect responses
            subplot(1,1,1);
            hold off;
            h1 = histogram(reaction(correct==1&action==1),0:50:2000);
            h1.Normalization= 'probability';
            hold on
            h2 = histogram(reaction(correct==0&action==1),0:50:2000);
            h2.Normalization= 'probability';
        case 'rt7' % proportion correct for each reaction time
            subplot(1,1,1);
            hold off;
            perCor=[];
            for u=0:19
                corr = reaction<(u+1)*100 & reaction >= u*100& action==1;
                num = sum(corr);
                numc = sum(corr & correct==1);
                perCor(u+1)=100*numc/num;
            end            
            plot(50:100:2000,perCor,'-g');
            hold on;
            num=[];
            numc=[];
            tot = sum(action==1);
            for u=0:19
                corr = reaction<(u+1)*100 & reaction >= u*100& action==1;
                num(u+1) = sum(corr&correct==1)/tot*100;
                numc(u+1) = sum(corr & correct==0)/tot*100;
            end
            plot(50:100:2000,num,'-xb')
            hold on;
            plot(50:100:2000,numc,'-xr');
            grid on;
            set(gca,'YTick',0:10:100)
            set(gca,'XTick',0:100:2000)
        case 'rt8'
            % show correctness and reaction time
            % for fast/fast, slow/slow, fast/slow, slow/fast
            % and level > 0
            lev = varargin{3};
            m = varargin{4};
            ffa = correct(visual==8 & audio==8 & level>=lev & mode==m);
            ssa = correct(visual==5 & audio==5 & level>=lev & mode==m);
            sfa = correct(visual==5 & audio==8 & level>=lev & mode==m);
            fsa = correct(visual==8 & audio==5 & level>=lev & mode==m);
            
            ffr = reaction(correct==1 & visual==8 & audio==8 & level>=lev & mode==m);
            ssr = reaction(correct==1 & visual==5 & audio==5 & level>=lev & mode==m);
            sfr = reaction(correct==1 & visual==5 & audio==8 & level>=lev & mode==m);
            fsr = reaction(correct==1 & visual==8 & audio==5 & level>=lev & mode==m);

            
            ff = [mean(ffa),mean(ffr),std(ffr)]';
            ss = [mean(ssa),mean(ssr),std(ssr)]';
            fs = [mean(fsa),mean(fsr),std(fsr)]';
            sf = [mean(sfa),mean(sfr),std(sfr)]';
            t = table(ff,ss,fs,sf,'RowNames',{'accuracy(mean)','rt(mean)','rt(stdev)'});
            disp(t);
            
                
        case 'learning'
            % for one user, look at how their accuracy and RT change by
            % game
            u = varargin{3};  % userId number
            m = varargin{4};  % mode visual=1 audio=2
            % get all of the game numbers of games played by the user
            [a,b] = histcounts(gameId(userId==u&mode==m));
            ugames = b(a>0); % this is the list of gamenumbers
            %histogram(ugames,1:1500)
            %display(ugames)
            rt=[];
            for g = 1:length(ugames)
                rt(g) = mean(reaction(gameId==ugames(g)&eventType==2&correct==1));
            end
            %display(rt);
            plot(rt);
            
        case 'learning2'
            m = varargin{3}; % mode
            hold off; 
            for u=users
                runDemo(fpdata,'learning',u,m); 
                hold on;
            end
            
        case 'timeCorrect'
            % here we calculate the percentage of responses at a given
            % speed that are correct. We expect that correctness rises as
            % the subject takes longer to respond
            
            pcorrect=[]; total=[];
            for i=1:41;
                keypresses =  eventType==2 & reaction >= 50*(i-1) & reaction < 50*i;
                pcorrect(i) = sum(keypresses&correct==1)/sum(keypresses);
                total(i) = sum(keypresses);
                if ((i<8) || (i>24)) pcorrect(i)=0; end
            end
            display(pcorrect);
            %plot(0:50:2000, pcorrect);
            plot(0:50:2000, pcorrect.*(total>10)); hold on;
            plot(0:50:2000, total/sum(total)); hold off;
            
        case 'timeCorrect2'
            % here we calculate the percentage of responses at a given
            % speed that are correct for a give user u.
            % We expect that correctness rises as
            % the subject takes longer to respond
            u = varargin{3}; % userId
            
            pcorrect=[];
            total = [];
            minBin=400; maxBin=1200; binSize=50;
            for i=1:21;
                keypresses = visual==8 & audio==8 & userId==u & eventType==2 & reaction >= 100*(i-1) & reaction < 100*i;
                pcorrect(i) = sum(keypresses&correct==1)/sum(keypresses);
                total(i) = sum(keypresses);
            end
            %display(pcorrect);
            density = total/sum(total);
            %display(total);
            %display(sum(total));
            %display([length(density),length(pcorrect)]);
            plot(0:100:2000, pcorrect.*(total>10)); hold on;
            
            plot(0:100:2000, total/sum(total)); hold off;
            
            
        case 'timeCorrect2a'
            % here we calculate the percentage of responses at a given
            % speed that are correct for a give user u.
            % We expect that correctness rises as
            % the subject takes longer to respond
            u = varargin{3}; % userId
            
            pcorrect=[];
            total = [];
            minBin=-400; maxBin=400; binSize=100;
            keypresses = userId==u & eventType==2 & level>=4 &  mode==1 & audio==0;
            bias = mean(reaction(keypresses & correct==1 ));
            for i=1:(maxBin-minBin)/binSize+1;
                
                bin = bias + minBin + binSize*(i-1);
                
                ukeypresses =  keypresses & reaction >= bin-binSize/2 & reaction < bin+binSize/2;
                pcorrect(i) = sum(ukeypresses&correct==1)/sum(ukeypresses);
                total(i) = sum(ukeypresses & correct==1);
            end
            %display(pcorrect);
            %display([u,pcorrect]);
            density = total/sum(total);
            %display(total);
            %display(sum(total));
            %display([length(density),length(pcorrect)]);
            s=sum(total);
            plot(minBin:binSize:maxBin, pcorrect.*(density>0.05).*(total>=3),'*'); hold on;
            
            plot(minBin:binSize:maxBin, total/sum(total),'-*'); hold off;
            out=[u,bias, pcorrect];
            
            
            
        case 'timeCorrect3'
            % this plots the RT/ProbCorrect curves for each user
            for u=users
                subplot(8,8,u);
                runDemo(fpdata,'timeCorrect2',u);
                %subplot(16,8,u+64);
                %runDemo(fpdata,'gameplay4',u);
            end
            

        case 'timeCorrectu7'
            % this plot shows the individual variation in accuracy
            % in terms of the reaction time, relative to the persons
            % mean reaction time.
            % we can tune timeCorrect2a to look at particular classes of
            % events
            subplot(2,2,1); subplot(1,1,1);
            ptable=[];
            for u=users7;
                hold on; pcorrect= runDemo(fpdata,'timeCorrect2a',u); 
                %display(pcorrect');
                ptable=[ptable;pcorrect]; 
            end
            z=sortrows(ptable,[2]);
            %display(ptable);
            display(z);
            display(mean(ptable));
            hold off;
            title('Accuracy wrt RT relative to subject''s mean RT: visual, level>=4');
            
            
            
            
            
                
        

    end






end

