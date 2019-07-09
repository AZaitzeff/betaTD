function [replacemap,newdict,K]=postpilchank(neighbors,sizeneigh,dict,IQmean,filled,K,betatol)

    %add mean of the betas
    replacemap=zeros(1,K);
    newdict=dict;
    claimed=zeros(1,K,'logical');
    if betatol>.75
        betatol=betatol/180*pi;
    end

    [~,order]=sort(IQmean,'descend');
    z=1;
    for k=1:K
        ind=order(k);
            if filled(ind) && ~claimed(ind)
                
                lenstack=100;
                stack=zeros(1,lenstack);
                coder.varsize('stack',[1 2000],[0,1]);
                stack(1)=ind;
                pop=1;
                next=2;
                mainbeta=dict(:,ind);
                replacemap(ind)=z;
                newdict(:,z)=mainbeta;
                claimed(ind)=1;
                while pop<next
                    cur=stack(pop);
                    pop=pop+1;
                    ns=sizeneigh(cur);
                    for i=1:ns
                        curind=neighbors(cur,i);
                        if ~claimed(curind) && (~filled(curind) || b2bmetric(mainbeta,dict(:,curind))<betatol)
                            replacemap(curind)=z;
                            stack(next)=curind;
                            claimed(curind)=1;
                            next=next+1;
                            if next>lenstack
                                templenstack=ceil(lenstack*1.5);
                                tempstack=zeros(1,templenstack);
                                tempstack(1:lenstack)=stack;
                                lenstack=templenstack;
                                stack=tempstack;
                            end
                        end
                    end
                end
            z=z+1;    
                
            end
    end
    newdict(:,z:K)=[];
    K=z-1;
end


