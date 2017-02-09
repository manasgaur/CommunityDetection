function [js]=jaccardsimilarity(A,B)

c=size(A,2);
insec=0;union=0;
for i=1:c
    if (A(1,i)==1) && (B(1,i)==1)
        insec=insec+1;union=union+1;
    else
        if (A(1,i)==0) && (B(1,i)==0)
        continue;
        else
            union=union+1;
        end
    end
end
js=insec/union;
end