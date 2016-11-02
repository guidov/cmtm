%This function detunes a core by using interp1 to linearly interpolate to a fixed number
%of control points.  Age, Dep, and Data all must agree, do not interpolate to an even age 1st.
%
%function [NAge] = detunePH(Age,Dep,Tfix)
%
%NAge   = New Age interpolated betweened fixed time periods Tfix.  
%Tfix   = the time step between control points, if blank defaults to 50*mean(diff(Age))

function [NAge] = detunePH(Age,Dep,Tfix)

     if length(Tfix)==0, Tfix=50*mean(diff(Age));  end;

     FAge = [Age(1)]; FDep = [Dep(1)];
     for t=Tfix:Tfix:Age(length(Age));
     place = find(Age>=t);
     FAge = [FAge Age(place(1))];
     FDep = [FDep Dep(place(1))]; 
     end;

     FAge(length(FAge)) = Age(length(Age));  %Make last fixed point the end of the data
     FDep(length(FDep)) = Dep(length(Dep));
     NAge = interp1(FDep,FAge,Dep);

display(['Control Points: ',num2str(length(FAge),2)])
display([FAge]);
