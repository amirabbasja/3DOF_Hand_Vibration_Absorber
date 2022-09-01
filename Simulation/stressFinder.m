function [stress] = stressFinder(prevPhase, strain, strainDot, T, eL, Gm, Ga, z0, zT0, S_M_finish, S_M_start, S_A_finish, S_A_start)
    syms q;
    fprintf("stressFinder; z0:%.5f ||strain:%.5f ||prevPhase:%.1f || \n", z0,strain,prevPhase);
    %In this Function we try to get the stress with strain as an input. the
    %logic is: First we assume we are not in transformation phase and try
    %to solve the equation, if no answers were found, our assumption was
    %wrong and we are in transformation phase
    
    %z0 is the z at the start of new phase, but z is the previouse step's
    %martnsite ratio. this comes handy when we are decreasing the strain
    zS0 = z0 - zT0;
    
    %Defining the phase transformation as function handlers
    %1.Martensite percentage for austenite to martensite
    Z_A2M = @(stress) ((1-zS0)/2*cos(pi/(S_M_start - S_M_finish) * (stress - S_M_finish)) + (1 + zS0)/2);
    %2.Martensite percentage for martensite to austenite
    Z_M2A = @(stress) (z0/2*(cos( pi/(S_A_start-S_A_finish)*(stress-S_A_start) )+1));

    %UPDATE: some times, we have negative strain. In order for this
    %function to handle that, we define This constant. if it is true, all
    %the calculations will be in positive numbers but at the end they will
    %be multiplied by -1
    NEGATIVESTRAIN = false;
    
    if(strain < 0)
        NEGATIVESTRAIN = true;
        strain = strain * -1;
    end

    if 0 <= strainDot
        % If stress/strain is increasing search like this
        %1.Assume no transformation
        if z0 == 0
            s = vpasolve( q == Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0) , q , [0 S_M_start]);
        elseif z0 == 1
            s = vpasolve( q == Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0) , q , [S_M_finish inf]);
        elseif 0 < z0 && z0 < 1
            s = vpasolve( q == Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0) , q , [0 S_M_start]);
        else
            fprintf("ERR, z0 = %.5f\n",z0)
        end

        
        if(size(s,1) == 0)
            fprintf("No answers found considering no transformation!\n")
            %2.No answers found, Try phase transformation

            % **Note**: I dont know why, but when we use (Func_SMA_Shear_Modulus(Z_A2M(q),Gm,Ga,1)^-1)
            % instead of (Z_M2A(q)/Gm + (1-Z_M2A(q))/Ga) in the relation below, vpasolve has some
            % convergance issues.
            s = vpasolve( strain == ((Z_A2M(q)/Gm + (1-Z_A2M(q))/Ga)) * q + eL * Z_A2M(q) , q, [S_M_start S_M_finish]);
    
            if(size(s,1) == 0)
                fprintf("No answers found considering transformation!\n -> The inputs are: prevPhase=%.1f ||  strain=%.5f ||  strainDot=%.1f ||  z0=%.5f ||\n ",prevPhase, strain, strainDot,z0)
                %No answers found. Some times this occurs when we finish phase
                % transformation and continue lineraly in the same
                %direction. At the point of transition two problems occur. 
                % first z never reaches exactly to 1 (despite how close it may be) 
                %,so we cant really know when the transformation ends. 
                % but when iteration  through e(i)'s we reach a point where no
                %answer is found from vpasolve. This is where we know we
                %have passed the transformation level. second is that we
                %have to update z0 to continue iteration, but we have to
                %do it without breaking the iteration loop. so we use
                %recursion to solve this matter.
                if(prevPhase == 2 || prevPhase == 222 )
                    fprintf("ERR1\n")
                    stress = stressFinder(prevPhase, strain, strainDot,T,eL,Gm,Ga,1,zT0, S_M_finish, S_M_start, S_A_finish, S_A_start);
                else
                    fprintf("ERR2\n")
                    dummy=1
                end
        
            elseif(size(s,1) == 1)
                %Only one answer found
                fprintf("One answer found by considering Transformation: stress = %.2f\n",s(1))
                stress = s(1);
            else
                %More than one answer found! 
                %Print an Error message!
                stress = -1;
            end
    
        elseif(size(s,1) == 1)
            %Only one answer found
            fprintf("One answer found by considering no Transformation: stress = %.2f\n",s(1))
            stress = s(1);
        else
            %More than one answer found! 
            %Print an Error message!
            stress = -1;
        end



    elseif strainDot < 0
        % If stress/strain is decreasing search like this
        %1.Assume no transformation
        if z0 == 0
            s = vpasolve( q == Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0) , q , [0 S_A_finish]);
        elseif z0 == 1
            s = vpasolve( q == Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0) , q , [S_A_start inf]);
        elseif 0 < z0 && z0 < 1
            s = vpasolve( q == Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0) , q , [S_A_start inf]);
        else
            fprintf("ERR, z0 = %.5f\n",z0)
        end
    
        if(size(s,1) == 0)
            fprintf("No answers found considering no transformation!\n")
            %2.No answers found, Try phase transformation
    
            % **Note**: I dont know why, but when we use (Func_SMA_Shear_Modulus(Z_A2M(q),Gm,Ga,1)^-1)
            % instead of (Z_M2A(q)/Gm + (1-Z_M2A(q))/Ga) in the relation below, vpasolve has some
            % convergance issues. 
            s = vpasolve( strain == ((Z_M2A(q)/Gm + (1-Z_M2A(q))/Ga)) * q + eL * Z_M2A(q) , q, [S_A_finish S_A_start]);
    
            if(size(s,1) == 0)
                fprintf("No answers found considering transformation!\n -> The inputs are: prevPhase=%.1f ||  strain=%.5f ||  strainDot=%.1f ||  z0=%.5f ||\n ",prevPhase, strain, strainDot,z0)
                %No answers found. Some times this occurs when we finish phase
                % transformation and continue lineraly in the same
                %direction. At the point of transition two problems occur. 
                % first z never reaches exactly to 1 (despite how close it may be) 
                %,so we cant really know when the transformation ends. 
                % but when iteration  through e(i)'s we reach a point where no
                %answer is found from vpasolve. This is where we know we
                %have passed the transformation level. second is that we
                %have to update z0 to continue iteration, but we have to
                %do it without breaking the iteration loop. so we use
                %recursion to solve this matter.
                if(prevPhase == 5)
                    stress = stressFinder(prevPhase, strain, strainDot,T,eL,Gm,Ga,0,zT0, S_M_finish, S_M_start, S_A_finish, S_A_start);
                else
                    dummy=1;
                end
        
            elseif(size(s,1) == 1)
                %Only one answer found
                fprintf("One answer found by considering Transformation: stress = %.2f\n",s(1))
                stress = s(1);

            else
                %More than one answer found! 
                %Print an Error message!
                fprintf("ERR3\n")
                stress = -1;
            end
    
        elseif(size(s,1) == 1)
            %Only one answer found
            fprintf("One answer found by considering no Transformation: stress = %.2f\n",s(1))
            stress = s(1);
        else
            %More than one answer found! 
            %Print an Error message!
            fprintf("ERR4\n")
            stress = -1;
        end

    end
    
    if NEGATIVESTRAIN
        stress = stress * -1;
    end
    
    stress = double(stress);
end