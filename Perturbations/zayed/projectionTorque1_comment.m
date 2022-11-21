function [T]=projectionTorque1a(G,vorbit)
% les 4 points qui delimitent le cadre 1
A=[0.17,0.11,0.05];%ok
B=[0.17,0.28,0.05];
C=[-0.17,0.28,0.05];
D=[-0.17,0.11,0.05];%ok


torque=[0,0,0];


%projection of B
pBy=0.11;
t=(pBy-B(2))/vorbit(2);
pBx=B(1)+t*vorbit(1);
pBz=B(3)+t*vorbit(3);

%projection of C
pCy=pBy;
pCx=C(1)+t*vorbit(1);%% because they have the same y
pCz=pBz;%%because of //

%vecteur directeur de notre droite
ux=A(1)-pBx;
uy=A(2)-pBy;%% zero we are in zx plan
uz=A(3)-pBz;
if pCz>-0.05                    %IN THE DIRECTION OF Z SHADOWING +Y FACE LIMIT BOTTOM
    %display("pCz>-0.05, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    if pBx<0.17                 %IN THE DIRECTION OF X SHADOWING +Y FACE LIMIT BACK
        %display("  pBx<0.17, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        if pBx>-0.17 %%we caculate(frcement pCx<-0.17)
                                %IN THE DIRECTION OF X SHADOWING +Y FACE
                                %LIMIT FRONT
                 %NOW WE ARE CREATING A SHADOW ON +Y LATERAL FACE
            %display("    pBx>-0.17, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            %pas
            N=100;%%%%%%%%%%%%%%%%%%%can be changed
            pas=ux/N;
            %integrale de rieman
            for i=0:N-1
                xr=pBx+i*pas;
                t=(xr-A(1))/ux;
                z=A(3)+t*uz;
                %geometric center for dA
                Gc=[xr+(pas/2),pBy,(z-0.05)/2];
                Rs=Gc-G;
                torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
            end
            %the other surface
            Gc=[(pBx-0.17)/2,pBy,(pBz-0.05)/2];
            Rs=Gc-G;
            torque=torque+(pBx+0.17)*(pBz+0.05)*cross(Rs,vorbit);
        else
            %display("    pBx<=-0.17, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            %%we starte (la droite (ApB) coupe [segment]
            N=200;
            pas=0.34/N;
            
            for i=0:N-1
                xr=-0.17+i*pas;
                t=(xr-A(1))/ux;
                z=A(3)+t*uz;
                Gc=[xr+(pas/2),pBy,(z-0.05)/2];
                Rs=Gc-G;
                torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
                
            end
            
            
            
        end
        
    elseif pCx<0.17 %%we konw already that pCx>-0.17
        %display("  pBx>=0.17 and pCx<0.17, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        %we calculate
        ux=D(1)-pCx;
        uy=D(2)-pCy;%%zero
        uz=D(3)-pCz;
        N=100;
        pas=ux/N;%negative always
        %%%%%%
        %S=((pCx+0.17)*(0.1+pCz+0.05))/2 +(0.17-pCx)*(pCz+0.05);
        %display(["pCx",pCx])
        %display(["pCz",pCz])
        %display(["calcul surface non cachée",S])

        for i=0:N-1
            xr=pCx+i*pas;
            t=(xr-D(1))/ux;
            z=D(3)+t*uz;
            %geometric center for dA
            Gc=[xr+(pas/2),pCy,(z-0.05)/2];
            Rs=Gc-G;
            torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
            
        end
        Gc=[(pCx+0.17)/2,pCy,(pCz-0.05)/2];
        Rs=Gc-G;
        torque=torque+(0.17-pCx)*(pCz+0.05)*cross(Rs,vorbit);
    else
        %display("  pBx>=0.17 and pCx>=0.17, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        ux=D(1)-pCx;
        uy=D(2)-pCy;
        uz=D(3)-pCz;
        %we calculate(coupe droite segment)
        N=200;
        pas=-0.34/N;
     
        for i=0:N-1
            xr=0.17+i*pas;
            t=(xr-D(1))/ux;
            z=D(3)+t*uz;
            Gc=[xr+(pas/2),pCy,(z-0.05)/2];
            Rs=Gc-G;
            torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
            
        end
        
    end
else
    %display("pCz<=-0.05, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    %let Dprim Aprim to be like
    Dprim=[-0.17,0.11,-0.05];
    Aprim=[0.17,0.11,-0.05];
    if pBx<0.17
        %display("  pBx<0.17, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        ux=A(1)-pBx;
        uy=A(2)-pBy;
        uz=A(3)-pBz;
        %%%%
        ADpv=A-Dprim;
        %%%%%
        tprim=(pBz-A(3))/ADpv(3);
        %on compare
        xc=A(1)+tprim*ADpv(1);
        if pBx<xc%%%derriere suivant les x negative
            %display("    pBx<xc, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            N=200;%%%%%%%%%%%%%%%%%%%can be changed
            pas=0.34/N;
            %integrale de rieman
            
            for i=0:N-1
                xr=-0.17+i*pas;%%%intersection DDprim
                t=(xr-A(1))/ux;
                z=A(3)+t*uz;
                %geometric center for dA
                Gc=[xr+(pas/2),pBy,(z-0.05)/2];
                Rs=Gc-G;
                torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
            end
            
        else
            %display("    pBx>=xc, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            %we calculate%%%%%%%%%%problem(we can use diagonale)
            t=(-0.05-A(3))/uz;%%intersection DprimAprim
            x=A(1)+t*ux;
            N=100;%%%%%%%%%%%%%%%%%%%can be changed
            pas=(0.17-x)/N;
            %integrale de rieman
          
            for i=0:N-1
                xr=x+i*pas;
                t=(xr-A(1))/ux;
                z=A(3)+t*uz;
                %geometric center for dA
                Gc=[xr+(pas/2),pBy,(z-0.05)/2];
                Rs=Gc-G;
                torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
            end
            
            
            
        end
        %%elseif pCx<0.17
        %we calculate
    else
        %display("  pBx>=0.17, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        ux=D(1)-pCx;
        uy=D(2)-pCy;
        uz=D(3)-pCz;
        %%%%
        DApv=D-Aprim;
        %%%%%
        tprim=(pCz-A(3))/DApv(3);
        %on compare
        xc=D(1)+tprim*DApv(1);
        
        if pCx>xc
            %display("    pCx>xc, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            N=200;%%%%%%%%%%%%%%%%%%%can be changed
            pas=-0.34/N;
            %integrale de rieman
            
            for i=0:N-1
                xr=0.17+i*pas;
                t=(xr-D(1))/ux;
                z=D(3)+t*uz;
                %geometric center for dA
                Gc=[xr+(pas/2),pCy,(z-0.05)/2];
                Rs=Gc-G;
                torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
            end
            
        else
            %display("    pCx<=xc, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            %we calculate%%%%%%%%%%problem(we can use diagonale)
            t=(-0.05-D(3))/uz;
            x=D(1)+t*ux;
            %display(["x",x])
            N=100;%%%%%%%%%%%%%%%%%%%can be changed
            pas=(-0.17-x)/N;

            %%%%%%%%%%%%
            %S=((x+0.17)*0.1)/2;
            %display(["calcul surface non cachée",S])

            %integrale de rieman
            for i=0:N-1
                xr=x+i*pas;
                t=(xr-D(1))/ux;
                z=D(3)+t*uz;
                %geometric center for dA
                Gc=[xr+(pas/2),pBy,(z-0.05)/2];
                Rs=Gc-G;
                torque= torque +cross(Rs,vorbit)*((z+0.05)*pas);
            end
            
            
            
        end
        
    end
    
end
T=torque;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
