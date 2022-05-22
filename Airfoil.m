classdef Airfoil < handle 
    
    properties
        inputs; %Inputs to Constructor function
        
        b, e, beta, a, alpha, x0, y0, cosa, sina; % Calculated form Inputs
        
        delta_eta1, delta_eta2;  
        
        dx_deta1, dy_deta1, dx_deta2, dy_deta2; % Metrics of Transformation
        deta1_dx, deta1_dy, deta2_dx, deta2_dy;
        
        Jacobian, C11, C22, C12;
        
        C11_plusHalf_i, C11_negHalf_i, 
        C22_plusHalf_i, C22_negHalf_i,
        C12_plusHalf_i, C12_negHalf_i,
        
        C11_plusHalf_j, C11_negHalf_j, 
        C22_plusHalf_j, C22_negHalf_j,
        C12_plusHalf_j, C12_negHalf_j,
        
        s_i_j
        s_in1_j
        s_ip1_j
        s_i_jn1
        s_i_jp1
        s_in1_jn1, s_ip1_jp1
        s_in1_jp1, s_ip1_jn1
    end
    
    methods
        function obj = Airfoil(inputs_struct)
            obj.inputs = inputs_struct;
            obj.b = inputs_struct.chord/4;
            obj.e = inputs_struct.max_thickness/1.3;
            obj.beta = 2*inputs_struct.max_camber;
            obj.a = obj.b *(1+obj.e)/cos(obj.beta);
            obj.alpha = inputs_struct.alpha_deg * pi / 180; 
            obj.x0 = -obj.b * obj.e;
            obj.y0 = obj.a * sin(obj.beta);
            obj.cosa = cos(obj.alpha);
            obj.sina = sin(obj.alpha);
        end
        
        function [innerCircle, outerCircle, airfoil] = joukowskiAirfoil(this)
            theta = linspace(0, 2*pi, this.inputs.i_max);
            
            innerCircle = struct(...
                'x', this.inputs.chord/2*cos(theta),...
                'y', this.inputs.chord/2*sin(theta));
            
            outerCircle = struct(...
                'x', this.inputs.R*cos(theta),...
                'y', this.inputs.R*sin(theta));
            
            sign=(sin(theta)./abs(sin(theta)));
            sign(1)=1;
            
            joukowski_y = 2*this.b*this.e*(1-innerCircle.x/2/this.b)...
                .*sign.*sqrt(1-(innerCircle.x/2/this.b).^2)...
                + 2*this.b*this.beta*(1-(innerCircle.x/2/this.b).^2);
            airfoil = struct(...
                'x', innerCircle.x,...
                'y', joukowski_y);           
            
        end
        
        function [xGrid, yGrid, rGrid] = generatePhysicalGrid(this, outerCircle, airfoil)
            xGrid = zeros(this.inputs.i_max, this.inputs.j_max);
            yGrid = zeros(this.inputs.i_max, this.inputs.j_max);
            
            for i=1:this.inputs.i_max
                xGrid(i,:)=linspace(airfoil.x(i),outerCircle.x(i),this.inputs.j_max);
                yGrid(i,:)=linspace(airfoil.y(i),outerCircle.y(i),this.inputs.j_max);
            end
            
            rGrid = sqrt(xGrid.^2+yGrid.^2);
            colormap([0.95 0.95 0.95]);
            g = pcolor(xGrid, yGrid, rGrid);
            set(g, 'EdgeColor', [0.7 0.7 0.7]);
        end
        
         
        function [Eta1, Eta2] = generateComputationalGrid(this, eta1_limit, eta2_limit)
            eta1 = linspace(eta1_limit(1), eta1_limit(2), this.inputs.i_max);
            eta2 = linspace(eta2_limit(1), eta2_limit(2), this.inputs.j_max);
            [Eta2, Eta1] = meshgrid(eta2, eta1);
            plot(Eta1, Eta2, Eta1', Eta2', 'Color', 'b');
            axis equal;
            this.delta_eta1 = (eta1_limit(2)-eta1_limit(1))/(this.inputs.i_max-1);
            this.delta_eta2 = (eta2_limit(2)-eta2_limit(1))/(this.inputs.j_max-1);
        end
        
        function [] = transformationMetrics(this, xGrid, yGrid, eta1Grid, eta2Grid)
            this.dx_deta1 = this.zerosImaxJmax();
            this.dy_deta1 = this.zerosImaxJmax();
            for i = 1:this.inputs.i_max 
                if(i == 1)
                    this.dx_deta1(i,:) =(xGrid(i+1,:)-xGrid(i,:))./(eta1Grid(i+1,:)-eta1Grid(i,:));
                    this.dy_deta1(i,:) =(yGrid(i+1,:)-yGrid(i,:))./(eta1Grid(i+1,:)-eta1Grid(i,:));
                elseif (i == this.inputs.i_max)       
                    this.dx_deta1(i,:) =(xGrid(i,:)-xGrid(i-1,:))./(eta1Grid(i,:)-eta1Grid(i-1,:));
                    this.dy_deta1(i,:) =(yGrid(i,:)-yGrid(i-1,:))./(eta1Grid(i,:)-eta1Grid(i-1,:));
                else
                    this.dx_deta1(i,:) =(xGrid(i+1,:)-xGrid(i-1,:))./(eta1Grid(i+1,:)-eta1Grid(i-1,:));
                    this.dy_deta1(i,:) =(yGrid(i+1,:)-yGrid(i-1,:))./(eta1Grid(i+1,:)-eta1Grid(i-1,:));
                end
            end
            
            this.dx_deta2 = this.zerosImaxJmax();
            this.dy_deta2 = this.zerosImaxJmax();
            for j = 1:this.inputs.j_max 
                if(j == 1)
                    this.dx_deta2(:,j) =(xGrid(:,j+1)-xGrid(:,j))./(eta2Grid(:,j+1)-eta2Grid(:,j));
                    this.dy_deta2(:,j) =(yGrid(:,j+1)-yGrid(:,j))./(eta2Grid(:,j+1)-eta2Grid(:,j));
                elseif (j == this.inputs.j_max )       
                    this.dx_deta2(:,j) =(xGrid(:,j)-xGrid(:,j-1))./(eta2Grid(:,j)-eta2Grid(:,j-1));
                    this.dy_deta2(:,j) =(yGrid(:,j)-yGrid(:,j-1))./(eta2Grid(:,j)-eta2Grid(:,j-1));
                else
                    this.dx_deta2(:,j) =(xGrid(:,j+1)-xGrid(:,j-1))./(eta2Grid(:,j+1)-eta2Grid(:,j-1));
                    this.dy_deta2(:,j) =(yGrid(:,j+1)-yGrid(:,j-1))./(eta2Grid(:,j+1)-eta2Grid(:,j-1));
                end
            end
                     
            this.Jacobian = this.dx_deta1.*this.dy_deta2-this.dx_deta2.*this.dy_deta1;
            this.deta1_dx = this.dy_deta2./this.Jacobian;
            this.deta1_dy = -this.dx_deta2./this.Jacobian;
            this.deta2_dx = -this.dy_deta1./this.Jacobian;
            this.deta2_dy = this.dx_deta1./this.Jacobian;
            
            
            this.C11 = (this.dx_deta2.^2+this.dy_deta2.^2)./this.Jacobian;
            this.C22 = (this.dx_deta1.^2+this.dy_deta1.^2)./this.Jacobian;
            this.C12 = -(this.dx_deta1.*this.dx_deta2+this.dy_deta1.*this.dy_deta2)./this.Jacobian;
            
            this.C11_plusHalf_i = this.zerosImaxJmax();
            this.C11_negHalf_i = this.zerosImaxJmax();
            this.C22_plusHalf_i = this.zerosImaxJmax();
            this.C22_negHalf_i = this.zerosImaxJmax();
            this.C12_plusHalf_i = this.zerosImaxJmax();
            this.C12_negHalf_i = this.zerosImaxJmax();
            
            for i=1:this.inputs.i_max
               if(i ~= this.inputs.i_max)
                   this.C11_plusHalf_i(i, :) = (this.C11(i, :)+this.C11(i+1, :))/2; 
                   this.C22_plusHalf_i(i, :) = (this.C22(i, :)+this.C22(i+1, :))/2; 
                   this.C12_plusHalf_i(i, :) = (this.C12(i, :)+this.C12(i+1, :))/2; 
               end
               if(i ~= 1)
                   this.C11_negHalf_i(i, :) = (this.C11(i, :)+this.C11(i-1, :))/2; 
                   this.C22_negHalf_i(i, :) = (this.C22(i, :)+this.C22(i-1, :))/2; 
                   this.C12_negHalf_i(i, :) = (this.C12(i, :)+this.C12(i-1, :))/2; 
               end
            end
            this.C11_plusHalf_i(this.inputs.i_max, :) = this.C11_plusHalf_i(1, :);
            this.C22_plusHalf_i(this.inputs.i_max, :) = this.C22_plusHalf_i(1, :);
            this.C12_plusHalf_i(this.inputs.i_max, :) = this.C12_plusHalf_i(1, :);
            
            this.C11_negHalf_i(1, :) = this.C11_negHalf_i(this.inputs.i_max, :);
            this.C22_negHalf_i(1, :) = this.C22_negHalf_i(this.inputs.i_max, :);
            this.C12_negHalf_i(1, :) = this.C12_negHalf_i(this.inputs.i_max, :);
            
            this.C11_plusHalf_j = this.zerosImaxJmax();
            this.C11_negHalf_j = this.zerosImaxJmax();
            this.C22_plusHalf_j = this.zerosImaxJmax();
            this.C22_negHalf_j = this.zerosImaxJmax();
            this.C12_plusHalf_j = this.zerosImaxJmax();
            this.C12_negHalf_j = this.zerosImaxJmax();
            
            for j=2:this.inputs.j_max-1 
                this.C11_plusHalf_j(:,j) = (this.C11(:,j)+this.C11(:,j+1))/2;
                this.C11_negHalf_j(:,j) = (this.C11(:,j)+this.C11(:,j-1))/2;
                this.C22_plusHalf_j(:,j) = (this.C22(:,j)+this.C22(:,j+1))/2;
                this.C22_negHalf_j(:,j) = (this.C22(:,j)+this.C22(:,j-1))/2;
                this.C12_plusHalf_j(:,j) = (this.C12(:,j)+this.C12(:,j+1))/2;
                this.C12_negHalf_j(:,j) = (this.C12(:,j)+this.C12(:,j-1))/2;
            end 
            
            deta1_deta2 = this.delta_eta1/this.delta_eta2;
            
            this.s_i_j = this.C11_plusHalf_i + this.C11_negHalf_i + (this.C22_plusHalf_j + this.C22_negHalf_j)*(deta1_deta2)^2;

            this.s_in1_j = this.C11_negHalf_i - (this.C12_plusHalf_j - this.C12_negHalf_j)*(deta1_deta2/2)^2;
            this.s_ip1_j = this.C11_plusHalf_i + (this.C12_plusHalf_j - this.C12_negHalf_j)*(deta1_deta2/2)^2;

            this.s_i_jn1 = (this.C22_negHalf_j)*(deta1_deta2)^2 - (this.C12_plusHalf_i - this.C12_negHalf_i)*(deta1_deta2/4);
            this.s_i_jp1 = (this.C22_plusHalf_j)*(deta1_deta2)^2 + (this.C12_plusHalf_i - this.C12_negHalf_i)*(deta1_deta2/4);
            
            this.s_in1_jn1 = (this.C12_negHalf_i + this.C12_negHalf_j)*(deta1_deta2/4);
            this.s_ip1_jp1 = (this.C12_plusHalf_i + this.C12_plusHalf_j)*(deta1_deta2/4);
            
            this.s_in1_jp1 = -(this.C12_negHalf_i + this.C12_plusHalf_j)*(deta1_deta2/4);
            this.s_ip1_jn1 = -(this.C12_plusHalf_i + this.C12_negHalf_j)*(deta1_deta2/4);
            
            
        end
        
        function [psi] = calculateDirichletBoundary(this)
            theta = linspace(0, 2*pi, this.inputs.i_max);
            u = this.inputs.Vinf * this.cosa;
            v = this.inputs.Vinf * this.sina;
            
            x = this.inputs.R * cos(theta);
            y = this.inputs.R * sin(theta);
            
            psiBoundary = zeros(this.inputs.i_max, 1);
            
            for i=2:this.inputs.i_max
                psiBoundary(i,1) = psiBoundary(i-1, 1) - v*(x(i) - x(i-1)) + u*(y(i) - y(i-1));
            end
            
            psi = this.zerosImaxJmax();
            psi(:, end) = psiBoundary;
            
            for i=1:this.inputs.i_max
                psi(i, 1:end) = linspace(psi(i,1),psi(i,end),this.inputs.j_max);
            end
        end
        
        function [psi] = iterate(this, psi_pre)
            
            psi = this.zerosImaxJmax();
            psi(:, 1) = psi_pre(:, 1);
            psi(:,end) = psi_pre(:, end);
            
            for i=1:this.inputs.i_max
                for j=2:this.inputs.j_max-1
                    if i==1
                        psi(i,j)=(this.s_in1_j(i,j) * psi_pre(this.inputs.i_max-1,j) + this.s_ip1_j(i,j) * psi_pre(i+1,j)...
                            + this.s_i_jp1(i,j) * psi_pre(i,j+1) +this.s_i_jn1(i,j)*psi(i,j-1)...
                            + this.s_in1_jn1(i,j)*psi_pre(this.inputs.i_max-1,j-1)+this.s_in1_jp1(i,j)*psi_pre(this.inputs.i_max-1,j+1)+this.s_ip1_jn1(i,j)*psi_pre(i+1,j-1)+this.s_ip1_jp1(i,j)*psi_pre(i+1,j+1))...
                            /this.s_i_j(i,j);
                    elseif i== this.inputs.i_max
                        psi(i,j)=(this.s_in1_j(i,j) * psi(i-1,j) + this.s_ip1_j(i,j) * psi_pre(1+1,j)...
                            + this.s_i_jp1(i,j) * psi_pre(i,j+1) +this.s_i_jn1(i,j)*psi(i,j-1)...
                            + this.s_in1_jn1(i,j)*psi(i-1,j-1)+this.s_in1_jp1(i,j)*psi(i-1,j+1)+this.s_ip1_jn1(i,j)*psi(1+1,j-1)+this.s_ip1_jp1(i,j)*psi_pre(1+1,j+1))...
                            /this.s_i_j(i,j);
                    else
                        psi(i,j)=(this.s_in1_j(i,j) * psi(i-1,j) + this.s_ip1_j(i,j) * psi_pre(i+1,j)...
                            + this.s_i_jp1(i,j) * psi_pre(i,j+1) +this.s_i_jn1(i,j)*psi(i,j-1)...
                            + this.s_in1_jn1(i,j)*psi(i-1,j-1)+this.s_in1_jp1(i,j)*psi(i-1,j+1)+this.s_ip1_jn1(i,j)*psi_pre(i+1,j-1)+this.s_ip1_jp1(i,j)*psi_pre(i+1,j+1))...
                            /this.s_i_j(i,j);
                    end
                end
            end
        end
        
        function [m] = zerosImaxJmax(this)
            m = zeros(this.inputs.i_max, this.inputs.j_max);
        end

    end
end

