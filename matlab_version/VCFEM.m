classdef VCFEM < handle
    properties
        E_m
        E_c
        pr_m
        pr_c
        G_mm
        G_mc
        G_cc
        H_m
        H_c
        temporary_ke11
        temporary_ke12
        temporary_ke21
        temporary_ke22
        %phi1
        %phi2
        K
        F
        phi_m
        phi_c
        G
        d_m
        total_area_m
        total_area_c
    end
    
    methods
        function obj = VCFEM(E_m, E_c, pr_m, pr_c)
            obj.E_m = E_m;
            obj.E_c = E_c;
            obj.pr_m = pr_m;
            obj.pr_c = pr_c;
        end

        function G = matrix_G_on_single_edge(obj,n1, n2, l, x1, y1, x2, y2,G_type)
            if G_type==7
                G = [0.5*n1, 0, 0.5*n1, 0;
                    1/3*n1*y1+1/6*n1*y2, 0, 1/6*n1*y1+1/3*n1*y2, 0;
                    0, 0.5*n2, 0, 0.5*n2;
                    0, 1/3*n2*x1+1/6*n2*x2, 0, 1/6*n2*x1+1/3*n2*x2;
                    0.5*n2, 0.5*n1, 0.5*n2, 0.5*n1;
                    1/3*(n1*x1-n2*y1)+1/6*(n1*x2-n2*y2), -1/3*n1*y1-1/6*n1*y2, 1/6*(n1*x1-n2*y1)+1/3*(n1*x2-n2*y2), -1/6*n1*y1-1/3*n1*y2;
                    -1/3*n2*x1-1/6*n2*x2, 1/3*(n2*y1-n1*x1)+1/6*(n2*y2-n1*x2), -1/6*n2*x1-1/3*n2*x2, 1/6*(n2*y1-n1*x1)+1/3*(n2*y2-n1*x2)] * l;
            elseif G_type==12
                    G = [l*n1/2, 0, l*n1/2, 0;
                        l*n1*(2*y1+y2)/6, 0, l*n1*(y1+2*y2)/6, 0;
                        0, l*n2/2, 0, l*n2/2;
                        0, l*n2*(2*x1+x2)/6, 0, l*n2*(x1+2*x2)/6;
                        l*n2/2, l*n1/2, l*n2/2, l*n1/2;
                        l*(2*n1*x1+n1*x2-2*n2*y1-n2*y2)/6, -l*n1*(2*y1+y2)/6, l*(n1*x1+2*n1*x2-n2*y1-2*n2*y2)/6, -l*n1*(y1+2*y2)/6;
                        -l*n2*(2*x1+x2)/6, l*(-2*n1*x1-n1*x2+2*n2*y1+n2*y2)/6, -l*n2*(x1+2*x2)/6, l*(-n1*x1-2*n1*x2+n2*y1+2*n2*y2)/6;
                        l*n1*(3*y1^2+2*y1*y2+y2^2)/12, 0, l*n1*(y1^2+2*y1*y2+3*y2^2)/12, 0;
                        0, l*n2*(3*x1^2+2*x1*x2+x2^2)/12, 0, l*n2*(x1^2+2*x1*x2+3*x2^2)/12;
                        l*(3*n1*x1^2+2*n1*x1*x2+n1*x2^2-6*n2*x1*y1-2*n2*x1*y2-2*n2*x2*y1-2*n2*x2*y2)/12, l*(-6*n1*x1*y1-2*n1*x1*y2-2*n1*x2*y1-2*n1*x2*y2+3*n2*y1^2+2*n2*y1*y2+n2*y2^2)/12, l*(n1*x1^2+2*n1*x1*x2+3*n1*x2^2-2*n2*x1*y1-2*n2*x1*y2-2*n2*x2*y1-6*n2*x2*y2)/12, l*(-2*n1*x1*y1-2*n1*x1*y2-2*n1*x2*y1-6*n1*x2*y2+n2*y1^2+2*n2*y1*y2+3*n2*y2^2)/12;
                        l*(6*n1*x1*y1+2*n1*x1*y2+2*n1*x2*y1+2*n1*x2*y2-3*n2*y1^2-2*n2*y1*y2-n2*y2^2)/24, l*n1*(-3*y1^2-2*y1*y2-y2^2)/24, l*(2*n1*x1*y1+2*n1*x1*y2+2*n1*x2*y1+6*n1*x2*y2-n2*y1^2-2*n2*y1*y2-3*n2*y2^2)/24, l*n1*(-y1^2-2*y1*y2-3*y2^2)/24;
                        l*n2*(-3*x1^2-2*x1*x2-x2^2)/24, l*(-3*n1*x1^2-2*n1*x1*x2-n1*x2^2+6*n2*x1*y1+2*n2*x1*y2+2*n2*x2*y1+2*n2*x2*y2)/24, l*n2*(-x1^2-2*x1*x2-3*x2^2)/24, l*(-n1*x1^2-2*n1*x1*x2-3*n1*x2^2+2*n2*x1*y1+2*n2*x1*y2+2*n2*x2*y1+6*n2*x2*y2)/24;
                        ];
             elseif G_type==18
                 G = [l*n1/2, 0, l*n1/2, 0;
                        l*n1*(2*y1+y2)/6, 0, l*n1*(y1+2*y2)/6, 0;
                        0, l*n2/2, 0, l*n2/2;
                        0, l*n2*(2*x1+x2)/6, 0, l*n2*(x1+2*x2)/6;
                        l*n2/2, l*n1/2, l*n2/2, l*n1/2;
                        l*(2*n1*x1+n1*x2-2*n2*y1-n2*y2)/6, -l*n1*(2*y1+y2)/6, l*(n1*x1+2*n1*x2-n2*y1-2*n2*y2)/6, -l*n1*(y1+2*y2)/6;
                        -l*n2*(2*x1+x2)/6, l*(-2*n1*x1-n1*x2+2*n2*y1+n2*y2)/6, -l*n2*(x1+2*x2)/6, l*(-n1*x1-2*n1*x2+n2*y1+2*n2*y2)/6;
                        l*n1*(3*y1^2+2*y1*y2+y2^2)/12, 0, l*n1*(y1^2+2*y1*y2+3*y2^2)/12, 0;
                        0, l*n2*(3*x1^2+2*x1*x2+x2^2)/12, 0, l*n2*(x1^2+2*x1*x2+3*x2^2)/12;
                        l*(3*n1*x1^2+2*n1*x1*x2+n1*x2^2-6*n2*x1*y1-2*n2*x1*y2-2*n2*x2*y1-2*n2*x2*y2)/12, l*(-6*n1*x1*y1-2*n1*x1*y2-2*n1*x2*y1-2*n1*x2*y2+3*n2*y1^2+2*n2*y1*y2+n2*y2^2)/12, l*(n1*x1^2+2*n1*x1*x2+3*n1*x2^2-2*n2*x1*y1-2*n2*x1*y2-2*n2*x2*y1-6*n2*x2*y2)/12, l*(-2*n1*x1*y1-2*n1*x1*y2-2*n1*x2*y1-6*n1*x2*y2+n2*y1^2+2*n2*y1*y2+3*n2*y2^2)/12;
                        l*(6*n1*x1*y1+2*n1*x1*y2+2*n1*x2*y1+2*n1*x2*y2-3*n2*y1^2-2*n2*y1*y2-n2*y2^2)/24, l*n1*(-3*y1^2-2*y1*y2-y2^2)/24, l*(2*n1*x1*y1+2*n1*x1*y2+2*n1*x2*y1+6*n1*x2*y2-n2*y1^2-2*n2*y1*y2-3*n2*y2^2)/24, l*n1*(-y1^2-2*y1*y2-3*y2^2)/24;
                        l*n2*(-3*x1^2-2*x1*x2-x2^2)/24, l*(-3*n1*x1^2-2*n1*x1*x2-n1*x2^2+6*n2*x1*y1+2*n2*x1*y2+2*n2*x2*y1+2*n2*x2*y2)/24, l*n2*(-x1^2-2*x1*x2-3*x2^2)/24, l*(-n1*x1^2-2*n1*x1*x2-3*n1*x2^2+2*n2*x1*y1+2*n2*x1*y2+2*n2*x2*y1+6*n2*x2*y2)/24;
                        (l*(4*n1*x1^3 + 3*n1*x1^2*x2 + 2*n1*x1*x2^2 + n1*x2^3 - 12*n2*x1^2*y1 - 3*n2*x1^2*y2 - 6*n2*x1*x2*y1 - 4*n2*x1*x2*y2 - 2*n2*x2^2*y1 - 3*n2*x2^2*y2))/20, (l*(-12*n1*x1^2*y1 - 3*n1*x1^2*y2 - 6*n1*x1*x2*y1 - 4*n1*x1*x2*y2 - 2*n1*x2^2*y1 - 3*n1*x2^2*y2 + 12*n2*x1*y1^2 + 6*n2*x1*y1*y2 + 2*n2*x1*y2^2 + 3*n2*x2*y1^2 + 4*n2*x2*y1*y2 + 3*n2*x2*y2^2))/20, (l*(n1*x1^3 + 2*n1*x1^2*x2 + 3*n1*x1*x2^2 + 4*n1*x2^3 - 3*n2*x1^2*y1 - 2*n2*x1^2*y2 - 4*n2*x1*x2*y1 - 6*n2*x1*x2*y2 - 3*n2*x2^2*y1 - 12*n2*x2^2*y2))/20, (l*(-3*n1*x1^2*y1 - 2*n1*x1^2*y2 - 4*n1*x1*x2*y1 - 6*n1*x1*x2*y2 - 3*n1*x2^2*y1 - 12*n1*x2^2*y2 + 3*n2*x1*y1^2 + 4*n2*x1*y1*y2 + 3*n2*x1*y2^2 + 2*n2*x2*y1^2 + 6*n2*x2*y1*y2 + 12*n2*x2*y2^2))/20;
                        (l*n1*(4*y1^3 + 3*y1^2*y2 + 2*y1*y2^2 + y2^3))/20, 0, (l*n1*(y1^3 + 2*y1^2*y2 + 3*y1*y2^2 + 4*y2^3))/20, 0;
                        (l*(12*n1*x1^2*y1 + 3*n1*x1^2*y2 + 6*n1*x1*x2*y1 + 4*n1*x1*x2*y2 + 2*n1*x2^2*y1 + 3*n1*x2^2*y2 - 12*n2*x1*y1^2 - 6*n2*x1*y1*y2 - 2*n2*x1*y2^2 - 3*n2*x2*y1^2 - 4*n2*x2*y1*y2 - 3*n2*x2*y2^2))/20, (l*(-12*n1*x1*y1^2 - 6*n1*x1*y1*y2 - 2*n1*x1*y2^2 - 3*n1*x2*y1^2 - 4*n1*x2*y1*y2 - 3*n1*x2*y2^2 + 4*n2*y1^3 + 3*n2*y1^2*y2 + 2*n2*y1*y2^2 + n2*y2^3))/20, (l*(3*n1*x1^2*y1 + 2*n1*x1^2*y2 + 4*n1*x1*x2*y1 + 6*n1*x1*x2*y2 + 3*n1*x2^2*y1 + 12*n1*x2^2*y2 - 3*n2*x1*y1^2 - 4*n2*x1*y1*y2 - 3*n2*x1*y2^2 - 2*n2*x2*y1^2 - 6*n2*x2*y1*y2 - 12*n2*x2*y2^2))/20, (l*(-3*n1*x1*y1^2 - 4*n1*x1*y1*y2 - 3*n1*x1*y2^2 - 2*n1*x2*y1^2 - 6*n1*x2*y1*y2 - 12*n1*x2*y2^2 + n2*y1^3 + 2*n2*y1^2*y2 + 3*n2*y1*y2^2 + 4*n2*y2^3))/20;
                        (l*(12*n1*x1*y1^2 + 6*n1*x1*y1*y2 + 2*n1*x1*y2^2 + 3*n1*x2*y1^2 + 4*n1*x2*y1*y2 + 3*n1*x2*y2^2 - 4*n2*y1^3 - 3*n2*y1^2*y2 - 2*n2*y1*y2^2 - n2*y2^3))/60, (l*n1*(-4*y1^3 - 3*y1^2*y2 - 2*y1*y2^2 - y2^3))/60, (l*(3*n1*x1*y1^2 + 4*n1*x1*y1*y2 + 3*n1*x1*y2^2 + 2*n1*x2*y1^2 + 6*n1*x2*y1*y2 + 12*n1*x2*y2^2 - n2*y1^3 - 2*n2*y1^2*y2 - 3*n2*y1*y2^2 - 4*n2*y2^3))/60, (l*n1*(-y1^3 - 2*y1^2*y2 - 3*y1*y2^2 - 4*y2^3))/60;
                        (l*n2*(-4*x1^3 - 3*x1^2*x2 - 2*x1*x2^2 - x2^3))/60,(l*(-4*n1*x1^3 - 3*n1*x1^2*x2 - 2*n1*x1*x2^2 - n1*x2^3 + 12*n2*x1^2*y1 + 3*n2*x1^2*y2 + 6*n2*x1*x2*y1 + 4*n2*x1*x2*y2 + 2*n2*x2^2*y1 + 3*n2*x2^2*y2))/60, (l*n2*(-x1^3 - 2*x1^2*x2 - 3*x1*x2^2 - 4*x2^3))/60, (l*(-n1*x1^3 - 2*n1*x1^2*x2 - 3*n1*x1*x2^2 - 4*n1*x2^3 + 3*n2*x1^2*y1 + 2*n2*x1^2*y2 + 4*n2*x1*x2*y1 + 6*n2*x1*x2*y2 + 3*n2*x2^2*y1 + 12*n2*x2^2*y2))/60;
                        0, (l*n2*(4*x1^3 + 3*x1^2*x2 + 2*x1*x2^2 + x2^3))/20, 0, (l*n2*(x1^3 + 2*x1^2*x2 + 3*x1*x2^2 + 4*x2^3))/20];
             end
        end
        function [G_mm, G_mc, G_cc] = matrix_G(obj, element)
            n_V_E = element.edge_m_num * 2;
            n_V_C = element.edge_c_num * 2;
            G_type_m = element.n_beta_M;
            G_type_c = element.n_beta_C;
            n_beta_M = element.n_beta_M;
            n_beta_C = element.n_beta_C;
            G_mm = zeros(n_beta_M, n_V_E);
            G_mc = zeros(n_beta_M, n_V_C);
            G_cc = zeros(n_beta_C,n_V_C);
            for i = 1:element.edge_m_num
                n1 = element.edge_m{i}.n1;
                n2 = element.edge_m{i}.n2;
                l = element.edge_m{i}.length;
                x1 = element.node_m(2 * i - 1);
                y1 = element.node_m(2 * i);
                x2 = element.node_m(2 * mod(i, element.edge_m_num) + 1);
                y2 = element.node_m(2 * mod(i, element.edge_m_num) + 2);
                G = obj.matrix_G_on_single_edge(n1, n2, l, x1, y1, x2, y2,G_type_m);
                G_mm(:, 2 * i - 1 : 2 * i) = G_mm(:, 2 * i - 1 : 2 * i) + G(:, 1:2);
                G_mm(:, 2 * mod(i, element.edge_m_num) + 1 : 2 * mod(i, element.edge_m_num) + 2) = G_mm(:, 2 * mod(i, element.edge_m_num) + 1 : 2 * mod(i, element.edge_m_num) + 2) + G(:, 3:4);
            end
            % Assuming edge_c, node_c and edge_c_num are defined similarly to edge_m, node_m and edge_m_num
            for i = 1:element.edge_c_num
                n1 = element.edge_c{i}.n1;
                n2 = element.edge_c{i}.n2;
                l = element.edge_c{i}.length;
                x1 = element.node_c(2 * i - 1);
                y1 = element.node_c(2 * i);
                x2 = element.node_c(2 * mod(i, element.edge_c_num) + 1);
                y2 = element.node_c(2 * mod(i, element.edge_c_num) + 2);
                G = obj.matrix_G_on_single_edge(n1, n2, l, x1, y1, x2, y2,G_type_m);
                G_mc(:, 2 * i - 1 : 2 * i) = G_mc(:, 2 * i - 1 : 2 * i) + G(:, 1:2);
                G_mc(:, 2 * mod(i, element.edge_c_num) + 1 : 2 * mod(i, element.edge_c_num) + 2) = G_mc(:, 2 * mod(i, element.edge_c_num) + 1 : 2 * mod(i, element.edge_c_num) + 2) + G(:, 3:4);
            end
            for i = 1:element.edge_c_num
                n1 = element.edge_c{i}.n1;
                n2 = element.edge_c{i}.n2;
                l = element.edge_c{i}.length;
                x1 = element.node_c(2 * i - 1);
                y1 = element.node_c(2 * i);
                x2 = element.node_c(2 * mod(i, element.edge_c_num) + 1);
                y2 = element.node_c(2 * mod(i, element.edge_c_num) + 2);
                G = obj.matrix_G_on_single_edge(n1, n2, l, x1, y1, x2, y2,G_type_c);
                G_cc(:, 2 * i - 1 : 2 * i) = G_cc(:, 2 * i - 1 : 2 * i) + G(:, 1:2);
                G_cc(:, 2 * mod(i, element.edge_c_num) + 1 : 2 * mod(i, element.edge_c_num) + 2) = G_cc(:, 2 * mod(i, element.edge_c_num) + 1 : 2 * mod(i, element.edge_c_num) + 2) + G(:, 3:4);
            end
            element.G_mm = G_mm;
            element.G_mc = G_mc;
            element.G_cc = G_cc;
        end
        function PTSP = get_PTSP(obj, x, y, type,PTSP_type)
            if strcmp(type, 'm')
                E = obj.E_m;
                pr = obj.pr_m;
                v = pr;
            elseif strcmp(type, 'c')
                E = obj.E_c;
                pr = obj.pr_c;
                v = pr;
            end
            if PTSP_type ==7
            PTSP = (1/E) * [1, y, -pr, -pr*x, 0, x, -pr*y;
                            y, y^2, -pr*y, -pr*x*y, 0, x*y, -pr*y^2;
                            -pr, -pr*y, 1, x, 0, -pr*x, y;
                            -pr*x, -pr*x*y, x, x^2, 0, -pr*x^2, x*y;
                            0, 0, 0, 0, 2*(1+pr), -2*(1+pr)*y, -2*(1+pr)*x;
                            x, x*y, -pr*x, -pr*x^2, -2*(1+pr)*y, x^2+2*(1+pr)*y^2, (2+pr)*x*y;
                            -pr*y, -pr*y^2, y, x*y, -2*(1+pr)*x, (2+pr)*x*y, y^2+2*(1+pr)*x^2];
            elseif PTSP_type==12
            PTSP = (1/E) * [1, y, -v, -v*x, 0, x, -v*y, y^2, -v*x^2, -v*y^2+x^2, x*y, -v*x*y;
                            y, y^2, -v*y, -v*x*y, 0, x*y, -v*y^2, y^3, -v*x^2*y, -v*y^3+x^2*y, x*y^2, -v*x*y^2;
                            -v, -v*y, 1, x, 0, -v*x, y, -v*y^2, x^2, -v*x^2+y^2, -v*x*y, x*y;
                            -v*x, -v*x*y, x, x^2, 0, -v*x^2, x*y, -v*x*y^2, x^3, -v*x^3+x*y^2, -v*x^2*y, x^2*y;
                            0, 0, 0, 0, 2*v+2, -y*(2*v+2), -x*(2*v+2), 0, 0, -2*x*y*(2*v+2), -(y^2*(2*v+2))/2, -(x^2*(2*v+2))/2;
                            x, x*y, -v*x, -v*x^2, -y*(2*v+2), x^2+y^2*(2*v+2), -v*x*y+x*y*(2*v+2), x*y^2, -v*x^3, -v*x*y^2+x^3+2*x*y^2*(2*v+2), x^2*y+(y^3*(2*v+2))/2, -v*x^2*y+(x^2*y*(2*v+2))/2;
                            -v*y, -v*y^2, y, x*y, -x*(2*v+2), -v*x*y+x*y*(2*v+2), x^2*(2*v+2)+y^2, -v*y^3, x^2*y, -v*x^2*y+2*x^2*y*(2*v+2)+y^3, -v*x*y^2+(x*y^2*(2*v+2))/2, (x^3*(2*v+2))/2+x*y^2;
                            y^2, y^3, -v*y^2, -v*x*y^2, 0, x*y^2, -v*y^3, y^4, -v*x^2*y^2, -v*y^4+x^2*y^2, x*y^3, -v*x*y^3;
                            -v*x^2, -v*x^2*y, x^2, x^3, 0, -v*x^3, x^2*y, -v*x^2*y^2, x^4, -v*x^4+x^2*y^2, -v*x^3*y, x^3*y;
                            -v*y^2+x^2, y*(-v*y^2+x^2), -v*x^2+y^2, x*(-v*x^2+y^2), -2*x*y*(2*v+2), 2*x*y^2*(2*v+2)+x*(-v*y^2+x^2), 2*x^2*y*(2*v+2)+y*(-v*x^2+y^2), y^2*(-v*y^2+x^2), x^2*(-v*x^2+y^2), 4*x^2*y^2*(2*v+2)+x^2*(-v*y^2+x^2)+y^2*(-v*x^2+y^2), x*y^3*(2*v+2)+x*y*(-v*y^2+x^2), x^3*y*(2*v+2)+x*y*(-v*x^2+y^2);
                            x*y, x*y^2, -v*x*y, -v*x^2*y, -(y^2*(2*v+2))/2, x^2*y+(y^3*(2*v+2))/2, -v*x*y^2+(x*y^2*(2*v+2))/2, x*y^3, -v*x^3*y, -v*x*y^3+x^3*y+x*y^3*(2*v+2), x^2*y^2+(y^4*(2*v+2))/4, -v*x^2*y^2+(x^2*y^2*(2*v+2))/4;
                            -v*x*y, -v*x*y^2, x*y, x^2*y, -(x^2*(2*v+2))/2, -v*x^2*y+(x^2*y*(2*v+2))/2, (x^3*(2*v+2))/2+x*y^2, -v*x*y^3, x^3*y, -v*x^3*y+x^3*y*(2*v+2)+x*y^3, -v*x^2*y^2+(x^2*y^2*(2*v+2))/4, (x^4*(2*v+2))/4+x^2*y^2
                            ];
            elseif PTSP_type==18
            PTSP = (1/E) *  [1, y, -v, -v*x, 0, x, -v*y, y^2, -v*x^2, -v*y^2 + x^2, x*y, -v*x*y, x*(-3*v*y^2 + x^2), y^3, y*(-v*y^2 + 3*x^2), x*y^2, -v*x^2*y, -v*x^3;
                             y, y^2, -v*y, -v*x*y, 0, x*y, -v*y^2, y^3, -v*x^2*y, y*(-v*y^2 + x^2), x*y^2, -v*x*y^2, x*y*(-3*v*y^2 + x^2), y^4, y^2*(-v*y^2 + 3*x^2), x*y^3, -v*x^2*y^2, -v*x^3*y;
                             -v, -v*y, 1, x, 0, -v*x, y, -v*y^2, x^2, -v*x^2 + y^2, -v*x*y, x*y, x*(-v*x^2 + 3*y^2), -v*y^3, y*(-3*v*x^2 + y^2), -v*x*y^2, x^2*y, x^3;
                             -v*x, -v*x*y, x, x^2, 0, -v*x^2, x*y, -v*x*y^2, x^3, x*(-v*x^2 + y^2), -v*x^2*y, x^2*y, x^2*(-v*x^2 + 3*y^2), -v*x*y^3, x*y*(-3*v*x^2 + y^2), -v*x^2*y^2, x^3*y, x^4;
                             0, 0, 0, 0, 2*v + 2, -2*y*(v + 1), -2*x*(v + 1), 0, 0, -4*x*y*(v + 1), y^2*(-v - 1), x^2*(-v - 1), -6*x^2*y*(v + 1), 0, -6*x*y^2*(v + 1), (2*y^3*(-v - 1))/3, (2*x^3*(-v - 1))/3, 0;
                             x, x*y, -v*x, -v*x^2, -2*y*(v + 1), x^2 + 2*y^2*(v + 1), x*y*(v + 2), x*y^2, -v*x^3, x*(3*v*y^2 + x^2 + 4*y^2), y*(x^2 + y^2*(v + 1)), x^2*y, x^2*(3*v*y^2 + x^2 + 6*y^2), x*y^3, x*y*(5*v*y^2 + 3*x^2 + 6*y^2), y^2*(x^2 + (2*y^2*(v + 1))/3), (x^3*y*(2 - v))/3, -v*x^4;
                             -v*y, -v*y^2, y, x*y, -2*x*(v+1), x*y*(v+2), 2*x^2*(v+1) + y^2, -v*y^3, x^2*y, y*(3*v*x^2 + 4*x^2 + y^2), x*y^2, x*(x^2*(v+1) + y^2), x*y*(5*v*x^2 + 6*x^2 + 3*y^2), -v*y^4, y^2*(3*v*x^2 + 6*x^2 + y^2), (x*y^3*(2-v))/3, x^2*((2*x^2*(v+1))/3 + y^2), x^3*y;
                             y^2, y^3, -v*y^2, -v*x*y^2, 0, x*y^2, -v*y^3, y^4, -v*x^2*y^2, y^2*(-v*y^2 + x^2), x*y^3, -v*x*y^3, x*y^2*(-3*v*y^2 + x^2), y^5, y^3*(-v*y^2 + 3*x^2), x*y^4, -v*x^2*y^3, -v*x^3*y^2;
                             -v*x^2, -v*x^2*y, x^2, x^3, 0, -v*x^3, x^2*y, -v*x^2*y^2, x^4, x^2*(-v*x^2 + y^2), -v*x^3*y, x^3*y, x^3*(-v*x^2 + 3*y^2), -v*x^2*y^3, x^2*y*(-3*v*x^2 + y^2), -v*x^3*y^2, x^4*y, x^5;
                             -v*y^2 + x^2, y*(-v*y^2 + x^2), -v*x^2 + y^2, x*(-v*x^2 + y^2), -4*x*y*(v + 1), x*(3*v*y^2 + x^2 + 4*y^2), y*(3*v*x^2 + 4*x^2 + y^2), y^2*(-v*y^2 + x^2), x^2*(-v*x^2 + y^2), 6*v*x^2*y^2 + x^4 + 8*x^2*y^2 + y^4, x*y*(v*y^2 + x^2 + 2*y^2), x*y*(v*x^2 + 2*x^2 + y^2), x*(8*v*x^2*y^2 + x^4 + 12*x^2*y^2 + 3*y^4), y^3*(-v*y^2 + x^2), y*(8*v*x^2*y^2 + 3*x^4 + 12*x^2*y^2 + y^4), (x*y^2*(v*y^2 + 3*x^2 + 4*y^2))/3, (x^2*y*(v*x^2 + 4*x^2 + 3*y^2))/3, x^3*(-v*x^2 + y^2);
                             x*y, x*y^2, -v*x*y, -v*x^2*y, y^2*(-v-1), y*(x^2 + y^2*(v+1)), x*y^2, x*y^3, -v*x^3*y, x*y*(v*y^2 + x^2 + 2*y^2), y^2*(x^2 + y^2*(v+1))/2, x^2*y*(1 - v)/2, x^2*y*(x^2 + 3*y^2), x*y^4, x*y^2*(2*v*y^2 + 3*x^2 + 3*y^2), y^3*(x^2 + y^2*(v+1))/3, x^3*y^2*(1 - 2*v)/3, -v*x^4*y;
                             -v*x*y, -v*x*y^2, x*y, x^2*y, x^2*(-v-1), x^2*y, x*(x^2*(v+1) + y^2), -v*x*y^3, x^3*y, x*y*(v*x^2 + 2*x^2 + y^2), x^2*y^2*(1 - v)/2, x^2*((x^2*(v+1))/2 + y^2), x^2*y*(2*v*x^2 + 3*x^2 + 3*y^2), -v*x*y^4, x*y^2*(3*x^2 + y^2), x^2*y^3*(1 - 2*v)/3, x^3*((x^2*(v+1))/3 + y^2), x^4*y;
                             x*(-3*v*y^2 + x^2), x*y*(-3*v*y^2 + x^2), x*(-v*x^2 + 3*y^2), x^2*(-v*x^2 + 3*y^2), -6*x^2*y*(v+1), x^2*(3*v*y^2 + x^2 + 6*y^2), x*y*(5*v*x^2 + 6*x^2 + 3*y^2), x*y^2*(-3*v*y^2 + x^2), x^3*(-v*x^2 + 3*y^2), x*(8*v*x^2*y^2 + x^4 + 12*x^2*y^2 + 3*y^4), x^2*y*(x^2 + 3*y^2), x^2*y*(2*v*x^2 + 3*x^2 + 3*y^2), x^2*(12*v*x^2*y^2 + x^4 + 18*x^2*y^2 + 9*y^4), x*y^3*(-3*v*y^2 + x^2), x*y*(8*v*x^2*y^2 + 3*x^4 + 18*x^2*y^2 + 3*y^4), x^2*y^2*(-v*y^2 + x^2 + 2*y^2), x^3*y*(v*x^2 + 2*x^2 + 3*y^2), x^4*(-v*x^2 + 3*y^2);
                             y^3, y^4, -v*y^3, -v*x*y^3, 0, x*y^3, -v*y^4, y^5, -v*x^2*y^3, y^3*(-v*y^2 + x^2), x*y^4, -v*x*y^4, x*y^3*(-3*v*y^2 + x^2), y^6, y^4*(-v*y^2 + 3*x^2), x*y^5, -v*x^2*y^4, -v*x^3*y^3;
                             y*(-v*y^2 + 3*x^2), y^2*(-v*y^2 + 3*x^2), y*(-3*v*x^2 + y^2), x*y*(-3*v*x^2 + y^2), -6*x*y^2*(v+1), x*y*(5*v*y^2 + 3*x^2 + 6*y^2), y^2*(3*v*x^2 + 6*x^2 + y^2), y^3*(-v*y^2 + 3*x^2), x^2*y*(-3*v*x^2 + y^2), y*(8*v*x^2*y^2 + 3*x^4 + 12*x^2*y^2 + y^4), x*y^2*(2*v*y^2 + 3*x^2 + 3*y^2), x*y^2*(3*x^2 + y^2), x*y*(8*v*x^2*y^2 + 3*x^4 + 18*x^2*y^2 + 3*y^4), y^4*(-v*y^2 + 3*x^2), y^2*(12*v*x^2*y^2 + 9*x^4 + 18*x^2*y^2 + y^4), x*y^3*(v*y^2 + 3*x^2 + 2*y^2), x^2*y^2*(-v*x^2 + 2*x^2 + y^2), x^3*y*(-3*v*x^2 + y^2);
                             x*y^2, x*y^3, -v*x*y^2, -v*x^2*y^2, (2*y^3*(-v-1))/3, y^2*(x^2 + (2*y^2*(v+1))/3), (x*y^3*(2-v))/3, x*y^4, -v*x^3*y^2, (x*y^2*(v*y^2 + 3*x^2 + 4*y^2))/3, y^3*(x^2 + (y^2*(v+1))/3), (x^2*y^3*(1-2*v))/3, x^2*y^2*(-v*y^2 + x^2 + 2*y^2), x*y^5, x*y^3*(v*y^2 + 3*x^2 + 2*y^2), y^4*(x^2 + (2*y^2*(v+1))/9), (x^3*y^3*(2-7*v))/9, -v*x^4*y^2;
                             -v*x^2*y, -v*x^2*y^2, x^2*y, x^3*y, (2*x^3*(-v-1))/3, (x^3*y*(2-v))/3, x^2*((2*x^2*(v+1))/3 + y^2), -v*x^2*y^3, x^4*y, (x^2*y*(v*x^2 + 4*x^2 + 3*y^2))/3, (x^3*y^2*(1-2*v))/3, x^3*((x^2*(v+1))/3 + y^2), x^3*y*(v*x^2 + 2*x^2 + 3*y^2), -v*x^2*y^4, x^2*y^2*(-v*x^2 + 2*x^2 + y^2), (x^3*y^3*(2-7*v))/9, x^4*((2*x^2*(v+1))/9 + y^2), x^5*y;
                             -v*x^3, -v*x^3*y, x^3, x^4, 0, -v*x^4, x^3*y, -v*x^3*y^2, x^5, x^3*(-v*x^2 + y^2), -v*x^4*y, x^4*y, x^4*(-v*x^2 + 3*y^2), -v*x^3*y^3, x^3*y*(-3*v*x^2 + y^2), -v*x^4*y^2, x^5*y, x^6];
            end
        end
        function J = Jacobian(obj, x1, y1, x2, y2, x3, y3)
            J = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        end

        function x = transform(obj, x1, x2, x3, r, s)
            x = (1 - r - s) * x1 + r * x2 + s * x3;
        end
        function integral = tri_integral(obj,x1,y1,x2,y2,x3,y3,J,type,PTSP_type)%type为相的类型
            tri_integral = zeros(PTSP_type,PTSP_type);
            if PTSP_type ==7
                F1 = obj.get_PTSP(obj.transform(x1, x2, x3, 1/3, 1/3), obj.transform(y1, y2, y3, 1/3, 1/3), type,PTSP_type);
                F2 = obj.get_PTSP(obj.transform(x1, x2, x3, 3/5, 1/5), obj.transform(y1, y2, y3, 3/5, 1/5), type,PTSP_type);
                F3 = obj.get_PTSP(obj.transform(x1, x2, x3, 1/5, 1/5), obj.transform(y1, y2, y3, 1/5, 1/5), type,PTSP_type);
                F4 = obj.get_PTSP(obj.transform(x1, x2, x3, 1/5, 3/5), obj.transform(y1, y2, y3, 1/5, 3/5), type,PTSP_type);
                tri_integral = 0.5 * J * (-27/48 * F1 + 25/48 * F2 + 25/48 * F3 + 25/48 * F4);  % 3阶代数精度
            elseif PTSP_type ==12
                r = [0.44594849, 0.09157621, 0.10810302, 0.44594849, 0.81684757,0.09157621];
                s = [0.44594849, 0.09157621, 0.44594849, 0.10810302, 0.09157621,0.81684757];
                half_weights = [0.11169079, 0.05497587, 0.11169079, 0.11169079, 0.05497587,0.05497587];
                for i=1:6
                    tri_integral = tri_integral + half_weights(i)*J*obj.get_PTSP(obj.transform(x1, x2, x3, r(i), s(i)), obj.transform(y1, y2, y3, r(i), s(i)), type,PTSP_type);
                end
            elseif PTSP_type==18
                r = [0.33333333, 0.26034597, 0.0651301 , 0.3128655 , 0.47930807,...
                    0.26034597, 0.86973979, 0.0651301 , 0.63844419, 0.63844419,...
                    0.04869032, 0.04869032, 0.3128655 ];
                s = [0.33333333, 0.26034597, 0.0651301 , 0.63844419, 0.26034597,...
                    0.47930807, 0.0651301 , 0.86973979, 0.3128655 , 0.04869032,...
                     0.3128655 , 0.63844419, 0.04869032];
                half_weights = [-0.07478502,  0.08780763,  0.02667362,  0.03855688,  0.08780763,...
                                 0.08780763,  0.02667362,  0.02667362,  0.03855688,  0.03855688,...
                                 0.03855688,  0.03855688,  0.03855688];
                for i=1:13
                    tri_integral = tri_integral + half_weights(i)*J*obj.get_PTSP(obj.transform(x1, x2, x3, r(i), s(i)), obj.transform(y1, y2, y3, r(i), s(i)), type,PTSP_type);
                end              
            end 
            integral = tri_integral;
        end
        function [H_m,H_c] = matrix_H(obj, element)
            if element.edge_m_num <= 5
                PTSP_type_m = 7;
                H_m = zeros(7,7);
            elseif element.edge_m_num <= 7
                PTSP_type_m = 12;
                H_m = zeros(12,12);
            else
                PTSP_type_m = 18;
                H_m = zeros(18,18);
            end

            if element.edge_c_num <= 5
                PTSP_type_c = 7;
                H_c = zeros(7,7);
            elseif element.edge_c_num <= 7
                PTSP_type_c = 12;
                H_c = zeros(12,12);
            else
                PTSP_type_c = 18;
                H_c = zeros(18,18);
            end
            %PTSP_type_c =PTSP_type_m;
            x1 = element.node_m(1);
            y1 = element.node_m(2);
            for i = 2:element.edge_m_num-1
                x2 = element.node_m(2 * (i - 1) + 1);
                y2 = element.node_m(2 * (i - 1) + 2);
                x3 = element.node_m(2 * i + 1);
                y3 = element.node_m(2 * i + 2);
                J = obj.Jacobian(x1, y1, x2, y2, x3, y3);
                if J == 0
                    continue;
                end
                triangle_integral = obj.tri_integral(x1,y1,x2,y2,x3,y3,J,'m',PTSP_type_m);
                H_m = H_m + triangle_integral;
            end
            x1 = element.node_c(1);
            y1 = element.node_c(2);
            for i = 2:element.edge_c_num-1
                x2 = element.node_c(2 * (i - 1) + 1);
                y2 = element.node_c(2 * (i - 1) + 2);
                x3 = element.node_c(2 * i + 1);
                y3 = element.node_c(2 * i + 2);
                J = obj.Jacobian(x1, y1, x2, y2, x3, y3);
                if J == 0
                    continue;
                end
                triangle_integral = obj.tri_integral(x1,y1,x2,y2,x3,y3,J,'m',PTSP_type_m);
                H_m = H_m - triangle_integral;
            end
            for i = 2:element.edge_c_num-1
                x2 = element.node_c(2 * (i - 1) + 1);
                y2 = element.node_c(2 * (i - 1) + 2);
                x3 = element.node_c(2 * i + 1);
                y3 = element.node_c(2 * i + 2);
                J = obj.Jacobian(x1, y1, x2, y2, x3, y3);
                if J == 0
                    continue;
                end
                triangle_integral = obj.tri_integral(x1,y1,x2,y2,x3,y3,J,'c',PTSP_type_c);
                H_c = H_c + triangle_integral;
            element.H_m = H_m;
            element.H_c = H_c;
            end
        end
        function [temporary_ke11, temporary_ke12, temporary_ke21, temporary_ke22] = temporary_keij(obj, H_m, H_c, G_mm, G_cc, G_mc)
            temporary_ke11 = G_mm' / H_m * G_mm;
            temporary_ke12 = - G_mm' / H_m * G_mc;
            temporary_ke21 = temporary_ke12';
            temporary_ke22 = G_cc' * (inv(H_c)) * G_cc + G_mc' * (inv(H_m)) * G_mc;
            obj.temporary_ke11 = temporary_ke11;
            obj.temporary_ke12 = temporary_ke12;
            obj.temporary_ke21 = temporary_ke21;
            obj.temporary_ke22 = temporary_ke22;
        end
        function [phi1, phi2] = phi(obj, element)
            phi_m = zeros(element.edge_m_num * 2, 3);
            phi_c = zeros(element.edge_c_num * 2, 3);
            col1 = double(mod((0:element.edge_m_num * 2 - 1), 2) == 0);
            col2 = double(mod((0:element.edge_m_num * 2 - 1), 2) == 1);
            x = element.node_m(1:2:end);
            y = element.node_m(2:2:end);
            col3 = zeros(size(col1));
            col3(1:2:end) = -y;
            col3(2:2:end) = x;
            phi_m(:, 1) = col1;
            phi_m(:, 2) = col2;
            phi_m(:, 3) = col3;
            obj.phi_m = phi_m;
            col1 = double(mod((0:element.edge_c_num * 2 - 1), 2) == 0);
            col2 = double(mod((0:element.edge_c_num * 2 - 1), 2) == 1);
            x = element.node_c(1:2:end);
            y = element.node_c(2:2:end);
            col3 = zeros(size(col1));
            col3(1:2:end) = -y;
            col3(2:2:end) = x;
            phi_c(:, 1) = col1;
            phi_c(:, 2) = col2;
            phi_c(:, 3) = col3;
            obj.phi_c = phi_c;
            phi1 = (phi_m/(phi_m'*phi_m))';%inv(phi_m' * phi_m) * phi_m';
            phi2 = -(phi_c/(phi_c'*phi_c))';%inv(phi_c' * phi_c) * phi_c';
            %obj.phi1 = phi1;
            %obj.phi2 = phi2;
        end
        function [Ke,Ke11,Ke12,Ke21,Ke22] = Keij(obj, phi1, phi2, temporary_ke11, temporary_ke12, temporary_ke21, temporary_ke22)
            Ke11 = temporary_ke11;
            Ke12 = [temporary_ke12, phi1'];
            %Ke21 = Ke12';
            Ke21 = [temporary_ke21; phi1];
            zero_matrix_temp = zeros(size(phi2, 1));
            Ke22 = [temporary_ke22, phi2'; phi2, zero_matrix_temp];
            Ke = Ke11 - Ke12 / Ke22 * Ke21;
        end
        function K = assembly_global_stiffness_matrix(obj, mesh)
            K = zeros(mesh.node_num * 2);
            for i = 1:mesh.element_num
                element = mesh.elements{i};
                [G_mm, G_mc, G_cc] = obj.matrix_G(element);
                [H_m, H_c] = obj.matrix_H(element);
                [temporary_ke11, temporary_ke12, temporary_ke21, temporary_ke22] = obj.temporary_keij(H_m, H_c, G_mm, G_cc, G_mc);
                [phi1, phi2] = obj.phi(element);
                [Ke,Ke11,Ke12,Ke21,Ke22] = obj.Keij(phi1, phi2, temporary_ke11, temporary_ke12, temporary_ke21, temporary_ke22);
                B = zeros(element.edge_m_num * 2, mesh.node_num * 2);
                for j = 1:element.edge_m_num
                    B(2 * j - 1, 2 * element.node_m_id(j) - 1) = 1;
                    B(2 * j, 2 * element.node_m_id(j)) = 1;
                end
                element.Ke = Ke;
                element.Ke11 = Ke11;
                element.Ke12 = Ke12;
                element.Ke21 = Ke21;
                element.Ke22 = Ke22;
                K = K + B' * Ke * B;
            end
            obj.K = K;
        end
        function F = calculate_global_nodal_load(obj, mesh, load)
            F = zeros(mesh.node_num * 2, 1);
            for k = 1:load.load_num
                i = load.node_i(k);
                j = load.node_j(k);
                B = zeros(4, mesh.node_num * 2);
                B(1:2, 2 * i - 1:2 * i) = eye(2);
                B(3:4, 2 * j - 1:2 * j) = eye(2);
                xi = mesh.nodes{i}.x;
                yi = mesh.nodes{i}.y;
                xj = mesh.nodes{j}.x;
                yj = mesh.nodes{j}.y;
                l = sqrt((xi - xj)^2 + (yi - yj)^2);
                if strcmp(load.type{k}, '作用均布载荷')
                    Pq = load.q(k) / 2 * [yi - yj; xj - xi; yi - yj; xj - xi];
                elseif strcmp(load.type{k},'作用均布切向载荷')
                    Pq = load.q(k) / 2 * [xj - xi; yj - yi; xj - xi; yj - yi];
                elseif strcmp(load.type{k}, 'x方向三角形均布载荷')
                    Pq = load.q(k) / 2 * l * [2 / 3; 0; 1 / 3; 0];
                elseif strcmp(load.type{k}, 'x方向均布载荷')
                    Pq = load.q(k) / 2 * l * [1; 0; 1; 0];
                end
                F = F + B' * Pq;
            end
            obj.F = F;
        end
        function displacement_condition(obj, dbc, a)
            if nargin < 3
                a = 1e10;
            end
            for i = 1:dbc.num
                if strcmp(dbc.type{i}, 'x位移')
                    obj.K(2 * dbc.node(i) - 1, 2 * dbc.node(i) - 1) = obj.K(2 * dbc.node(i) - 1, 2 * dbc.node(i) - 1) * a;
                    obj.F(2 * dbc.node(i) - 1) = obj.K(2 * dbc.node(i) - 1, 2 * dbc.node(i) - 1) * dbc.d(i);
                elseif strcmp(dbc.type{i}, 'y位移')
                    obj.K(2 * dbc.node(i), 2 * dbc.node(i)) = obj.K(2 * dbc.node(i), 2 * dbc.node(i)) * a;
                    obj.F(2 * dbc.node(i)) = obj.K(2 * dbc.node(i), 2 * dbc.node(i)) * dbc.d(i);
                end
            end
        end

        function d_m = solve_displacement_external_node(obj)
                d_m = obj.K \ obj.F;
                obj.d_m = d_m;
        end
        function solve_displacement_internal_node(obj,mesh)
            for i = 1:mesh.element_num
                element = mesh.elements{i};
                element.d_m = zeros(element.edge_m_num*2,1);
                element.d_c = zeros(element.edge_c_num*2,1);
                for j = 1:element.edge_m_num
                    element.d_m(2*j-1:2*j) = obj.d_m(2*element.node_m_id(j)-1:2*element.node_m_id(j));
                end
                temp_d_c = -inv(element.Ke22)*(element.Ke21*element.d_m);
                element.d_c = temp_d_c(1:2*element.edge_c_num);
            end
        end
        function sigma = get_sigma(obj,x,y,beta,P_type)
            if P_type==7
                P7 = [1,y,0,0,0,x,0;
                      0,0,1,x,0,0,y;
                      0,0,0,0,1,-y,-x];
                sigma = P7*beta;
            elseif P_type==12
                P12 = [1,y,0,0,0,x,0,y^2,0,x^2,x*y,0;
                       0,0,1,x,0,0,y,0,x^2,y^2,0,x*y;
                       0,0,0,0,1,-y,-x,0,0,-2*x*y,-y^2/2,-x^2/2];
                sigma = P12*beta;
            elseif P_type==18
                P18 = [1,y,0,0,0,x,0,y^2,0,x^2,x*y,0,x^3,y^3,3*x^2*y,x*y^2,0,0;
                       0,0,1,x,0,0,y,0,x^2,y^2,0,x*y,3*x*y^2,0,y^3,0,x^2*y,x^3;
                       0,0,0,0,1,-y,-x,0,0,-2*x*y,-y^2/2,-x^2/2,-3*x^2*y,0,-3*x*y^2,-y^3/3,-x^3/3,0];
                sigma = P18*beta;
            end
        end

        function strain = get_strain(obj,x,y,beta,type,P_type)
            sigma = obj.get_sigma(x,y,beta,P_type);
            if strcmp(type,'m')
               E = obj.E_m;
               pr = obj.pr_m;
            elseif strcmp(type,'c')
               E = obj.E_c;
               pr = obj.pr_c;
            end
            S = 1/E*[1,-pr,0;
                -pr,1,0;
                 0,0,2*(1+pr)];
            strain = S*sigma;
        end

        function integral = sigma_integral(obj,x1,y1,x2,y2,x3,y3,J,beta,P_type)
            sigma_integral = zeros(3,1);
            if P_type ==7
                F1 = obj.get_sigma(obj.transform(x1, x2, x3, 1/3, 1/3), obj.transform(y1, y2, y3, 1/3, 1/3), beta,P_type);
                F2 = obj.get_sigma(obj.transform(x1, x2, x3, 3/5, 1/5), obj.transform(y1, y2, y3, 3/5, 1/5), beta,P_type);
                F3 = obj.get_sigma(obj.transform(x1, x2, x3, 1/5, 1/5), obj.transform(y1, y2, y3, 1/5, 1/5), beta,P_type);
                F4 = obj.get_sigma(obj.transform(x1, x2, x3, 1/5, 3/5), obj.transform(y1, y2, y3, 1/5, 3/5), beta,P_type);
                sigma_integral = 0.5 * J * (-27/48 * F1 + 25/48 * F2 + 25/48 * F3 + 25/48 * F4);  % 3阶代数精度
            elseif P_type ==12
                r = [0.44594849, 0.09157621, 0.10810302, 0.44594849, 0.81684757,0.09157621];
                s = [0.44594849, 0.09157621, 0.44594849, 0.10810302, 0.09157621,0.81684757];
                half_weights = [0.11169079, 0.05497587, 0.11169079, 0.11169079, 0.05497587,0.05497587];
                for i=1:6
                    sigma_integral = sigma_integral + half_weights(i)*J*obj.get_sigma(obj.transform(x1, x2, x3, r(i), s(i)), obj.transform(y1, y2, y3, r(i), s(i)), beta,P_type);
                end
            elseif P_type==18
                r = [0.33333333, 0.26034597, 0.0651301 , 0.3128655 , 0.47930807,...
                    0.26034597, 0.86973979, 0.0651301 , 0.63844419, 0.63844419,...
                    0.04869032, 0.04869032, 0.3128655 ];
                s = [0.33333333, 0.26034597, 0.0651301 , 0.63844419, 0.26034597,...
                    0.47930807, 0.0651301 , 0.86973979, 0.3128655 , 0.04869032,...
                     0.3128655 , 0.63844419, 0.04869032];
                half_weights = [-0.07478502,  0.08780763,  0.02667362,  0.03855688,  0.08780763,...
                                 0.08780763,  0.02667362,  0.02667362,  0.03855688,  0.03855688,...
                                 0.03855688,  0.03855688,  0.03855688];
                for i=1:13
                    sigma_integral = sigma_integral + half_weights(i)*J*obj.get_sigma(obj.transform(x1, x2, x3, r(i), s(i)), obj.transform(y1, y2, y3, r(i), s(i)), beta,P_type);
                end              
            end 
            integral = sigma_integral;
        end

        function [total_sigma_integral,total_strain_integral,total_area] = calculate_average_stress(obj,mesh)
            total_sigma_integral = zeros(3,1);
            total_strain_integral = zeros(3,1);
            obj.total_area_m=0;
            obj.total_area_c=0;
            total_area = 0;
            S_m = 1/obj.E_m*[1,-obj.pr_m,0;
                -obj.pr_m,1,0;
                 0,0,2*(1+obj.pr_m)];
            S_c = 1/obj.E_c*[1,-obj.pr_c,0;
                -obj.pr_c,1,0;
                 0,0,2*(1+obj.pr_c)];
            for i = 1:mesh.element_num
                element = mesh.elements{i};
                element.beta_m = element.H_m\element.G_mm*element.d_m-element.H_m\element.G_mc*element.d_c;
                element.beta_c = element.H_c\element.G_cc*element.d_c;
                x1 = element.node_m(1);%按照变形前的位型积分
                y1 = element.node_m(2);
                element.area_m = 0;
                element.area_c = 0;
                element.sigma_m_integral = zeros(3,1);
                element.sigma_c_integral = zeros(3,1);
                element.strain_m_integral = zeros(3,1);
                element.strain_c_integral = zeros(3,1);
                for j = 2:element.edge_m_num-1
                    x2 = element.node_m(2 * (j - 1) + 1);
                    y2 = element.node_m(2 * (j - 1) + 2);
                    x3 = element.node_m(2 * j + 1);
                    y3 = element.node_m(2 * j + 2);
                    J = obj.Jacobian(x1, y1, x2, y2, x3, y3);
                    element.area_m = element.area_m + J*0.5;
                    if J == 0
                        continue;
                    end
                    triangle_integral = obj.sigma_integral(x1,y1,x2,y2,x3,y3,J,element.beta_m,element.n_beta_M);
                    element.sigma_m_integral = element.sigma_m_integral + triangle_integral;
                    element.strain_m_integral = element.strain_m_integral + S_m*triangle_integral;
                end
                x1 = element.node_c(1);
                y1 = element.node_c(2);
                for j = 2:element.edge_c_num-1
                    x2 = element.node_c(2 * (j - 1) + 1);
                    y2 = element.node_c(2 * (j - 1) + 2);
                    x3 = element.node_c(2 * j + 1);
                    y3 = element.node_c(2 * j + 2);
                    J = obj.Jacobian(x1, y1, x2, y2, x3, y3);
                    element.area_m = element.area_m-J*0.5;
                    element.area_c = element.area_c+J*0.5;
                    if J == 0
                        continue;
                    end
                    triangle_integral = obj.sigma_integral(x1,y1,x2,y2,x3,y3,J,element.beta_m,element.n_beta_M);
                    element.sigma_m_integral = element.sigma_m_integral - triangle_integral;
                    element.strain_m_integral = element.strain_m_integral - S_m*triangle_integral;
                    triangle_integral = obj.sigma_integral(x1,y1,x2,y2,x3,y3,J,element.beta_c,element.n_beta_C);
                    element.sigma_c_integral = element.sigma_c_integral + triangle_integral;
                    element.strain_c_integral = element.strain_c_integral + S_c*triangle_integral;
                end
                total_area = total_area + element.area_m + element.area_c;
                obj.total_area_m = obj.total_area_m+element.area_m;
                obj.total_area_c = obj.total_area_c+element.area_c;
                total_sigma_integral = total_sigma_integral + element.sigma_m_integral + element.sigma_c_integral;
                total_strain_integral = total_strain_integral + element.strain_m_integral + element.strain_c_integral;
            end
        end
        function [effective_E,effective_pr] = get_effective_modulus(obj,average_stress,average_strain)
            sigma_x = average_stress(1);
            sigma_y = average_stress(2);
            epsilon_x = average_strain(1);
            epsilon_y = average_strain(2);
            effective_E = (sigma_x^2-sigma_y^2)/(sigma_x*epsilon_x-sigma_y*epsilon_y);
            effective_pr = (sigma_y*epsilon_x-sigma_x*epsilon_y)/(sigma_x*epsilon_x-sigma_y*epsilon_y);
        end

    end
end

