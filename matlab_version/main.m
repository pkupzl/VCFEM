run_example = 0;
if run_example==0
   opt = Initialize();
elseif run_example==1
    opt=Initialize_from_csv('D:/VCFEM_dataset/original_image/0/',1/1024,'上边左切');
end
mesh = Mesh(opt.node_num,...
            opt.element_num,...
            opt.nodes,...
            opt.edge_m_nums,...
            opt.edge_c_nums,...
            opt.node_ms,...
            opt.node_cs,...
            opt.node_m_ids,...
            opt.node_c_ids);
load = Load(opt.load_num,opt.load_type,opt.q,opt.node_i,opt.node_j);
dbc = Displacement_BC(opt.dbc_num,opt.dbc_node,opt.dbc_type,opt.d);
vcfem = VCFEM(opt.E_m,opt.E_c,opt.pr_m,opt.pr_c);
Ke = vcfem.assembly_global_stiffness_matrix(mesh);
F = vcfem.calculate_global_nodal_load(mesh,load);
vcfem.displacement_condition(dbc,10^10);
K = vcfem.K;
d_m = vcfem.solve_displacement_external_node();
vcfem.solve_displacement_internal_node(mesh);
[total_sigma_integral,total_strain_integral,total_area] = vcfem.calculate_average_stress(mesh);
fprintf('average sigma_x = %f\n', total_sigma_integral(1)/total_area);
fprintf('average sigma_y = %f\n', total_sigma_integral(2)/total_area);
fprintf('average tau_xy = %f\n', total_sigma_integral(3)/total_area);
fprintf('average epsilon_x = %f\n', total_strain_integral(1)/total_area);
fprintf('average epsilon_y = %f\n', total_strain_integral(2)/total_area);
fprintf('average epsilon_xy = %f\n', total_strain_integral(3)/total_area);
[effective_E,effective_pr] = vcfem.get_effective_modulus(total_sigma_integral/total_area,total_strain_integral/total_area);
fprintf('effective E = %f\n', effective_E);
fprintf('effective pr = %f\n', effective_pr);
if run_example==0
    Visualize('g--',mesh, opt.nodes, opt.topPoints, opt.bottomPoints, opt.leftPoints, opt.rightPoints);
    Visualize('b-',mesh, opt.nodes, opt.topPoints, opt.bottomPoints, opt.leftPoints, opt.rightPoints,d_m,opt.particle_nodes);
end
if run_example==1
    Visualize('g--',mesh, opt.nodes, opt.topPoints, opt.bottomPoints, opt.leftPoints, opt.rightPoints);
    Visualize('b-',mesh, opt.nodes, opt.topPoints, opt.bottomPoints, opt.leftPoints, opt.rightPoints,d_m,opt.particle_nodes);
end