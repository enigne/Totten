steps = [3];

%Load some necessary codes {{{
glacier = 'Totten'; hem='s';
clustername = oshostname();

if strcmpi(oshostname(), 'totten')
	prepath = '/totten_1/chenggong/';
else
	prepath = '/Users/gongcheng/Research/';
end
addpath([prepath, glacier, '/PostProcessing/']);
addpath([prepath, glacier, '/src/']);
projPath = [prepath, glacier, '/'];
% }}}
%Cluster parameters{{{
if strcmpi(clustername, 'andes')
	cluster=andes('numnodes',1,'cpuspernode',64, 'memory', 32);
	cluster.time = jobTime;
	waitonlock = 0;
elseif strcmpi(clustername, 'frontera')
	%cluster=frontera('numnodes', 1,'cpuspernode',56,'queue','flex');
	cluster=frontera('numnodes',3,'cpuspernode',56);
	cluster.time = jobTime;
	waitonlock = 0;
else
	cluster=generic('name',oshostname(),'np', 64);
	waitonlock = Inf;
end
clear clustername
org=organizer('repository',[projPath, '/Models'],'prefix',['Model_' glacier '_'],'steps',steps); clear steps;
fprintf(['\n  ========  ' upper(glacier) '  ========\n\n']);
%}}}

% Step 1--5
if perform(org, 'Mesh')% {{{
	%refine mesh using surface velocities as metric
	if strcmp(hem,'n')
		md=triangle(model,['Exp/' glacier '.exp'],500);
		[velx, vely]=interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
		vel  = sqrt(velx.^2+vely.^2);

		md=bamg(md,'hmin',100,'hmax',2500,'field',vel,'err',5);
		[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
		md.mesh.epsg=3413;
	else
		md=triangle(model,['Exp/Totten.exp'],5e3);
		coarse=15e3;
		fine_vel=1500;

		for i=1:2,
			disp(['--- Performing static mesh adaptation. Step ' num2str(i) '/4']);
			% using a priori analysis (observed velocity)
			disp('   -- Interpolating some data');
			[velx vely] = interpMouginotAnt2017(md.mesh.x,md.mesh.y);
			surface=interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'surface');
			ocean_levelset=-ones(size(md.mesh.x));% all floating
			ocean_levelset(find(interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'mask')==2))=1; % grounded from BedMachine
			ice_levelset=-ones(size(md.mesh.x));
			ice_levelset(find(interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'mask')==0))=1; % grounded from BedMachine

			pos=find(isnan(velx) | isnan(vely) | ice_levelset>0);% | ocean_levelset<0);
			velx(pos)=0; vely(pos)=0; vel=sqrt(velx.^2+vely.^2);

			hVertices = NaN(md.mesh.numberofvertices,1);
			hVertices(find(vel>200)) = fine_vel;
			md=bamg(md,'gradation',1.6180,'anisomax',1.e6,'KeepVertices',0,'Hessiantype',0,'Metrictype',0,...
				'hmax',coarse,'hmin',fine_vel,'hVertices',hVertices,'field',vel,'err',3);

			md.private.bamg=struct();
		end

		[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,-1);
		md.mesh.epsg=3031;
		md.mesh.scale_factor=(1+sin(md.mesh.lat*pi/180))/(1+sin(-71*pi/180));
	end

	savemodel(org,md);
end %}}}
if perform(org, ['Param'])% {{{

	md=loadmodel(org,'Mesh');
	md=setflowequation(md,'SSA','all');

	if strcmp(hem,'n')
		md=setmask(md,'','');
		md=parameterize(md,'Greenland.par');
	else
		md=setmask(md,'','');
		md=parameterize(md,'Antarctica.par');
	end
	savemodel(org,md);
end%}}}
if perform(org, ['Inversion_drag_Weertman'])% {{{

	md=loadmodel(org, ['Param']);

	costcoeffs = [20000, 2, 0.5e-8];
	% Set the friction law to Weertman
	md.friction=frictionweertman();
	md.friction.m = 3.0*ones(md.mesh.numberofelements,1);
	md.friction.C=EstimateFric_Weertman(md); % initial guess from Driving Stress (using m=3)

	%No friction on PURELY ocean element
	pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
	flags=ones(md.mesh.numberofvertices,1);
	flags(md.mesh.elements(pos_e,:))=0;
	md.friction.C(find(flags))=0.0;

	% also set floating ice friction to 0.0
	pos=find(md.mask.ocean_levelset<0);
	md.friction.C(pos) = 0.0;

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);
	md.transient.amr_frequency = 0;

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=zeros(md.mesh.numberofvertices,numel(md.inversion.cost_functions));
	md.inversion.cost_functions_coefficients(:,1)=costcoeffs(1);
	md.inversion.cost_functions_coefficients(:,2)=costcoeffs(2);
	md.inversion.cost_functions_coefficients(:,3)=costcoeffs(3);
	pos=find(md.mask.ice_levelset>0);
	md.inversion.cost_functions_coefficients(pos,1:2)=0;

	%Controls
	md.inversion.control_parameters={'FrictionC'};
	md.inversion.maxsteps=400;
	md.inversion.maxiter =400;
	md.inversion.min_parameters=0*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=5e4*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;
	md.inversion.dxmin = 0.01;
	md.inversion.gttol = 1e-6;
	%Additional parameters
	md.stressbalance.restol=1e-3;
	md.stressbalance.reltol=1e-2;
	md.stressbalance.abstol=10;

	md.settings.solver_residue_threshold=NaN;

	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	%Go solve
	md.cluster=cluster;
	md=solve(md,'sb');

	%Put results back into the model
	md.friction.C=md.results.StressbalanceSolution.FrictionC;
	md.initialization.vx=md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;

	savemodel(org,md);
end%}}}
