function varargout=runme(varargin)

	%Check inputs {{{
	if nargout>1
		help runme
		error('runme error message: bad usage');
	end
	%recover options
	options=pairoptions(varargin{:});
	% }}}
	%GET cluster name: totten{{{
	clustername = getfieldvalue(options,'cluster name','totten');
	% }}}
	%GET steps {{{
	steps = getfieldvalue(options,'steps',[1]);
	% }}}
	%GET Friction type: 'Schoof' {{{
	friction = getfieldvalue(options,'friction','Schoof');
	if ~ismember(friction,{'Schoof','Weertman','Budd'})
		disp('runme warning: friction law not supported, defaulting to use ''Schoof law''')
		friction = 'Schoof';
	end
	flagFriction = 0; % 0-by default,schoof, 1-Weertman, 2-Budd

	if strcmp(friction, 'Schoof')
		flagFriction = 0;
	elseif strcmp(friction, 'Weertman')
		flagFriction = 1;
	elseif strcmp(friction, 'Budd')
		flagFriction = 2;
	end
	% }}}
	%GET calving law: Obs{{{
	calving = getfieldvalue(options,'calving','Obs');
	if ~ismember(calving,{'VM', 'Obs', 'Greene', 'CALFIN', 'TermPicks'})
		disp('runme warning: calving law not supported, defaulting to use ''Observed ice front positions''')
		calving = 'Obs';
	end

	% set default obs
	if strcmp(calving, 'Obs')
		calving = 'Greene';
	end

	flagObsCalving = 1;
	if strcmp(calving, 'VM')
		flagObsCalving = 0;
	end
	% }}}
	%GET SMB model: MAR or RACMO{{{
	smb_model = getfieldvalue(options,'smb model','MAR');
	% }}}
	%GET rerun_inversion: 0{{{
	rerun_inversion = getfieldvalue(options,'rerun inversion', 0);
	% }}}
	%GET sigma: 1.0{{{
	sigma = getfieldvalue(options,'sigma', 1.0);
	% }}}
	%GET Cmax: 0.8 {{{
	Cmax = getfieldvalue(options,'Cmax',0.8);
	% }}}
	%GET savepath: '/'{{{
	savePath = getfieldvalue(options,'savePath', './testrun');
	% }}}
	%GET stabilization for levelset: 5 - SUPG {{{
	levelsetStabilization = getfieldvalue(options,'levelset stabilization', 5);
	% }}}
	%GET reinitialization for levelset: 50 {{{
	levelsetReinit = getfieldvalue(options,'levelset reinitialize', 50);
	% }}}
	%GET finalTime: 2020 {{{
	finalTime = getfieldvalue(options,'finalTime', 2022);
	% }}}
	%GET residue threshold for linear solver: 1e-6{{{
	residue_threshold = getfieldvalue(options,'residue threshold', 1e-6);
	% }}}
	%GET jobTime: 5 (hours){{{
	jobTime = getfieldvalue(options,'jobTime', 5);
	% }}}
	%GET cost coefficients: [1,1,1e-6] {{{
	costcoeffs = getfieldvalue(options,'cost coefficients', [1,1,1e-6]); 
	% }}}
	%GET select subregion: [2,2,1] {{{
	subregion = getfieldvalue(options,'select subregion', [1,1,1e-6]); 
	% }}}

	%Load some necessary codes {{{
	glacier = 'Totten'; hem='s';

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

		if (rerun_inversion)
			disp(['  Rerun inversion using previous results as initial guess'])
			md=loadmodel(org, ['Inversion_drag_Weertman']);
		else
			md=loadmodel(org, ['Param']);

			% Set the friction law to Weertman
			md.friction=frictionweertman();
			md.friction.m = 3.0*ones(md.mesh.numberofelements,1);
			md.friction.C=EstimateFric_Weertman(md); % initial guess from Driving Stress (using m=3)
		end

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
	if perform(org, ['Inversion_drag_Budd'])% {{{
		if (rerun_inversion)
			disp(['  Rerun inversion using previous results as initial guess'])
			md=loadmodel(org, ['Inversion_drag_Budd']);
		else
			md=loadmodel(org, 'InversionB');
		end
		% fixed after add this option to effective pressure coupling
		md.friction.coupling = 2;

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
		md.inversion.control_parameters={'FrictionCoefficient'};
		md.inversion.maxsteps = 200;
		md.inversion.maxiter = 200;
		md.inversion.min_parameters=0*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=1e4*ones(md.mesh.numberofvertices,1);
		md.inversion.control_scaling_factors=1;
		md.inversion.dxmin = 0.01;
		md.inversion.gttol = 1e-6;

		%Additional parameters
		md.stressbalance.restol=1e-6;
		md.stressbalance.reltol=1e-4;
		md.stressbalance.abstol=NaN;

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		%Go solve
		md.cluster=cluster;
		md=solve(md,'sb');

		%Put results back into the model
		md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
		md.initialization.vx=md.results.StressbalanceSolution.Vx;
		md.initialization.vy=md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end%}}}
	varargout{1} = md;
	return;

	if perform(org, ['Inversion_drag_ISMIP', damage_suffix, '_Schoof'])% {{{

		if (rerun_inversion)
			disp(['  Rerun inversion using previous results as initial guess'])
			md=loadmodel(org, ['Inversion_drag_ISMIP', damage_suffix, '_Schoof']);
		else
			md=loadmodel(org, ['Param_ISMIP', damage_suffix]);

			% Set the friction law to schoof's
			md.friction=frictionschoof();
			md.friction.m = 1.0/3.0*ones(md.mesh.numberofelements,1);
			md.friction.Cmax = Cmax*ones(md.mesh.numberofvertices,1);
			md.friction.C = 2000*ones(md.mesh.numberofvertices,1);
			md.friction.coupling = 2;
		end

		%No friction on PURELY ocean element
		pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
		flags=ones(md.mesh.numberofvertices,1);
		flags(md.mesh.elements(pos_e,:))=0;
		md.friction.C(find(flags))=1e-2;

		fill_pos = find(min(md.initialization.vel(md.mesh.elements),[],2)==0); % fill-in postions
		ids = zeros(md.mesh.numberofvertices, 1);
		ids(md.mesh.elements(fill_pos,:))=1;

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
		pos=find((md.mask.ice_levelset>0) );
		md.inversion.cost_functions_coefficients(pos,1:2)=0;

		% set C=1e-2 for no obs region
		md.friction.C(ids>0) = 1e-2;
		md.inversion.cost_functions_coefficients(ids>0,1:2)=0;

		%Controls
		md.inversion.control_parameters={'FrictionC'};
		md.inversion.maxsteps=500;
		md.inversion.maxiter =500;
		md.inversion.min_parameters=1e-9*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=1e5*ones(md.mesh.numberofvertices,1);
		md.inversion.control_scaling_factors=3000;
		md.inversion.dxmin = 1e-6;
		%Additional parameters
		md.stressbalance.restol=1e-5;
		md.stressbalance.reltol=1e-3;
		md.stressbalance.abstol=NaN;
		md.stressbalance.maxiter = 100;

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		%Go solve
		md.cluster=cluster;
		md=solve(md,'sb');

		% get a direction to update cost coeff
		newCoeff = updateCostCoeff(md, ratio);
		disp(sprintf('With the given ratio %d, %d, %d \n', ratio))
		disp(sprintf('The coefficients can be updated to %g,   %g,   %g \n', newCoeff))

		%Put results back into the model
		md.friction.C=md.results.StressbalanceSolution.FrictionC;
		md.initialization.vx=md.results.StressbalanceSolution.Vx;
		md.initialization.vy=md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end%}}}

	% step 6--10
	if perform(org, ['MARERA5smb_ISMIP', suffix])% {{{
		md=loadmodel(org,['Inversion_drag_ISMIP', suffix]);

		filelist = {...
			'MARv3.11-monthly-ERA5-2007.nc'
		'MARv3.11-monthly-ERA5-2008.nc'
		'MARv3.11-monthly-ERA5-2009.nc'
		'MARv3.11-monthly-ERA5-2010.nc'
		'MARv3.11-monthly-ERA5-2011.nc'
		'MARv3.11-monthly-ERA5-2012.nc'
		'MARv3.11-monthly-ERA5-2013.nc'
		'MARv3.11-monthly-ERA5-2014.nc'
		'MARv3.11-monthly-ERA5-2015.nc'
		'MARv3.11-monthly-ERA5-2016.nc'
		'MARv3.11-monthly-ERA5-2017.nc'
		'MARv3.11-monthly-ERA5-2018.nc'
		'MARv3.11-monthly-ERA5-2019.nc'
		};

		md.smb.mass_balance = [];
		T0 = 2007;
		for i=1:numel(filelist)
			filename = ['/totten_1/ModelData/Greenland/MARv3.11-ERA5/' filelist{i}];
			disp(['Loading ' filename]);

			%Get time
			T  = double(ncread(filename,'time'));
			time = T0 + (i-1) + (T+0.5)/12;

			%Coordinates
			x=double(ncread(filename,'x'));
			y=double(ncread(filename,'y'));

			%Load SMB for this year
			SMB=double(ncread(filename,'SMB'))/1000.0*12.0*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/month to mIE/yr

			%Interpolate
			for m=1:numel(time)
				smb_mesh=InterpFromGridToMesh(x, y, SMB(:,:,m)', md.mesh.x, md.mesh.y, 0);
				md.smb.mass_balance = [md.smb.mass_balance, [smb_mesh;time(m)]];
			end
		end
		%Clean up
		savemodel(org,md);
	end%}}}
	if perform(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]),% {{{

		md=loadmodel(org, ['MARERA5smb_ISMIP', suffix]);

		md.initialization.pressure = zeros(md.mesh.numberofvertices,1); %FIXME

		% spc thickness on the boundary
		md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);
		pos=find((md.mask.ice_levelset<0).*(md.mesh.vertexonboundary));
		md.masstransport.spcthickness(pos) = md.geometry.thickness(pos);

		% Set parameters
		md.inversion.iscontrol=0;
		md.timestepping.start_time = 2007;
		md.timestepping.time_step  = 0.005;
		md.timestepping.final_time = finalTime;
		md.settings.output_frequency = 1;

		md.transient.ismovingfront=1;
		md.transient.isslc = 0;
		md.transient.isthermal=0;
		md.transient.isstressbalance=1;
		md.transient.ismasstransport=1;
		md.transient.isgroundingline=1;
		md.groundingline.migration = 'SubelementMigration';

		% spclevelset
		disp('Assigning Calving front');
		% initial condition is from the geometry
		pos=find(md.mask.ice_levelset<0); md.mask.ice_levelset(pos)=-1;
		pos=find(md.mask.ice_levelset>0); md.mask.ice_levelset(pos)=+1;

		if flagObsCalving
			switch calving
				case 'CALFIN'
					disp('Using CALFIN ice front--->')
					initLevelset = reinitializelevelset(md, md.mask.ice_levelset);
					% at time step for the spc
					initLevelset = [initLevelset;md.timestepping.start_time+md.timestepping.time_step];
					% get observed calving front
					CFcontour = [projPath, '/DATA/Helheim_Calving_Front_2007_2020/termini_2007_2020.shp'];
					distance = ExpToLevelSet(md.mesh.x, md.mesh.y, CFcontour);
					% transient spc
					md.levelset.spclevelset = [initLevelset, distance];
				case 'Greene'
					disp('Using Greene monthly ice front--->')
					% get observed calving front
					mask = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [md.timestepping.start_time, md.timestepping.final_time]);
					% convert icemask to levelset distance
					distance = zeros(size(mask));
					for i = 1:size(mask,2)
						distance(1:end-1,i) = reinitializelevelset(md, mask(1:end-1,i));
					end
					distance(end,:) = mask(end,:);
					% transient spc
					md.levelset.spclevelset = distance;
				case 'TermPicks'
					disp('Using TermPicks ice front--->')
					% get observed calving front
					CFcontour = [projPath, '/DATA/TermPicks_closed.shp'];
					distance = ExpToLevelSet(md.mesh.x, md.mesh.y, CFcontour);
					% HOTFIX:remove id=591, time=2011.3945
					id = find(distance(end,:)<2011.394 | distance(end,:)>2011.395);
					% transient spc
					md.levelset.spclevelset = distance(:,id);
				otherwise
					error('Unknown obs calving data set')
			end
			disp(['Load observed calving front: ', num2str(size(md.levelset.spclevelset,2)), ...
				' records from ', datestr(md.levelset.spclevelset(end,1)*365.25, 'yyyy-mm-dd'), ' to ',...
				datestr(md.levelset.spclevelset(end,end)*365.25, 'yyyy-mm-dd')]);
		else
			% only set boundary conditions, so that levelset can solve for the calving front
			md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
			pos = find(md.mesh.vertexonboundary);
			md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
		end
		% only set boundary conditions, so that levelset can solve for the calving front
		md.levelset.stabilization = levelsetStabilization;
		disp(['  Levelset function uses stabilization ', num2str(md.levelset.stabilization)]);
		md.levelset.reinit_frequency = levelsetReinit;
		disp(['  Levelset function reinitializes every ', num2str(md.levelset.reinit_frequency), ' time steps']);

		savemodel(org,md);
	end%}}}
	if perform(org, ['Transient_ERA5_ISMIP', calving_suffix, suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_floatingice = 200*10^3;
		md.calving.stress_threshold_groundedice = sigma*10^6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:

		% meltingrate
		timestamps = [md.timestepping.start_time, md.timestepping.final_time];
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		md.cluster = cluster;
		md.verbose.solution = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		md.transient.requested_outputs={'default','SmbMassBalance', 'IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org, ['Transient_ERA5_ISMIP_spcvel', calving_suffix, suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_floatingice = 200*10^3;
		md.calving.stress_threshold_groundedice = sigma*10^6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:

		% meltingrate
		timestamps = [md.timestepping.start_time, md.timestepping.final_time];
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		% load obs vel to constrain the solutions
		obsdatafile = [projPath, './PostProcessing/Results/timeSeries_Obs_mapped.mat'];
		disp(['Loading obs vel from ', obsdatafile])
		load(obsdatafile);
		% save the current spcs for BC.
		bcvx = md.stressbalance.spcvx;
		bcvy = md.stressbalance.spcvy;
		pos = find(~isnan(md.stressbalance.spcvx));
		md.stressbalance.spcvx = [vx_obs;time];
		md.stressbalance.spcvy = [vy_obs;time];

		Nt = numel(time);
		md.stressbalance.spcvx(pos,:) = repmat(bcvx(pos,1), 1, Nt);
		md.stressbalance.spcvy(pos,:) = repmat(bcvy(pos,1), 1, Nt);

		% very important, to not interpolate in time
		md.timestepping.interp_forcing = 0;

		md.cluster = cluster;
		md.verbose.solution = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org, ['Transient_ERA5_ISMIP_Arate_driven', suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% set the first 0.5 years levelset to spc levelset, but the rest to be NaN except boundary
		relaxTid = find(md.levelset.spclevelset(end, :)>md.timestepping.start_time+relaxtime);
		disp(['  Set levelset to NaN after ', num2str(md.levelset.spclevelset(end, relaxTid(1)))]);
		md.levelset.spclevelset(1:end-1, relaxTid) = NaN;
		% only set boundary conditions, so that levelset can solve for the calving front
		pos = find(md.mesh.vertexonboundary);
		md.levelset.spclevelset(pos, :) = md.mask.ice_levelset(pos)*ones(1,length(md.levelset.spclevelset(end,:)));

		% only set boundary conditions, so that levelset can solve for the calving front
		md.levelset.stabilization = levelsetStabilization;
		disp(['  Levelset function uses stabilization ', num2str(md.levelset.stabilization)]);
		md.levelset.reinit_frequency = levelsetReinit;
		disp(['  Levelset function reinitializes every ', num2str(md.levelset.reinit_frequency), ' time steps']);

		% calving parameters
		md.calving = calvingparameterization();
		md.calving.use_param = -1;

		aRateFile = [projPath, '/PostProcessing/Results/Arates_', arateSource, '_projected_aver',num2str(meanwindow),'.mat'];
		disp(['Loading ablation rate from: ', aRateFile]);
		cm = load(aRateFile);
		timesteps = cm.time;
		aRate = cm.aRateOnMesh;

		% meltingrate
		md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices+1, length(cm.time));
		md.frontalforcings.meltingrate(end,:) = timesteps;

		md.frontalforcings.ablationrate = md.frontalforcings.meltingrate;
		md.frontalforcings.ablationrate(1:end-1,:) = aRate;

		md.cluster = cluster;
		md.verbose.solution = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}

	% step 11--15
	if perform(org, ['Transient_ISMIP_parameterization', suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% set the first 1 years levelset to spc levelset, but the rest to be NaN except boundary
		relaxTid = find(md.levelset.spclevelset(end, :)>md.timestepping.start_time+relaxtime);
		disp(['  Set levelset to NaN after ', num2str(md.levelset.spclevelset(end, relaxTid(1)))]);
		md.levelset.spclevelset(1:end-1, relaxTid) = NaN;
		% only set boundary conditions, so that levelset can solve for the calving front
		pos = find(md.mesh.vertexonboundary);
		md.levelset.spclevelset(pos, :) = md.mask.ice_levelset(pos)*ones(1,length(md.levelset.spclevelset(end,:)));

		% only set boundary conditions, so that levelset can solve for the calving front
		md.levelset.stabilization = levelsetStabilization;
		disp(['  Levelset function uses stabilization ', num2str(md.levelset.stabilization)]);
		md.levelset.reinit_frequency = levelsetReinit;
		disp(['  Levelset function reinitializes every ', num2str(md.levelset.reinit_frequency), ' time steps']);

		% calving parameters
		md.calving = calvingparameterization();
		md.calving.use_param = use_param;
		md.calving.theta = theta;
		md.calving.alpha = alpha;
		md.calving.xoffset = xoffset;
		md.calving.yoffset = yoffset;
		disp(md.calving)

		aRateFile = [projPath, '/PostProcessing/Results/Arates_thd0_FC_', arateSource, '_Isoline_aver',num2str(meanwindow),'.mat'];
		disp(['  Loading ablation rate from: ', aRateFile]);
		cm = load(aRateFile);
		% meltingrate
		timesteps = cm.time;
		aRate = max(cm.aRateC);

		% meltingrate
		md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices+1, length(cm.time));
		md.frontalforcings.meltingrate(end,:) = timesteps;

		md.frontalforcings.ablationrate = md.frontalforcings.meltingrate;
		md.frontalforcings.ablationrate(1:end-1,:) = ones(md.mesh.numberofvertices,1)*aRate;

		md.cluster = cluster;
		md.settings.waitonlock = waitonlock; % do not wait for complete

		md.verbose.solution = 1;
		md.miscellaneous.name = [savePath];
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate', ...
			'StrainRateparallel', 'StrainRateperpendicular'};

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org, ['Transient_ERA5_ISMIP_calvingvm', suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% only set boundary conditions, so that levelset can solve for the calving front
		md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
		pos = find(md.mesh.vertexonboundary);
		md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);

		% only set boundary conditions, so that levelset can solve for the calving front
		md.levelset.stabilization = levelsetStabilization;
		disp(['  Levelset function uses stabilization ', num2str(md.levelset.stabilization)]);
		md.levelset.reinit_frequency = levelsetReinit;
		disp(['  Levelset function reinitializes every ', num2str(md.levelset.reinit_frequency), ' time steps']);


		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_groundedice = sigma*1e6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:
		md.calving.stress_threshold_floatingice = 0.5e6;

		% load TF and qsg
		datafile = [projPath, '/PreProcessing/ThermalForcing/TF_and_qsg.mat'];
		disp(['  Loading TF and qsg from ', datafile])
		tqdata = load(datafile);
		timestamps =tqdata.time;
		TF = max(0, (md.geometry.bed<0 & md.mesh.x>2.5e5 & md.mesh.y<-2.5e6)*tqdata.TF);
		qsg = max(0, tqdata.qsg);

		% set the melting parameterization
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;
		A = 3e-4;
		B = 0.15;
		Alpha = 0.39;
		Beta  = 1.18;
		BED = -repmat(min(0,md.geometry.bed),[1 numel(timestamps)]);
		md.frontalforcings.meltingrate(1:end-1,:) = ((A*BED.*qsg.^Alpha + B).*TF.^Beta) * 365;

		md.cluster = cluster;
		md.verbose.solution = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org, ['Transient_ISMIP_medianRatio', suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% set the first 1 years levelset to spc levelset, but the rest to be NaN except boundary
		relaxTid = find(md.levelset.spclevelset(end, :)>md.timestepping.start_time+relaxtime);
		disp(['  Set levelset to NaN after ', num2str(md.levelset.spclevelset(end, relaxTid(1)))]);
		md.levelset.spclevelset(1:end-1, relaxTid) = NaN;
		% only set boundary conditions, so that levelset can solve for the calving front
		pos = find(md.mesh.vertexonboundary);
		md.levelset.spclevelset(pos, :) = md.mask.ice_levelset(pos)*ones(1,length(md.levelset.spclevelset(end,:)));

		% only set boundary conditions, so that levelset can solve for the calving front
		md.levelset.stabilization = levelsetStabilization;
		disp(['  Levelset function uses stabilization ', num2str(md.levelset.stabilization)]);
		md.levelset.reinit_frequency = levelsetReinit;
		disp(['  Levelset function reinitializes every ', num2str(md.levelset.reinit_frequency), ' time steps']);

		maxbed = min(md.geometry.bed);

		% calving parameters
		md.calving = calvingparameterization();
		md.calving.use_param = use_param;
		md.calving.theta = theta;
		md.calving.alpha = 1./maxbed;
		md.calving.xoffset = xoffset;
		md.calving.yoffset = yoffset;
		disp(md.calving)

		aRateFile = [projPath, '/PostProcessing/Results/Arates_Reanalysis_Model_Isoline_aver',num2str(meanwindow),'.mat'];
		disp(['  Loading ablation rate from: ', aRateFile]);
		cm = load(aRateFile);
		% meltingrate
		timesteps = cm.time;
		aRate = median(cm.aRateC./cm.BedC, 'omitnan');

		% meltingrate
		md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices+1, length(cm.time));
		md.frontalforcings.meltingrate(end,:) = timesteps;

		md.frontalforcings.ablationrate = md.frontalforcings.meltingrate;
		md.frontalforcings.ablationrate(1:end-1,:) = maxbed* ones(md.mesh.numberofvertices,1)*aRate;

		md.cluster = cluster;
		md.settings.waitonlock = waitonlock; % do not wait for complete

		md.verbose.solution = 1;
		md.miscellaneous.name = [savePath];
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate', ...
			'StrainRateparallel', 'StrainRateperpendicular'};

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org,'Transient_Calibration') % {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% reset the final time
		md.timestepping.final_time = finalTime;

		%Prepare VEL obs =========================================================
		disp('Preparing observed velocities');
		obsData = load([projPath, './PostProcessing/Results/velObs_onmesh.mat']);

		timestamps = 0.5*(obsData.TStart+obsData.TEnd);
		pos = find((timestamps>=md.timestepping.start_time)&(timestamps<=md.timestepping.final_time));
		[~, ia, ic] = unique(timestamps(pos), 'stable');
		dataPos = pos(ia);

		disp(['  Found ', num2str(numel(dataPos)), ' observations between ', num2str(md.timestepping.start_time), ' and ', num2str(md.timestepping.final_time)])

		timestamps = timestamps(dataPos);
		vx_obs = obsData.vx_onmesh(:, dataPos);
		vy_obs = obsData.vy_onmesh(:, dataPos);

		%obsData = interpFromMEaSUREsGeotiff(md.mesh.x,md.mesh.y, timestamps(1), timestamps(2), 'glacier', 'Helheim');
		%vx_obs = [obsData(:).vx];
		%vy_obs = [obsData(:).vy];

		%fill in NaNs
		mask_obs=(~isnan(vx_obs) & ~isnan(vy_obs));
		for i=1:numel(timestamps)
			flags = isnan(vx_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			vx_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), vx_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');

			flags = isnan(vy_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			vy_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), vy_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');
		end

		%Prepare cost functions
		count = 1;
		weights = ones(md.mesh.numberofvertices,2);
		for i = 1:length(timestamps)

			%Filter out observations that are not within our simulation period
			if timestamps(i)<md.timestepping.start_time ||  timestamps(i)>md.timestepping.final_time
				continue
			end

			% 0.5*(Vx-Vx_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['VxMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Vx', 'observation_string','VxObs', 'observation', vx_obs(:, i)./md.constants.yts,...
				'weights',weights(:,1).*mask_obs(:,i), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;

			% 0.5*(Vy-Vy_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['VyMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Vy', 'observation_string','VyObs', 'observation', vy_obs(:, i)./md.constants.yts,...
				'weights',weights(:,2).*mask_obs(:,i), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;
		end

		%		%Prepare S obs =================================================================
		%		disp('Preparing Surface observations');
		%
		%		timestamps = [2007; 2012; 2017];
		%		z_obs = repmat(md.geometry.surface,[1 numel(timestamps)]);
		%		dhdt=interpShepherd2019(md.mesh.x, md.mesh.y, 'dhdt_2002_2006');
		%		z_obs(:,1) = z_obs(:,1) + dhdt*2;
		%		dhdt=interpShepherd2019(md.mesh.x, md.mesh.y, 'dhdt_2007_2011');
		%		z_obs(:,2) = z_obs(:,1) + dhdt*5;
		%		dhdt=interpShepherd2019(md.mesh.x, md.mesh.y, 'dhdt_2012_2016');
		%		z_obs(:,3) = z_obs(:,2) + dhdt*5;
		%
		%		%fill in NaNs
		%		mask_obs=(~isnan(z_obs));
		%		for i=1:numel(timestamps)
		%			flags = isnan(z_obs(:,i));
		%			pos1  = find(flags); pos2 = find(~flags);
		%			z_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), z_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');
		%		end
		%
		%		%Prepare cost functions
		%		weights = 1/(10000*md.constants.yts)*ones(md.mesh.numberofvertices,1);
		%		for i = 1:length(timestamps)
		%
		%			%Filter out observations that are not within our simulation period
		%			if timestamps(i)<md.timestepping.start_time ||  timestamps(i)>md.timestepping.final_time
		%				continue
		%			end
		%			% 0.5*(S-S_obs).^2
		%			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['DEMMisfit' num2str(count)],...
		%				'definitionstring',['Outputdefinition' num2str(count)],...
		%				'model_string','Surface', 'observation_string','SurfaceObservation', 'observation', z_obs(:, i),...
		%				'weights',weights.*mask_obs(:,i), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
		%			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
		%			count = count +1;
		%		end

		%Prepare Tikhonov ==============================================================
		disp('Preparing Tikhonov regularization');

		if isprop(md.friction, 'C')
			frictionCoeffName = 'FrictionC';
			tikhonovCoeff = 1e-12; % For Budd
		else
			frictionCoeffName = 'FrictionCoefficient';
			tikhonovCoeff = 1e-7; % For Budd
		end

		md.outputdefinition.definitions{count}=cfdragcoeffabsgrad('name','TikhonovReg',...
			'definitionstring',['Outputdefinition' num2str(count)],...
			'model_string', frictionCoeffName,...
			'weights',tikhonovCoeff*ones(md.mesh.numberofvertices,1),...
			'weights_string','WeightsSurfaceObservation');
		md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],...
			'type','scalar','fos_reverse_index',1);
		count = count +1;

		% Controls
		for i=1:1
			switch i
				case 1;
					if flagFriction < 2 % Schoof or Weertman
						field =md.friction.C; name = 'FrictionC'; scaling = 3000;
					else
						field =md.friction.coefficient; name = 'FrictionCoefficient'; scaling = 10;
					end
				case 2;
					field =md.materials.rheology_B;
					name = 'MaterialsRheologyBbar';
					scaling = 1e8;
				case 3;
					field =md.smb.mass_balance;
					name = 'SmbMassBalance';
					scaling = 1/md.constants.yts;
				otherwise error('not supported');
			end
			md.autodiff.independents{i} = independent('name',name,'type','vertex','nods',md.mesh.numberofvertices,...
				'control_size', size(field,2), 'min_parameters',1e-5*field, 'max_parameters',100*field, 'control_scaling_factor',scaling);
		end

		% General settings
		md.autodiff.driver='fos_reverse';
		md.autodiff.isautodiff=true;
		md.settings.checkpoint_frequency = 5;
		md.inversion=adm1qn3inversion(md.inversion);
		md.inversion.iscontrol=1;
		md.inversion.maxiter=30;
		md.inversion.maxsteps=md.inversion.maxiter;
		md.inversion.dxmin=1e-15;
		md.verbose.convergence=0;
		md.verbose=verbose('solution',true,'control',true);
		md.transient.requested_outputs={'default', 'IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset' };

		md.stressbalance.maxiter=20;

		%		md.timestepping.final_time = 2019.75;
		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_floatingice = 200*10^3;
		md.calving.stress_threshold_groundedice = sigma*10^6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:

		% meltingrate
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		md.cluster = cluster;
		md.settings.waitonlock = 0; % not wait for output 
		md.toolkits.DefaultAnalysis=issmmumpssolver();
		md.toolkits.RecoveryAnalysis=issmmumpssolver();
		%md.settings.solver_residue_threshold = 1.e-2;
		md.verbose.solution = 0;
		if strcmp(md.cluster.name, 'discovery')
			md.cluster.memory=500;
			md.cluster.interactive = 0; %only needed if you are using the generic cluster
		else
			md.cluster.interactive = 0; %only needed if you are using the generic cluster
		end
		md.miscellaneous.name = [savePath];

		md=solve(md,'Transient','runtimename',false);

		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org, ['Transient_FromAD', calving_suffix, suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% load AD results

		switch flagFriction
			case 0 % Schoof
				%ADFolder = '20230712_154351_Helheim_AD_Schoof/';
				%ADFolder = '20230825_151616_Helheim_AD_Schoof/'; % Schoof vel+Benn DEM
				ADFolder = '20230906_145535_Helheim_AD_Schoof/'; % spc H
			case 1 % Weertman
				ADFolder = '20230713_202014_Helheim_AD_Weertman/';
			case 2 % Budd
				ADFolder = '20230703_145405_Helheim_AD_Budd/';
			otherwise
				error('The friction in the setting is not used.');
		end


		ADPath = [projPath, '/Models/', ADFolder, 'Model_Helheim_Big_Transient'];
		disp(['loading friction coefficients from AD results in ', ADPath])
		md_AD = loadmodel(ADPath);
		if isprop(md.friction, 'C')
			md.friction.C = md_AD.results.TransientSolution(1).FrictionC;
		else
			md.friction.coefficient = md_AD.results.TransientSolution(1).FrictionCoefficient;
		end

		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_floatingice = 200*10^3;
		md.calving.stress_threshold_groundedice = sigma*10^6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:

		% meltingrate
		timestamps = [md.timestepping.start_time, md.timestepping.final_time];
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		md.cluster = cluster;
		md.verbose.solution = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}

	% step 16--20
	if perform(org, ['StartFrom1981', suffix])% {{{
		md=loadmodel(org,['Inversion_drag_ISMIP', suffix]);

		% load 1981 DEM
		md.timestepping.start_time = 1981;
		%md.timestepping.final_time = finalTime;
		md.timestepping.final_time = 2020;

		%filename = [projPath, './DATA/HGL_150m_mo_2sep23.xy'];
		filename = ['./DATA/HGL_150m_mo_2sep23.xy'];
		md.geometry.surface = loaddataFromxy(md, filename);

		% load 1981 front position
		disp('Using Greene monthly ice front--->')
		% get observed calving front
		mask = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [md.timestepping.start_time, md.timestepping.final_time]);
		% convert icemask to levelset distance
		distance = zeros(size(mask));
		for i = 1:size(mask,2)
			distance(1:end-1,i) = reinitializelevelset(md, mask(1:end-1,i));
		end
		distance(end,:) = mask(end,:);
		% transient spc
		md.levelset.spclevelset = distance;


		% adjust thickness accordingly



		savemodel(org,md);
	end%}}}
	if perform(org,'Transient_Calibration_withBeenDEM') % {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% reset the final time
		md.timestepping.final_time = finalTime;
		disp(['Run AD from ', num2str(md.timestepping.start_time), ' to ', num2str(md.timestepping.final_time)])

		%Prepare VEL obs =========================================================
		disp('Preparing observed velocities');
		obsData = load([projPath, './PostProcessing/Results/velObs_onmesh.mat']);

		timestamps = 0.5*(obsData.TStart+obsData.TEnd);
		pos = find((timestamps>=md.timestepping.start_time)&(timestamps<=md.timestepping.final_time));
		[~, ia, ic] = unique(timestamps(pos), 'stable');
		dataPos = pos(ia);
		disp(['  Found ', num2str(numel(dataPos)), ' observations between ', num2str(md.timestepping.start_time), ' and ', num2str(md.timestepping.final_time)])

		timestamps = timestamps(dataPos);
		vx_obs = obsData.vx_onmesh(:, dataPos);
		vy_obs = obsData.vy_onmesh(:, dataPos);

		%fill in NaNs
		mask_obs=(~isnan(vx_obs) & ~isnan(vy_obs));
		for i=1:numel(timestamps)
			flags = isnan(vx_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			vx_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), vx_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');

			flags = isnan(vy_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			vy_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), vy_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');
		end

		%Prepare cost functions
		count = 1;
		weights = ones(md.mesh.numberofvertices,2);
		for i = 1:length(timestamps)

			%Filter out observations that are not within our simulation period
			if timestamps(i)<md.timestepping.start_time ||  timestamps(i)>md.timestepping.final_time
				continue
			end

			% 0.5*(Vx-Vx_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['VxMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Vx', 'observation_string','VxObs', 'observation', vx_obs(:, i)./md.constants.yts,...
				'weights',weights(:,1).*mask_obs(:,i), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;

			% 0.5*(Vy-Vy_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['VyMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Vy', 'observation_string','VyObs', 'observation', vy_obs(:, i)./md.constants.yts,...
				'weights',weights(:,2).*mask_obs(:,i), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;
		end

		%Prepare surface obs =================================================================
		disp('Preparing Surface observations');


		% Benn's DEM
		timestamps = [2009.75];
		bennDEM = interpFromGeotiff([projPath, '/DATA/z0_2019.75_smb_firn_corr.tif'], md.mesh.x, md.mesh.y, NaN); % This is WGS84
		geoid = interpBedmachineGreenland(md.mesh.x, md.mesh.y,'geoid');
		z_obs = bennDEM - geoid;
		thickness = z_obs - md.geometry.bed;

		% not include NaNs, and thickness< 100m in the inversion
		mask_obs= ((~isnan(z_obs)) & thickness > 300);
		%fill in NaNs
		for i=1:numel(timestamps)
			flags = isnan(z_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			z_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), z_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');
		end

		%Prepare cost functions
		weights = 1/(1e5*md.constants.yts)*ones(md.mesh.numberofvertices,1);
		for i = 1:length(timestamps)
			%Filter out observations that are not within our simulation period
			if timestamps(i)<md.timestepping.start_time ||  timestamps(i)>md.timestepping.final_time
				continue
			end
			% 0.5*(S-S_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['DEMMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Surface', 'observation_string','SurfaceObservation', 'observation', z_obs(:, i),...
				'weights',weights.*mask_obs(:,i), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;
		end

		%Prepare Tikhonov ==============================================================
		disp('Preparing Tikhonov regularization');

		if isprop(md.friction, 'C')
			frictionCoeffName = 'FrictionC';
			tikhonovCoeff = 1e-12; % For Budd
		else
			frictionCoeffName = 'FrictionCoefficient';
			tikhonovCoeff = 1e-7; % For Budd
		end

		md.outputdefinition.definitions{count}=cfdragcoeffabsgrad('name','TikhonovReg',...
			'definitionstring',['Outputdefinition' num2str(count)],...
			'model_string', frictionCoeffName,...
			'weights',tikhonovCoeff*ones(md.mesh.numberofvertices,1),...
			'weights_string','WeightsSurfaceObservation');
		md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],...
			'type','scalar','fos_reverse_index',1);
		count = count +1;

		% Controls
		for i=1:1
			switch i
				case 1;
					if flagFriction < 2 % Schoof or Weertman
						field =md.friction.C; name = 'FrictionC'; scaling = 3000;
					else
						field =md.friction.coefficient; name = 'FrictionCoefficient'; scaling = 10;
					end
				case 2;
					field =md.materials.rheology_B;
					name = 'MaterialsRheologyBbar';
					scaling = 1e8;
				case 3;
					field =md.smb.mass_balance;
					name = 'SmbMassBalance';
					scaling = 1/md.constants.yts;
				otherwise error('not supported');
			end
			md.autodiff.independents{i} = independent('name',name,'type','vertex','nods',md.mesh.numberofvertices,...
				'control_size', size(field,2), 'min_parameters',1e-5*field, 'max_parameters',100*field, 'control_scaling_factor',scaling);
		end

		% General settings
		md.autodiff.driver='fos_reverse';
		md.autodiff.isautodiff=true;
		md.settings.checkpoint_frequency = 5;
		md.inversion=adm1qn3inversion(md.inversion);
		md.inversion.iscontrol=1;
		md.inversion.maxiter=30;
		md.inversion.maxsteps=md.inversion.maxiter;
		md.inversion.dxmin=1e-15;
		md.verbose.convergence=0;
		md.verbose=verbose('solution',true,'control',true);
		md.transient.requested_outputs={'default', 'IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset'};

		md.stressbalance.maxiter=20;

		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_floatingice = 200*10^3;
		md.calving.stress_threshold_groundedice = sigma*10^6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:

		% meltingrate
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		md.cluster = cluster;
		md.settings.waitonlock = 0; % not wait for output 
		md.toolkits.DefaultAnalysis=issmmumpssolver();
		md.toolkits.RecoveryAnalysis=issmmumpssolver();
		%md.settings.solver_residue_threshold = 1.e-2;
		md.verbose.solution = 0;
		if strcmp(md.cluster.name, 'discovery')
			md.cluster.memory=500;
			md.cluster.interactive = 0; 
		else
			md.cluster.interactive = 0;
		end
		md.miscellaneous.name = [savePath];

		md=solve(md,'Transient','runtimename',false);

		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org,'Transient_Calibration_debug') % {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% reset the final time
		md.timestepping.final_time = finalTime;

		%Prepare VEL obs =========================================================
		disp('Preparing observed velocities');
		obsData = load([projPath, './PostProcessing/Results/velObs_onmesh.mat']);

		timestamps = 0.5*(obsData.TStart+obsData.TEnd);
		pos = find((timestamps>=md.timestepping.start_time)&(timestamps<=md.timestepping.final_time));
		[~, ia, ic] = unique(timestamps(pos), 'stable');
		dataPos = pos(ia);

		disp(['  Found ', num2str(numel(dataPos)), ' observations between ', num2str(md.timestepping.start_time), ' and ', num2str(md.timestepping.final_time)])

		timestamps = timestamps(dataPos);
		vx_obs = obsData.vx_onmesh(:, dataPos);
		vy_obs = obsData.vy_onmesh(:, dataPos);

		%obsData = interpFromMEaSUREsGeotiff(md.mesh.x,md.mesh.y, timestamps(1), timestamps(2), 'glacier', 'Helheim');
		%vx_obs = [obsData(:).vx];
		%vy_obs = [obsData(:).vy];

		%fill in NaNs
		mask_obs=(~isnan(vx_obs) & ~isnan(vy_obs));
		for i=1:numel(timestamps)
			flags = isnan(vx_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			vx_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), vx_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');

			flags = isnan(vy_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			vy_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), vy_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');
		end

		%Prepare cost functions
		% add obs (initial vel) every 10 timesteps
		timestamps = linspace(md.timestepping.start_time, md.timestepping.final_time, ceil((md.timestepping.final_time-md.timestepping.start_time)/0.05) +1);

		count = 1;
		weights = ones(md.mesh.numberofvertices,2);
		for i = 1:length(timestamps)
			%Filter out observations that are not within our simulation period
			if timestamps(i)<md.timestepping.start_time ||  timestamps(i)>md.timestepping.final_time
				continue
			end
			% 0.5*(Vx-Vx_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['VxMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Vx', 'observation_string','VxObs', 'observation', vx_obs(:, 1)./md.constants.yts,...
				'weights',weights(:,1).*mask_obs(:,1), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;

			% 0.5*(Vy-Vy_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['VyMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Vy', 'observation_string','VyObs', 'observation', vy_obs(:, 1)./md.constants.yts,...
				'weights',weights(:,2).*mask_obs(:,1), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;
		end

		%Prepare surface obs =================================================================
		disp('Preparing Surface observations');

		%	timestamps = [2009.75];
		%		md.timestepping.final_time = 2019.75;
		timestamps = md.timestepping.final_time;
		% Benn's DEM
		bennDEM = interpFromGeotiff([projPath, '/DATA/z0_2019.75_smb_firn_corr.tif'], md.mesh.x, md.mesh.y, NaN); % This is WGS84
		geoid = interpBedmachineGreenland(md.mesh.x, md.mesh.y,'geoid');
		z_obs = bennDEM - geoid;
		thickness = z_obs - md.geometry.bed;

		% not include NaNs, and thickness< 100m in the inversion
		mask_obs= ((~isnan(z_obs)) & thickness > 300);
		%fill in NaNs
		for i=1:numel(timestamps)
			flags = isnan(z_obs(:,i));
			pos1  = find(flags); pos2 = find(~flags);
			z_obs(pos1,i) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), z_obs(pos2,i), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');
		end

		%Prepare cost functions
		weights = 1/(1e5*md.constants.yts)*ones(md.mesh.numberofvertices,1);
		for i = 1:length(timestamps)
			%Filter out observations that are not within our simulation period
			if timestamps(i)<md.timestepping.start_time ||  timestamps(i)>md.timestepping.final_time
				continue
			end
			% 0.5*(S-S_obs).^2
			md.outputdefinition.definitions{count} = cfsurfacesquare('name',['DEMMisfit' num2str(count)],...
				'definitionstring',['Outputdefinition' num2str(count)],...
				'model_string','Surface', 'observation_string','SurfaceObservation', 'observation', z_obs(:, i),...
				'weights',weights.*mask_obs(:,i), 'weights_string','WeightsSurfaceObservation', 'datatime', timestamps(i));
			md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
			count = count +1;
		end

		%Prepare Tikhonov ==============================================================
		disp('Preparing Tikhonov regularization');

		if isprop(md.friction, 'C')
			frictionCoeffName = 'FrictionC';
			tikhonovCoeff = 1e-12; % For Budd
		else
			frictionCoeffName = 'FrictionCoefficient';
			tikhonovCoeff = 1e-7; % For Budd
		end

		md.outputdefinition.definitions{count}=cfdragcoeffabsgrad('name','TikhonovReg',...
			'definitionstring',['Outputdefinition' num2str(count)],...
			'model_string', frictionCoeffName,...
			'weights',tikhonovCoeff*ones(md.mesh.numberofvertices,1),...
			'weights_string','WeightsSurfaceObservation');
		md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],...
			'type','scalar','fos_reverse_index',1);
		count = count +1;

		% Controls
		for i=1:1
			switch i
				case 1;
					if flagFriction < 2 % Schoof or Weertman
						field =md.friction.C; name = 'FrictionC'; scaling = 3000;
					else
						field =md.friction.coefficient; name = 'FrictionCoefficient'; scaling = 10;
					end
				case 2;
					field =md.materials.rheology_B;
					name = 'MaterialsRheologyBbar';
					scaling = 1e8;
				case 3;
					field =md.smb.mass_balance;
					name = 'SmbMassBalance';
					scaling = 1/md.constants.yts;
				otherwise error('not supported');
			end
			md.autodiff.independents{i} = independent('name',name,'type','vertex','nods',md.mesh.numberofvertices,...
				'control_size', size(field,2), 'min_parameters',1e-5*field, 'max_parameters',100*field, 'control_scaling_factor',scaling);
		end

		% General settings
		md.autodiff.driver='fos_reverse';
		md.autodiff.isautodiff=true;
		md.settings.checkpoint_frequency = 5;
		md.inversion=adm1qn3inversion(md.inversion);
		md.inversion.iscontrol=1;
		md.inversion.maxiter=1;
		md.inversion.maxsteps=md.inversion.maxiter;
		md.inversion.dxmin=1e-15;
		md.verbose.convergence=0;
		md.verbose=verbose('solution',true,'control',true);
		md.transient.requested_outputs={'default', 'IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset'};

		md.stressbalance.maxiter=20;

		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_floatingice = 200*10^3;
		md.calving.stress_threshold_groundedice = sigma*10^6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:

		% meltingrate
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		md.cluster = cluster;
		md.settings.waitonlock = 0; % not wait for output 
		md.toolkits.DefaultAnalysis=issmmumpssolver();
		md.toolkits.RecoveryAnalysis=issmmumpssolver();
		%md.settings.solver_residue_threshold = 1.e-2;
		md.verbose.solution = 0;
		md.debug.valgrind = 0;
		if strcmp(md.cluster.name, 'discovery')
			md.cluster.memory=500;
			md.cluster.interactive = 0; %only needed if you are using the generic cluster
		else
			md.cluster.interactive = 0; %only needed if you are using the generic cluster
		end
		md.miscellaneous.name = [savePath];

		md=solve(md,'Transient','runtimename',false);

		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org, ['ResetFrontFrom_ADTransient', calving_suffix, suffix])% {{{

		md=loadmodel(org, ['Param_ISMIP', damage_suffix]);

		% load transient AD results, just use the first step surface DEM for inversion
		switch flagFriction
			case 0 % Schoof
				ADFolder = '20230911_Helheim_transient_useADinv_Schoof/'; % spc H
			otherwise
				error('The friction in the setting is not used.');
		end
		ADPath = [projPath, '/Models/', ADFolder, 'Model_Helheim_Big_Transient'];
		disp(['loading transient solution from ', ADPath])
		md_AD = loadmodel(ADPath);

		% set ice mask
		md.mask.ice_levelset = sign(md_AD.results.TransientSolution(1).MaskIceLevelset);

		% set geometry
		disp(['      reading surface from first time step of ', ADFolder]);
		md.geometry.surface = md_AD.results.TransientSolution(1).Surface;
		pos = find(md.mask.ice_levelset>0);
		md.geometry.surface(pos) = md.geometry.base(pos)+10; %Minimum thickness

		md.geometry.thickness = md.geometry.surface - md.geometry.bed;
		pos=find(md.geometry.thickness<=10);
		md.geometry.surface(pos) = md.geometry.base(pos)+10; %Minimum thickness
		md.geometry.thickness = md.geometry.surface - md.geometry.bed;

		md.masstransport.min_thickness = 10;

		disp('   Adjusting ice mask');
		%Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
		pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0);
		md.mask.ice_levelset(md.mesh.elements(pos,:)) = 1;
		% For the region where surface is NaN, set thickness to small value (consistency requires >0)
		pos=find((md.mask.ice_levelset<0).*(md.geometry.surface<0));
		md.mask.ice_levelset(pos)=1;
		pos=find((md.mask.ice_levelset<0).*(isnan(md.geometry.surface)));
		md.mask.ice_levelset(pos)=1;

		disp('      -- reconstruct thickness');
		md.geometry.thickness=md.geometry.surface-md.geometry.base;

		% Set the friction law to schoof's
		md.friction=frictionschoof();
		md.friction.m = 1.0/3.0*ones(md.mesh.numberofelements,1);
		md.friction.Cmax = Cmax*ones(md.mesh.numberofvertices,1);
		md.friction.C = 2000*ones(md.mesh.numberofvertices,1);
		md.friction.coupling = 2;

		%No friction on PURELY ocean element
		pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
		flags=ones(md.mesh.numberofvertices,1);
		flags(md.mesh.elements(pos_e,:))=0;
		md.friction.C(find(flags))=0.01;

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
		md.inversion.maxsteps=500;
		md.inversion.maxiter =500;
		md.inversion.min_parameters=1e-4*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=1e5*ones(md.mesh.numberofvertices,1);
		md.inversion.control_scaling_factors=1;
		md.inversion.dxmin = 1e-6;
		%Additional parameters
		md.stressbalance.restol=1e-5;
		md.stressbalance.reltol=1e-3;
		md.stressbalance.abstol=NaN;
		md.stressbalance.maxiter = 100;

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		%Go solve
		md.cluster=cluster;
		md=solve(md,'sb');

		% get a direction to update cost coeff
		newCoeff = updateCostCoeff(md, ratio);
		disp(sprintf('With the given ratio %d, %d, %d \n', ratio))
		disp(sprintf('The coefficients can be updated to %g,   %g,   %g \n', newCoeff))

		%Put results back into the model
		md.friction.C=md.results.StressbalanceSolution.FrictionC;
		md.initialization.vx=md.results.StressbalanceSolution.Vx;
		md.initialization.vy=md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end%}}}
	if perform(org, ['PrepareForPINN_flowline_transient'])% {{{

		md=loadmodel('./Models/20231018_Helheim_ERA5_ISMIP_ResetFront_Schoof/Model_Helheim_Big_Transient.mat');

		mu = md.materials.rheology_B(1);

		disp('  Loading the central flowline')
		flowline = load([projPath, 'PostProcessing/Results/flowlines_', glacier, '_3.mat']);
		x2d = flowline.flowlineList{2}.x;
		y2d = flowline.flowlineList{2}.y;
		x = flowline.flowlineList{2}.Xmain(:) .* 1e3;
		disp(['    Found ', num2str(size(x)), ' points on the flowline']);

		% time index
		% 200 time steps over the first year
		timeid = [1:200];
		t = [md.results.TransientSolution(timeid).time].*md.constants.yts;

		% get data from the results
		vel_obs = sqrt([md.results.TransientSolution(:).Vx].^2+[md.results.TransientSolution(:).Vy].^2);
		h_obs = [md.results.TransientSolution(:).Surface];
		H_obs = [md.results.TransientSolution(:).Thickness];
		smb_obs = [md.results.TransientSolution(:).SmbMassBalance];
		icemask_obs = [md.results.TransientSolution(:).MaskIceLevelset];

		% project to x-t plane
		vel = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, vel_obs(:,timeid)./md.constants.yts, x2d, y2d);
		h = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, h_obs(:,timeid), x2d, y2d);
		H = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, H_obs(:,timeid), x2d, y2d);
		C = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, md.friction.C, x2d, y2d);
		smb = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, smb_obs(:,timeid)./md.constants.yts, x2d, y2d); % need to convert to m/s
		ice_levelset = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, icemask_obs(:,timeid), x2d, y2d);

		DBC = logical(zeros(length(x2d),length(timeid)));
		DBC(1,:) = true;
		DBC(end,:) = true;
		DBC = DBC & (ice_levelset<0);
		icemask = (ice_levelset<0) & (~DBC);

		% find the calving front position
		cx = interpZeroPos(x, ice_levelset);
		ct = t(:);
		nx = ones(size(cx, 1), 1);

		% generate colocation points
		Xmin = [min(x), min(t)];
		Xmax = [max(x), max(t)];
		Nf = 15000;
		X_ = Xmin + (Xmax - Xmin) .* lhsdesign(Nf, 2);
		icemask_ = InterpFromGridToMesh(x(:), t(:), icemask'+0, X_(:,1), X_(:,2),0);
		X_f = X_(icemask_>0.8, :);
		icemask_f = icemask_(icemask_>0.8);

		% inversion for flowline as a 2d problem {{{
		x1d = x;
		Lx = max(x1d) - min(x1d);
		Ly = 4e4;
		Nx = length(x1d);
		ny = 50;
		yts = md.constants.yts;
		y1d = linspace(0, Ly, ny)';
		vx1d = repmat(vel(:,1)', ny, 1).*yts;
		vy1d = zeros(size(vx1d));
		h1d = repmat(h(:,1)', ny, 1);
		H1d = repmat(H(:,1)', ny, 1);
		C1d = repmat(C', ny, 1);
		bed = h1d - H1d;

		disp('   Creating 2D mesh');
		md = squaremesh(model(), Lx, Ly, Nx, ny);
		% project data
		disp('   Projecting flowline data to 2D mesh');
		md.initialization.vx = InterpFromGrid(x1d, y1d, vx1d, md.mesh.x, md.mesh.y);
		md.initialization.vy = InterpFromGrid(x1d, y1d, vy1d, md.mesh.x, md.mesh.y);
		md.geometry.bed = InterpFromGrid(x1d, y1d, bed, md.mesh.x, md.mesh.y);
		md.geometry.base = InterpFromGrid(x1d, y1d, bed, md.mesh.x, md.mesh.y);
		md.geometry.surface = InterpFromGrid(x1d, y1d, h1d, md.mesh.x, md.mesh.y);

		% set masks
		disp('   Setting ice mask');
		icemask1d = -1*ones(size(md.mesh.x));
		icemask1d(md.mesh.x>cx(1)) = 1;
		md.mask.ice_levelset = icemask1d;
		pos = find(md.mask.ice_levelset>0);
		md.geometry.surface(pos) = md.geometry.base(pos)+10; %Minimum thickness
		md.geometry.thickness = md.geometry.surface - md.geometry.bed;
		pos=find(md.geometry.thickness<=10);
		md.geometry.surface(pos) = md.geometry.base(pos)+10; %Minimum thickness
		md.geometry.thickness = md.geometry.surface - md.geometry.bed;
		md.masstransport.min_thickness = 10;

		disp('   Adjusting ice mask');
		%Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
		pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0);
		md.mask.ice_levelset(md.mesh.elements(pos,:)) = 1;
		% For the region where surface is NaN, set thickness to small value (consistency requires >0)
		pos=find((md.mask.ice_levelset<0).*(md.geometry.surface<0));
		md.mask.ice_levelset(pos)=1;
		pos=find((md.mask.ice_levelset<0).*(isnan(md.geometry.surface)));
		md.mask.ice_levelset(pos)=1;

		disp('      -- reconstruct thickness');
		md.geometry.thickness=md.geometry.surface-md.geometry.base;

		disp('      reading velocities ');
		md.inversion.vx_obs = md.initialization.vx;
		md.inversion.vy_obs = md.initialization.vy;
		md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
		md.initialization.vz  = zeros(md.mesh.numberofvertices,1);
		md.initialization.vel = md.inversion.vel_obs;

		%		disp('   Initialize basal friction using driving stress');
		%		disp('   -- Smooth the ice surface with 20 L2 projections and then compute the surface slopes');
		%		asurf    = averaging(md,md.geometry.surface,20); % maybe executing 20 L2 projection is ok
		%		[sx,sy,s]= slope(md,asurf); % slope 's' comes on elements
		%		sslope   = averaging(md,s,1); % average the slope once on the vertices, because 's' comes on elements, we need this data on vertices

		% set rheology
		disp('   Creating flow law parameters (assume ice is at -5°C for now)');
		md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);
		md.materials.rheology_B = mu .*ones(md.mesh.numberofvertices,1);

		%Deal with boundary conditions:
		disp('   Set Boundary conditions');
		md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
		md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
		md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
		md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
		md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);
		pos=find((md.mask.ice_levelset<0).*(md.mesh.vertexonboundary));
		md.stressbalance.spcvx(pos)=md.initialization.vx(pos);
		md.stressbalance.spcvy(pos)=md.initialization.vy(pos);
		md.stressbalance.spcvz(pos)=0;

		disp('   Initial basal friction ');
		md.friction = frictionweertman();
		md.friction.m = 3.0*ones(md.mesh.numberofelements,1);
		md.friction.C = 2000*ones(md.mesh.numberofvertices,1);

		%No friction on PURELY ocean element
		pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
		flags=ones(md.mesh.numberofvertices,1);
		flags(md.mesh.elements(pos_e,:))=0;
		md.friction.C(find(flags))=0.0;

		md=setflowequation(md,'SSA','all');

		md.mask.ocean_levelset = ones(md.mesh.numberofvertices,1);
		md.cluster=generic('name',oshostname(),'np',40);
		md.miscellaneous.name = ['test'];

		%Control general
		md.inversion=m1qn3inversion(md.inversion);
		md.inversion.iscontrol=1;
		md.verbose=verbose('solution',false,'control',true);
		md.transient.amr_frequency = 0;

		%Cost functions
		md.inversion.cost_functions=[101 103 501];
		md.inversion.cost_functions_coefficients=zeros(md.mesh.numberofvertices,numel(md.inversion.cost_functions));
		md.inversion.cost_functions_coefficients(:,1)=1000;
		md.inversion.cost_functions_coefficients(:,2)=0.7;
		md.inversion.cost_functions_coefficients(:,3)=3.525e-8;
		pos=find(md.mask.ice_levelset>0);
		md.inversion.cost_functions_coefficients(pos,1:2)=0;

		%Controls
		md.inversion.control_parameters={'FrictionC'};
		md.inversion.maxsteps=500;
		md.inversion.maxiter =500;
		md.inversion.min_parameters=0.01*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=5e4*ones(md.mesh.numberofvertices,1);
		md.inversion.control_scaling_factors=1;
		md.inversion.dxmin = 1e-6;
		%Additional parameters
		md.stressbalance.restol=0.01;
		md.stressbalance.reltol=0.1;
		md.stressbalance.abstol=NaN;

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md=solve(md,'Stressbalance');

		ratio = [10,5,1];
		% get a direction to update cost coeff
		newCoeff = updateCostCoeff(md, ratio);
		disp(sprintf('With the given ratio %d, %d, %d \n', ratio))
		disp(sprintf('The coefficients can be updated to %g,   %g,   %g \n', newCoeff))

		%Put results back into the model
		md.friction.C=md.results.StressbalanceSolution.FrictionC;
		md.initialization.vx=md.results.StressbalanceSolution.Vx;
		md.initialization.vy=md.results.StressbalanceSolution.Vy;

		% project C back to flowline
		y = Ly/2*ones(size(x));
		C = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, md.friction.C, x1d, y);
		%}}}


		% create x-t grid 
		[T, X] = meshgrid(t, x); % x-column, t-row

		% only B.C.
		X_bc = [X(DBC), T(DBC)];
		u_bc = [vel(DBC), h(DBC), H(DBC), smb(DBC)];

		% TODO: add C at the boundary

		% reshape to 1D array before saving
		x = X(:); t = T(:);
		vel = vel(:);
		h = h(:);
		H = H(:);
		smb = smb(:);
		DBC = DBC(:);
		icemask = icemask(:);

		% save
		savefile = [projPath, '/DATA/', glacier, suffix,'_PINN_flowline_transient_recomputeC.mat'];
		disp(['  Saving to ', savefile]);
		save(savefile, ...
			'x', 't', 'vel', 'h', 'H', 'smb', 'mu', 'DBC', 'icemask', ...      % obs data
			'x1d', 'C', ... % for validation 
			'cx', 'ct', 'nx', ...
			'X_bc', 'u_bc',...
			'X_f', 'icemask_f');                                                 % collocation points
	end%}}}

	% step 21--25
	if perform(org, ['PrepareForPINN_2DTransient'])% {{{

		% make a ref model at the fastflow region
		mdRef = triangle(model,['Exp/fastflow_CF.exp'],200);
		x = mdRef.mesh.x;
		y = mdRef.mesh.y;

		% load the model results from a transient simulation
		md=loadmodel('./Models/20231018_Helheim_ERA5_ISMIP_ResetFront_Schoof/Model_Helheim_Big_Transient.mat');

		% visocsity
		mu = md.materials.rheology_B(1);

		% time index
		% 200 time steps over the first year
		timeid = [1:200];
		t = [md.results.TransientSolution(timeid).time].*md.constants.yts;

		% get data from the results
		vel_obs = sqrt([md.results.TransientSolution(:).Vx].^2+[md.results.TransientSolution(:).Vy].^2);
		h_obs = [md.results.TransientSolution(:).Surface];
		H_obs = [md.results.TransientSolution(:).Thickness];
		smb_obs = [md.results.TransientSolution(:).SmbMassBalance];
		icemask_obs = [md.results.TransientSolution(:).MaskIceLevelset];

		% project to x-t plane
		vel = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, vel_obs(:,timeid)./md.constants.yts, x, y);
		h = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, h_obs(:,timeid), x, y);
		H = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, H_obs(:,timeid), x, y);
		C = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, md.friction.C, x, y);
		smb = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, smb_obs(:,timeid)./md.constants.yts, x, y); % need to convert to m/s
		ice_levelset = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, icemask_obs(:,timeid), x, y);

		DBC = logical(zeros(length(x2d),length(timeid)));
		DBC(1,:) = true;
		DBC(end,:) = true;
		DBC = DBC & (ice_levelset<0);
		icemask = (ice_levelset<0) & (~DBC);

		% find the calving front position
		cx = interpZeroPos(x, ice_levelset);
		ct = t(:);
		nx = ones(size(cx, 1), 1);

		% generate colocation points
		Xmin = [min(x), min(y), min(t)];
		Xmax = [max(x), min(y), max(t)];
		Nf = 15000;
		X_ = Xmin + (Xmax - Xmin) .* lhsdesign(Nf, 2);
		icemask_ = InterpFromGridToMesh(x(:), t(:), icemask'+0, X_(:,1), X_(:,2),0);
		X_f = X_(icemask_>0.8, :);
		icemask_f = icemask_(icemask_>0.8);

		% create x-y-t grid
		[T, X] = meshgrid(t, x); % x-column, t-row

		% only B.C.
		X_bc = [X(DBC), T(DBC)];
		u_bc = [vel(DBC), h(DBC), H(DBC), smb(DBC)];

		% TODO: add C at the boundary

		% reshape to 1D array before saving
		x = X(:); t = T(:);
		vel = vel(:);
		h = h(:);
		H = H(:);
		smb = smb(:);
		DBC = DBC(:);
		icemask = icemask(:);

		% save
		savefile = [projPath, '/DATA/', glacier, suffix,'_PINN_transient2D.mat'];
		disp(['  Saving to ', savefile]);
		save(savefile, ...
			'x', 't', 'vel', 'h', 'H', 'smb', 'mu', 'DBC', 'icemask', ...      % obs data
			'x1d', 'C', ... % for validation
			'cx', 'ct', 'nx', ...
			'X_bc', 'u_bc',...
			'X_f', 'icemask_f');                                                 % collocation points
	end%}}}
	if perform(org, ['Inversion_drag_forPINN_fastflowCF'])% {{{

		md_ref=loadmodel(org, ['Inversion_drag_ISMIP', damage_suffix, '_Weertman']);

		md=triangle(model,['Exp/fastflow_CF.exp'],200);
		x = md.mesh.x;
		y = md.mesh.y;

		% laod vel obs from 1-ITSLIVE 00, 2-ITSLIVE 2007, 3-Rignot2012, 3-Joughin
		velSource = 4;
		if velSource == 1
			[vx, vy] = interpFromITSLIVE(x, y, 0, 0);
			vx = vx(1:end-1, 1);
			vy = vy(1:end-1, 1);
			velName = 'ITSLIVE0';
		elseif velSource == 2
			[vx, vy] = interpFromITSLIVE(x, y, 2007, 2007);
			vx = vx(1:end-1, 1);
			vy = vy(1:end-1, 1);
			velName = 'ITSLIVE2007';
		elseif velSource == 3
			[vx, vy] = interpRignot2012(x, y);
			velName = 'Rignot2012';
		elseif velSource == 4
			[vx, vy] = interpJoughinCompositeGreenland(x, y);
			velName = 'JoughinComposite';
		else
			error('Unknown source of velocity observation')
		end

		%Important: scale velocity
		vx = vx ./ md.constants.yts;
		vy = vy ./ md.constants.yts;

		% load from 1-Bedmachine, TODO:2-ICESat2
		demSource = 1;
		if demSource == 1
			disp('   Interpolating mask');
			mask = int8(interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask'));
			md.mask.ice_levelset= -1*ones(md.mesh.numberofvertices,1);
			pos = find(mask<1);
			md.mask.ice_levelset(pos)=1;

			disp('      reading MC bed (assumes no floating ice)');
			bed = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed');

			disp('      reading Howat surface');
			h = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'surface');
			pos = find(md.mask.ice_levelset>0);
			h(pos) = bed(pos)+10; %Minimum thickness

			H = h - bed;
			pos=find(H<=10);
			h(pos) = bed(pos)+10; %Minimum thickness
			H = h - bed;

			disp('   Adjusting ice mask');
			%Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
			pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0);
			md.mask.ice_levelset(md.mesh.elements(pos,:)) = 1;
			% For the region where surface is NaN, set thickness to small value (consistency requires >0)
			pos=find((md.mask.ice_levelset<0).*(h<0));
			md.mask.ice_levelset(pos)=1;
			pos=find((md.mask.ice_levelset<0).*(isnan(h)));
			md.mask.ice_levelset(pos)=1;
		else
			error('Unknown source of the DEM')
		end

		% not used anymore for frictionNN, just have this here to not change the python code
		C = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.friction.C, x, y);
		mu = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.materials.rheology_B, x, y);

		DBC = (md.mesh.vertexonboundary>0) & (md.mask.ice_levelset<0);
		icemask = (md.mask.ice_levelset<0) & (~DBC);

		Xmin = min([x, y]);
		Xmax = max([x, y]);
		Nf = 15000;
		X_ = Xmin + (Xmax - Xmin) .* lhsdesign(Nf, 2);
		icemask_ = InterpFromMeshToMesh2d(md.mesh.elements, x, y, icemask+0, X_(:,1), X_(:,2));
		X_f = X_(icemask_>0.8, :);
		icemask_f = icemask_(icemask_>0.8);

		% calving front boundary
		contourline = isoline(md, md.mask.ice_levelset, 'value', 0);
		tcx = contourline.x;
		tcy = contourline.y;
		% find the unique points only
		utc = unique([tcx, tcy], 'rows', 'stable');
		tcx = utc(:,1);
		tcy = utc(:,2);

		cx = 0.5 * (tcx(1:end-1) + tcx(2:end));
		cy = 0.5 * (tcy(1:end-1) + tcy(2:end));
		nx = (tcy(1:end-1) - tcy(2:end));
		ny = -(tcx(1:end-1) - tcx(2:end));
		nn = sqrt(nx.^2+ny.^2);
		nx = nx ./ nn;
		ny = ny ./ nn;

		% apply a filter
		windowSize = 10;
		wd = (1/windowSize)*ones(1,windowSize);
		smoothnx = filter(wd,1,nx,[],1);
		smoothny = filter(wd,1,ny,[],1);
		% need to normalize again
		nn = sqrt(smoothnx.^2+smoothny.^2);
		smoothnx = smoothnx ./ nn;
		smoothny = smoothny ./ nn;

		% set mask for only covered area
		icemask = (md.mask.ice_levelset<0) & (~DBC);

		% B.C.
		X_bc = [x(DBC), y(DBC)];
		u_bc = [vx(DBC), vy(DBC), h(DBC), H(DBC), C(DBC), mu(DBC)];

		% check the collocation points
		figure
		plot(X_f(:,1), X_f(:,2), 'o')
		hold on
		plot(x(DBC), y(DBC), 'ko')
		plot(cx, cy, 'r*')


		plotmodel(md, 'data', vx, 'figure', 2)
		hold on
		for i = 1:length(cx)
			plot([cx(i), cx(i)+smoothnx(i)*1000],[cy(i), cy(i)+smoothny(i)*1000], 'r')
			plot([cx(i), cx(i)+nx(i)*1000],[cy(i), cy(i)+ny(i)*1000],'b--')
		end
		legend('off')

		% remove nans in the data
		nanflag = isnan(sum(u_bc, 2));
		disp(['Found ', num2str(sum(nanflag)), ' NaNs in the boundary condition'])
		X_bc = X_bc(nanflag<1, :);
		u_bc = u_bc(nanflag<1, :);

		nanflag = isnan(vx+vy);
		disp(['Found ', num2str(sum(nanflag)), ' NaNs in the data'])
		x = x(nanflag<1);
		y = y(nanflag<1);
		h = h(nanflag<1);
		H = H(nanflag<1);
		C = C(nanflag<1);
		mu = mu(nanflag<1);
		DBC = DBC(nanflag<1);
		icemask = icemask(nanflag<1);
		vx = vx(nanflag<1);
		vy = vy(nanflag<1);
		% save
		if damageType == 0
			mu = md_ref.materials.rheology_B(1);
		end
		save([projPath, '/DATA/', glacier, '_PINN_fastflow', damage_suffix, '_', velName, '.mat'], ...
			'x', 'y', 'vx', 'vy', 'C', 'h', 'H', 'mu', 'DBC', 'icemask', ...		% obs data
			'cx', 'cy', 'nx', 'ny', 'smoothnx', 'smoothny', ...						% calving front
			'X_bc', 'u_bc',...
			'X_f', 'icemask_f');																	% collocation points
	end%}}}
	if perform(org, ['Inversion_drag_forPINN'])% {{{

		md_ref=loadmodel(org, ['Inversion_drag_ISMIP', damage_suffix, '_Weertman']);
		mu = md_ref.materials.rheology_B(1);

		md=triangle(model,['Exp/Helheim_Big.exp'], 250);

		x = md.mesh.x;
		y = md.mesh.y;

		% laod vel obs from 1-ITSLIVE 00, 2-ITSLIVE 2007, 3-Rignot2012, 3-Joughin
		if velSource == 1
			[vx, vy] = interpFromITSLIVE(x, y, 0, 0);
			vx = vx(1:end-1, 1);
			vy = vy(1:end-1, 1);
			velName = 'ITSLIVE0';
		elseif velSource == 2
			[vx, vy] = interpFromITSLIVE(x, y, 2007, 2007);
			vx = vx(1:end-1, 1);
			vy = vy(1:end-1, 1);
			velName = 'ITSLIVE2007';
		elseif velSource == 3
			[vx, vy] = interpRignot2012(x, y);
			velName = 'Rignot2012';
		elseif velSource == 4
			[vx, vy] = interpJoughinCompositeGreenland(x, y);
			velName = 'JoughinComposite';
		else
			error('Unknown source of velocity observation')
		end

		%Important: scale velocity
		vx = vx ./ md.constants.yts;
		vy = vy ./ md.constants.yts;

		% geometry
		disp('   Interpolating mask');
		mask = int8(interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask'));
		md.mask.ice_levelset= -1*ones(md.mesh.numberofvertices,1);
		pos = find(mask<1);
		md.mask.ice_levelset(pos)=1;

		disp('      reading MC bed (assumes no floating ice)');
		bed = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed');

		disp('      reading Howat surface');
		h = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'surface');
		pos = find(md.mask.ice_levelset>0);
		h(pos) = bed(pos)+10; %Minimum thickness

		H = h - bed;
		pos=find(H<=10);
		h(pos) = bed(pos)+10; %Minimum thickness
		H = h - bed;

		disp('   Adjusting ice mask');
		%Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
		pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0);
		md.mask.ice_levelset(md.mesh.elements(pos,:)) = 1;
		% For the region where surface is NaN, set thickness to small value (consistency requires >0)
		pos=find((md.mask.ice_levelset<0).*(h<0));
		md.mask.ice_levelset(pos)=1;
		pos=find((md.mask.ice_levelset<0).*(isnan(h)));
		md.mask.ice_levelset(pos)=1;

		% not used anymore for frictionNN, just have this here to not change the python code
		C = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.friction.C, x, y);

		DBC = (md.mesh.vertexonboundary>0) & (md.mask.ice_levelset<0);
		icemask = (md.mask.ice_levelset<0) & (~DBC);

		Xmin = min([x, y]);
		Xmax = max([x, y]);
		Nf = floor(1.5*md.mesh.numberofvertices);
		X_ = Xmin + (Xmax - Xmin) .* lhsdesign(Nf, 2);
		icemask_ = InterpFromMeshToMesh2d(md.mesh.elements, x, y, icemask+0, X_(:,1), X_(:,2));
		X_f = X_(icemask_>0.8, :);
		icemask_f = icemask_(icemask_>0.8);

		% calving front boundary
		contourline = isoline(md, md.mask.ice_levelset, 'value', 0);
		tcx = contourline.x;
		tcy = contourline.y;
		% find the unique points only
		utc = unique([tcx, tcy], 'rows', 'stable');
		tcx = utc(:,1);
		tcy = utc(:,2);

		cx = 0.5 * (tcx(1:end-1) + tcx(2:end));
		cy = 0.5 * (tcy(1:end-1) + tcy(2:end));
		nx = (tcy(1:end-1) - tcy(2:end));
		ny = -(tcx(1:end-1) - tcx(2:end));
		nn = sqrt(nx.^2+ny.^2);
		nx = nx ./ nn;
		ny = ny ./ nn;

		% apply a filter
		windowSize = 10;
		wd = (1/windowSize)*ones(1,windowSize);
		smoothnx = filter(wd,1,nx,[],1);
		smoothny = filter(wd,1,ny,[],1);
		% need to normalize again
		nn = sqrt(smoothnx.^2+smoothny.^2);
		smoothnx = smoothnx ./ nn;
		smoothny = smoothny ./ nn;

		% set mask for only covered area
		icemask = (md.mask.ice_levelset<0) & (~DBC);

		% B.C.
		X_bc = [x(DBC), y(DBC)];
		u_bc = [vx(DBC), vy(DBC), h(DBC), H(DBC), C(DBC)];

		% check the collocation points
		figure
		plot(X_f(:,1), X_f(:,2), 'o')
		hold on
		plot(x(DBC), y(DBC), 'ko')
		plot(cx, cy, 'r*')


		plotmodel(md, 'data', vx, 'figure', 2)
		hold on
		for i = 1:length(cx)
			plot([cx(i), cx(i)+smoothnx(i)*1000],[cy(i), cy(i)+smoothny(i)*1000], 'r')
			plot([cx(i), cx(i)+nx(i)*1000],[cy(i), cy(i)+ny(i)*1000],'b--')
		end
		legend('off')

		% remove nans in the data
		nanflag = isnan(sum(u_bc, 2));
		disp(['Found ', num2str(sum(nanflag)), ' NaNs in the boundary condition'])
		X_bc = X_bc(nanflag<1, :);
		u_bc = u_bc(nanflag<1, :);

		nanflag = isnan(vx+vy);
		disp(['Found ', num2str(sum(nanflag)), ' NaNs in the data'])
		x = x(nanflag<1);
		y = y(nanflag<1);
		h = h(nanflag<1);
		H = H(nanflag<1);
		C = C(nanflag<1);
		DBC = DBC(nanflag<1);
		icemask = icemask(nanflag<1);
		vx = vx(nanflag<1);
		vy = vy(nanflag<1);
		% save
		save([projPath, '/DATA/', glacier, '_PINN_obs_', velName, '.mat'], ...
			'x', 'y', 'vx', 'vy', 'C', 'h', 'H', 'mu', 'DBC', 'icemask', ...		% obs data
			'cx', 'cy', 'nx', 'ny', 'smoothnx', 'smoothny', ...						% calving front
			'X_bc', 'u_bc',...
			'X_f', 'icemask_f');																	% collocation points
	end%}}}
	if perform(org, ['Inversion_drag_forPINN_anisotropicMesh'])% {{{

		md=loadmodel(org, ['Inversion_drag_ISMIP', damage_suffix, '_Weertman']);
		mu = md.materials.rheology_B(1);

		x = md.mesh.x;
		y = md.mesh.y;

		% laod vel obs from 1-ITSLIVE 00, 2-ITSLIVE 2007, 3-Rignot2012, 3-Joughin
		if velSource == 1
			[vx, vy] = interpFromITSLIVE(x, y, 0, 0);
			vx = vx(1:end-1, 1);
			vy = vy(1:end-1, 1);
			velName = 'ITSLIVE0';
		elseif velSource == 2
			[vx, vy] = interpFromITSLIVE(x, y, 2007, 2007);
			vx = vx(1:end-1, 1);
			vy = vy(1:end-1, 1);
			velName = 'ITSLIVE2007';
		elseif velSource == 3
			[vx, vy] = interpRignot2012(x, y);
			velName = 'Rignot2012';
		elseif velSource == 4
			[vx, vy] = interpJoughinCompositeGreenland(x, y);
			velName = 'JoughinComposite';
		else
			error('Unknown source of velocity observation')
		end

		%Important: scale velocity
		vx = vx ./ md.constants.yts;
		vy = vy ./ md.constants.yts;

		% geometry
		H = md.geometry.thickness;
		h = md.geometry.surface;

		% not used anymore for frictionNN, just have this here to not change the python code
		C = md.friction.C;

		DBC = (md.mesh.vertexonboundary>0) & (md.mask.ice_levelset<0);
		icemask = (md.mask.ice_levelset<0) & (~DBC);

		% create 3 random points for each element
		[rx,ry,r]= getCenterPoints(md.mesh.elements, md.mesh.x, md.mesh.y);
		repeat = 3;
		Nf = md.mesh.numberofelements;
		% 1,3,5 scaling of r, 2,4,6 angles
		xrand = lhsdesign(Nf*repeat, 2);
		R = repmat(r, repeat, 1);
		x0 = repmat(rx, repeat, 1);
		y0 = repmat(ry, repeat, 1);
		X_ = [x0+R.*xrand(:,1).*cos(2*pi*xrand(:,2)), y0+R.*xrand(:,1).*sin(2*pi*xrand(:,2))];

		icemask_ = InterpFromMeshToMesh2d(md.mesh.elements, x, y, icemask+0, X_(:,1), X_(:,2));
		X_f = X_(icemask_>0.8, :);
		icemask_f = icemask_(icemask_>0.8);

		% calving front boundary
		contourline = isoline(md, md.mask.ice_levelset, 'value', 0);
		tcx = contourline.x;
		tcy = contourline.y;
		% find the unique points only
		utc = unique([tcx, tcy], 'rows', 'stable');
		tcx = utc(:,1);
		tcy = utc(:,2);

		cx = 0.5 * (tcx(1:end-1) + tcx(2:end));
		cy = 0.5 * (tcy(1:end-1) + tcy(2:end));
		nx = (tcy(1:end-1) - tcy(2:end));
		ny = -(tcx(1:end-1) - tcx(2:end));
		nn = sqrt(nx.^2+ny.^2);
		nx = nx ./ nn;
		ny = ny ./ nn;

		% apply a filter
		windowSize = 10;
		wd = (1/windowSize)*ones(1,windowSize);
		smoothnx = filter(wd,1,nx,[],1);
		smoothny = filter(wd,1,ny,[],1);
		% need to normalize again
		nn = sqrt(smoothnx.^2+smoothny.^2);
		smoothnx = smoothnx ./ nn;
		smoothny = smoothny ./ nn;

		% set mask for only covered area
		icemask = (md.mask.ice_levelset<0) & (~DBC);

		% B.C.
		X_bc = [x(DBC), y(DBC)];
		u_bc = [vx(DBC), vy(DBC), h(DBC), H(DBC), C(DBC)];

		% check the collocation points
		figure
		plot(X_f(:,1), X_f(:,2), 'o')
		hold on
		plot(x(DBC), y(DBC), 'ko')
		plot(cx, cy, 'r*')


		plotmodel(md, 'data', vx, 'figure', 2)
		hold on
		for i = 1:length(cx)
			plot([cx(i), cx(i)+smoothnx(i)*1000],[cy(i), cy(i)+smoothny(i)*1000], 'r')
			plot([cx(i), cx(i)+nx(i)*1000],[cy(i), cy(i)+ny(i)*1000],'b--')
		end
		legend('off')

		% remove nans in the data
		nanflag = isnan(sum(u_bc, 2));
		disp(['Found ', num2str(sum(nanflag)), ' NaNs in the boundary condition'])
		X_bc = X_bc(nanflag<1, :);
		u_bc = u_bc(nanflag<1, :);

		nanflag = isnan(vx+vy);
		disp(['Found ', num2str(sum(nanflag)), ' NaNs in the data'])
		x = x(nanflag<1);
		y = y(nanflag<1);
		h = h(nanflag<1);
		H = H(nanflag<1);
		C = C(nanflag<1);
		DBC = DBC(nanflag<1);
		icemask = icemask(nanflag<1);
		vx = vx(nanflag<1);
		vy = vy(nanflag<1);
		% save
		save([projPath, '/DATA/', glacier, '_PINN_obs_', velName, '.mat'], ...
			'x', 'y', 'vx', 'vy', 'C', 'h', 'H', 'mu', 'DBC', 'icemask', ...		% obs data
			'cx', 'cy', 'nx', 'ny', 'smoothnx', 'smoothny', ...						% calving front
			'X_bc', 'u_bc',...
			'X_f', 'icemask_f');																	% collocation points
	end%}}}
	if perform(org, ['prepare_subdomain'])% {{{

		md_ref=loadmodel(org, ['Inversion_drag_ISMIP', damage_suffix, '_Weertman']);

		% get the range of the larger domain
		xmin = min(md_ref.mesh.x); xmax = max(md_ref.mesh.x);
		ymin = min(md_ref.mesh.y); ymax = max(md_ref.mesh.y);

		if 0 
			row = ceil(subregion(3)/subregion(2));
			column = subregion(3) - (row-1)*subregion(2);
			dx = (xmax - xmin)/subregion(2);
			dy = (ymax - ymin)/subregion(1);
			%topleft = [xmin+(column-1)*dx, ymax-(row-1)*dy];
			%bottomleft = [xmin+(column-1)*dx, ymax-(row)*dy];
			%bottomright = [xmin+(column)*dx, ymax-(row)*dy];
			%topright = [xmin+(column)*dx, ymax-(row-1)*dy];
			
			% write to exp
			out = struct();
			out.density = 1;
			out.closed = 1;
			out.name = ['subdomain_',num2str(subregion(1)), '_',num2str(subregion(2)), '_', num2str(subregion(3)),'.exp'];
			outFile = [projPath, '/Exp/', out.name];
			out.x = [xmin+(column-1)*dx, xmin+(column-1)*dx, xmin+(column)*dx, xmin+(column)*dx ];
			out.y = [ymax-(row-1)*dy, ymax-(row)*dy, ymax-(row)*dy, ymax-(row-1)*dy];
			% close the exp
			out.x(end+1) = out.x(1);
			out.y(end+1) = out.y(1);

			out.nods = numel(out.x);

			disp(['  Save to exp file: ', outFile])
			expwrite(out, outFile);
		else
			subdomain = 'fastflow_CF';
		end
		%subdomain = 'fast_moving_large';
		md=triangle(model,['Exp/', subdomain, '.exp'],500);

		velx = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.inversion.vx_obs, md.mesh.x, md.mesh.y);
		vely = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.inversion.vy_obs, md.mesh.x, md.mesh.y);
		vel  = sqrt(velx.^2+vely.^2);
		md=bamg(md,'hmin',100,'hmax', 500,'field',vel,'err',5);

		x = md.mesh.x;
		y = md.mesh.y;
		md.materials.rheology_B = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.materials.rheology_B, x, y);
		md.inversion.vx_obs = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.inversion.vx_obs, x, y);
		md.inversion.vy_obs = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.inversion.vy_obs, x, y);
		md.initialization.vx = md.inversion.vx_obs;
		md.initialization.vy = md.inversion.vy_obs;
		md.geometry.surface = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.geometry.surface, x, y);
		md.geometry.thickness = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.geometry.thickness, x, y);
		md.friction = frictionweertman();
		md.friction.C = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, md_ref.friction.C, x, y);
		md.mask.ice_levelset = InterpFromMeshToMesh2d(md_ref.mesh.elements, md_ref.mesh.x, md_ref.mesh.y, (md_ref.mask.ice_levelset), x, y);


		% SMB
		md.smb.mass_balance =interpRACMO1km(md.mesh.x,md.mesh.y);
		md.balancethickness.thickening_rate =  interpSmith2020(md.mesh.x,md.mesh.y,'gris_filt');
		pos = find(~isnan(md.balancethickness.thickening_rate));
		tr = griddata(md.mesh.x(pos), md.mesh.y(pos), md.balancethickness.thickening_rate(pos), md.mesh.x, md.mesh.y,'natural');

		pos = find(~isnan(tr));
		F = scatteredInterpolant(md.mesh.x(pos),md.mesh.y(pos),tr(pos),'linear','linear');
		fillnan = F(md.mesh.x, md.mesh.y);

		md.balancethickness.thickening_rate = fillnan;
		md.balancethickness.thickening_rate(pos) = tr(pos);

		% save
		savemodel(org,md);
		% save as a struct
		saveasstruct(md, ['Helheim_', subdomain, '.mat'])
	end%}}}

	% step 26--30
	if perform(org, ['Prepare_for_GNSS'])% {{{

		md=loadmodel(org,'Mesh');

		% get the mesh
		elements=md.mesh.elements;
		x=md.mesh.x;
		y=md.mesh.y;

		%compute areas;
		eleAreas=GetAreas(elements,x,y);
		disp(['Minimum element size: ', num2str(min(eleAreas))])
		disp(['Maximum element size: ', num2str(max(eleAreas))])

		% get the center
		cx = 1/3*(x(elements(:,1)) + x(elements(:,2)) + x(elements(:,3)));
		cy = 1/3*(y(elements(:,1)) + y(elements(:,2)) + y(elements(:,3)));

		% convert to ll
		[lat, lon] = xy2ll(cx, cy, 1);

		% write to a txt file
		filename = [projPath, '/DATA/Helheim.txt'];

		data = [lat, lon, eleAreas];
		% Write matrix to text file
		writematrix(data, filename, 'delimiter', '\t');

		disp(['Matrix has been written to ' filename]);

		% add damage for shear margin
		savemodel(org,md);
	end%}}}
	if perform(org, ['Transient_Prep_RACMO_FromAD_', calving_suffix, suffix])% {{{

		md=loadmodel(org, ['Transient_ERA5_ISMIP_prep', calving_suffix, suffix]);

		% load AD results

		switch flagFriction
			case 0 % Schoof
				ADFolder = '20230906_145535_Helheim_AD_Schoof/'; % spc H
			case 1 % Weertman
				ADFolder = '20230713_202014_Helheim_AD_Weertman/';
			case 2 % Budd
				ADFolder = '20230703_145405_Helheim_AD_Budd/';
			otherwise
				error('The friction in the setting is not used.');
		end

		ADPath = [projPath, '/Models/', ADFolder, 'Model_Helheim_Big_Transient'];
		disp(['loading friction coefficients from AD results in ', ADPath])
		md_AD = loadmodel(ADPath);
		if isprop(md.friction, 'C')
			md.friction.C = md_AD.results.TransientSolution(1).FrictionC;
		else
			md.friction.coefficient = md_AD.results.TransientSolution(1).FrictionCoefficient;
		end

		% calving parameters
		md.calving = calvingvonmises();
		md.calving.stress_threshold_floatingice = 200*10^3;
		md.calving.stress_threshold_groundedice = sigma*10^6; %1: a little much retreat, 0.9,0.95: too much retreat 1.2, 1.1, 1.05,1.02:less retreat:

		% meltingrate
		timestamps = [md.timestepping.start_time, md.timestepping.final_time];
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		% load RACMO
		md.smb.mass_balance = [];
		for i= 2007:2019
			filename = ['/totten_1/ModelData/Greenland/RACMO23p2_2022_Greenland/smb/smb_rec.', num2str(i), '.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.nc'];
			disp(['Loading ' filename]);

			%Get time
			T  = double(ncread(filename,'time')); % in days
			if leapyear(i)
				time = i + (T)/366;
			else
				time = i + (T)/365;
			end

			%Coordinates
			x=double(ncread(filename,'x'));
			y=double(ncread(filename,'y'));

			%Load SMB for this year
			SMB=double(ncread(filename,'smb_rec'))/1000.0*12.0*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/month to mIE/yr

			%Interpolate
			for m=1:numel(time)
				smb_mesh=InterpFromGridToMesh(x, y, SMB(:,:,m)', md.mesh.x, md.mesh.y, 0);
				md.smb.mass_balance = [md.smb.mass_balance, [smb_mesh;time(m)]];
			end
		end

		savemodel(org,md);
	end%}}}
	if perform(org, ['Transient_RACMO_FromAD_', calving_suffix, suffix])% {{{

		md=loadmodel(org, ['Transient_Prep_RACMO_FromAD_', calving_suffix, suffix]);

		md.cluster = cluster;
		md.verbose.solution = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
	if perform(org, ['Transient_prep_extend_', smb_model, calving_suffix, suffix])% {{{

		md=loadmodel(org,['Transient_Prep_RACMO_FromAD_', calving_suffix, suffix]);
		% Step 1: set final_time to 2022
		md.timestepping.final_time = 2022;

		% Step 2: extend the experiment to 2022, since Calving Front from Greene goes to 2022
		% get observed calving front
		mask = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [md.timestepping.start_time, md.timestepping.final_time]);
		% convert icemask to levelset distance
		distance = zeros(size(mask));
		for i = 1:size(mask,2)
			distance(1:end-1,i) = reinitializelevelset(md, mask(1:end-1,i));
		end
		distance(end,:) = mask(end,:);
		% transient spc
		md.levelset.spclevelset = distance;

		% Step 3: RACMO also goes to 2022, but not MAR(only to 2020), for MAR we will repeat the average of the last three years
		md.smb.mass_balance = [];
		% meltingrate
		timestamps = [md.timestepping.start_time, md.timestepping.final_time];
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		if strcmp(smb_model, 'RACMO')
			% load RACMO
			for i= 2007:2021
				filename = ['/totten_1/ModelData/Greenland/RACMO23p2_2022_Greenland/smb/smb_rec.', num2str(i), '.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.nc'];
				disp(['Loading ' filename]);

				%Get time
				T  = double(ncread(filename,'time')); % in days
				if leapyear(i)
					time = i + (T)/366;
				else
					time = i + (T)/365;
				end

				%Coordinates
				x=double(ncread(filename,'x'));
				y=double(ncread(filename,'y'));

				%Load SMB for this year
				SMB=double(ncread(filename,'smb_rec'))/1000.0*12.0*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/month to mIE/yr

				%Interpolate
				for m=1:numel(time)
					smb_mesh=InterpFromGridToMesh(x, y, SMB(:,:,m)', md.mesh.x, md.mesh.y, 0);
					md.smb.mass_balance = [md.smb.mass_balance, [smb_mesh;time(m)]];
				end
			end
		elseif strcmp(smb_model, 'MAR')
			md.smb.mass_balance = [];
			for i=2007:2019
				filename = ['/totten_1/ModelData/Greenland/MARv3.11-ERA5/MARv3.11-monthly-ERA5-', num2str(i), '.nc'];
				disp(['Loading ' filename]);

				%Get time
				T  = double(ncread(filename,'time'));
				time = i + (T+0.5)/12;

				%Coordinates
				x=double(ncread(filename,'x'));
				y=double(ncread(filename,'y'));

				%Load SMB for this year
				SMB=double(ncread(filename,'SMB'))/1000.0*12.0*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/month to mIE/yr

				%Interpolate
				for m=1:numel(time)
					smb_mesh=InterpFromGridToMesh(x, y, SMB(:,:,m)', md.mesh.x, md.mesh.y, 0);
					md.smb.mass_balance = [md.smb.mass_balance, [smb_mesh;time(m)]];
				end
			end

			% Take average of the last three years from the MAR, repeat to the finalTime
			last3pos = find((md.smb.mass_balance(end,:)>=2017) & (md.smb.mass_balance(end,:)<2020));
			monthlypos = reshape(last3pos, 3, 12);

			averageSMB = zeros(md.mesh.numberofvertices, 12);
			for i = 1:12
				averageSMB(:,i) = mean(md.smb.mass_balance(1:end-1, monthlypos(:,i)),2);
			end
			repeatAverageSMB = repmat(averageSMB,1,2);

			% put this to md.smb.mass_balance for 2020-2022
			md.smb.mass_balance = [md.smb.mass_balance, [repeatAverageSMB;md.smb.mass_balance(end,1:24)+(2020-2007)]];
		elseif strcmp(smb_model, 'MAR_6k')
			md.smb.mass_balance = [];
			for i= md.timestepping.start_time: md.timestepping.final_time
				if (i>=1950) & (i<=2021)
					filename = ['/totten_1/ModelData/Greenland/MAR_6km/MARv3.11.5-Greenland-6km-daily-ERA-', num2str(i), '.nc'];
					disp(['Loading ' filename]);
				else
					disp(['Year ', num2str(i), ' is not coverd in the MAR data, skip for now!'])
					continue;
				end

				%Get time
				T0 = ncreadatt(filename, 'TIME', 'time_origin');
				T  = ncread(filename,'TIME');
				time = date2decyear(double(datenum(T0)+T));
				if ( (min(time) <i ) || (max(time)>i+1) )
					error(['TIME in ', filename, 'exceeds the year ', num2str(i)])
				end

				%Coordinates
				LAT=ncread(filename,'LAT');
				LON=ncread(filename,'LON');

				%Load SMB for this year
				SMB=double(ncread(filename,'SMB'))/1000*365.25*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/day to mIE/yr

				%Interpolate
				INDEX=BamgTriangulate(LON(:),LAT(:));
				SMB_reshaped = reshape(squeeze(SMB(:,:,1,:)),[numel(LAT) numel(time)]);
				smb_mar=InterpFromMeshToMesh2d(INDEX,LON(:),LAT(:),SMB_reshaped,md.mesh.long,md.mesh.lat);

				%Now use a monthly average
				smb_monthly = zeros(md.mesh.numberofvertices+1,1);
				for month=1:12
					pos = find(time>=time(1)+(month-1)/12  & time<time(1)+month/12);
					if ~isnan(mean(time(pos)))
						smb_monthly(1:end-1,month) = mean(smb_mar(:,pos),2);
						smb_monthly(end    ,month) = mean(time(pos));
					else
						disp(['Skipping year ' num2str(time(1)+month/12)]);
					end
				end

				md.smb.mass_balance = [md.smb.mass_balance smb_monthly];
			end
		else
			error(['Unknown SMB model: ', smb_model])
		end

		%Clean up
		savemodel(org,md);
	end%}}}
	if perform(org, ['Transient_extend_', smb_model, calving_suffix, suffix])% {{{

		md=loadmodel(org, ['Transient_prep_extend_', smb_model, calving_suffix, suffix]);

		md.cluster = cluster;
		md.verbose.solution = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end

		%Clean up
		savemodel(org,md);
	end%}}}

	% step 31--35
	if perform(org, ['Inversion_MOLHO_Budd'])% {{{
		md=loadmodel(org,['Transient_Prep_RACMO_FromAD_', calving_suffix, suffix]);
		md=setflowequation(md,'MOLHO','all');
      md.stressbalance.spcvx_base = md.stressbalance.spcvx;
      md.stressbalance.spcvy_base = md.stressbalance.spcvy;
      md.stressbalance.spcvx_shear = nan * md.stressbalance.spcvx;
      md.stressbalance.spcvy_shear = nan * md.stressbalance.spcvy;
		% remove unused fields
		md.smb.mass_balance = NaN;
		md.levelset.spclevelset = NaN;

		% fixed after add this option to effective pressure coupling
		md.friction.coupling = 2;

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
		md.inversion.control_parameters={'FrictionCoefficient'};
		md.inversion.maxsteps=400;
		md.inversion.maxiter =400;
		md.inversion.min_parameters=1e-5*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=1e4*ones(md.mesh.numberofvertices,1);
		md.inversion.control_scaling_factors=1;
		md.inversion.dxmin = 1e-6;

		%Additional parameters
		md.stressbalance.restol=1e-6;
		md.stressbalance.reltol=1e-4;
		md.stressbalance.abstol=NaN;

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
      md.settings.solver_residue_threshold = 1e-5;

		%Go solve
		md.cluster=cluster;
		md=solve(md,'sb');

		%Put results back into the model
		md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
		md.initialization.vx=md.results.StressbalanceSolution.Vx;
		md.initialization.vy=md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end%}}}
	if perform(org, ['Inversion_MOLHO_Schoof'])% {{{
      md=loadmodel(org,['Transient_Prep_RACMO_FromAD_', calving_suffix, suffix]);
      md=setflowequation(md,'MOLHO','all');
      md.stressbalance.spcvx_base = md.stressbalance.spcvx;
      md.stressbalance.spcvy_base = md.stressbalance.spcvy;
      md.stressbalance.spcvx_shear = nan * md.stressbalance.spcvx;
      md.stressbalance.spcvy_shear = nan * md.stressbalance.spcvy;
      % remove unused fields
      md.smb.mass_balance = NaN;
      md.levelset.spclevelset = NaN;

		%No friction on PURELY ocean element
		pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
		flags=ones(md.mesh.numberofvertices,1);
		flags(md.mesh.elements(pos_e,:))=0;
		md.friction.C(find(flags))=1e-2;

		fill_pos = find(min(md.initialization.vel(md.mesh.elements),[],2)==0); % fill-in postions
		ids = zeros(md.mesh.numberofvertices, 1);
		ids(md.mesh.elements(fill_pos,:))=1;

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
		pos=find((md.mask.ice_levelset>0) );
		md.inversion.cost_functions_coefficients(pos,1:2)=0;

		% set C=1e-2 for no obs region
		md.friction.C(ids>0) = 1e-2;
		md.inversion.cost_functions_coefficients(ids>0,1:2)=0;

		%Controls
		md.inversion.control_parameters={'FrictionC'};
		md.inversion.maxsteps=500;
		md.inversion.maxiter =500;
		md.inversion.min_parameters=1e-9*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=1e5*ones(md.mesh.numberofvertices,1);
		md.inversion.control_scaling_factors=3000;
		md.inversion.dxmin = 1e-6;
		%Additional parameters
		md.stressbalance.restol=1e-5;
		md.stressbalance.reltol=1e-3;
		md.stressbalance.abstol=NaN;
		md.stressbalance.maxiter = 100;

		%md.toolkits.DefaultAnalysis=mumpsoptions();
		%md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
		md.settings.solver_residue_threshold = 1e-5;
		%Go solve
		md.cluster=cluster;
		md=solve(md,'sb');

		%Put results back into the model
		md.friction.C=md.results.StressbalanceSolution.FrictionC;
		md.initialization.vx=md.results.StressbalanceSolution.Vx;
		md.initialization.vy=md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end%}}}
	if perform(org, ['Inversion_MOLHO_Weertman'])% {{{

      md=loadmodel(org,['Transient_Prep_RACMO_FromAD_', calving_suffix, suffix]);
      md=setflowequation(md,'MOLHO','all');
      md.stressbalance.spcvx_base = md.stressbalance.spcvx;
      md.stressbalance.spcvy_base = md.stressbalance.spcvy;
      md.stressbalance.spcvx_shear = nan * md.stressbalance.spcvx;
      md.stressbalance.spcvy_shear = nan * md.stressbalance.spcvy;
      % remove unused fields
      md.smb.mass_balance = NaN;
      md.levelset.spclevelset = NaN;

		%No friction on PURELY ocean element
		pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
		flags=ones(md.mesh.numberofvertices,1);
		flags(md.mesh.elements(pos_e,:))=0;
		md.friction.C(find(flags))=0.0;

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
		md.inversion.min_parameters=1e-6*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=5e4*ones(md.mesh.numberofvertices,1);
		md.inversion.control_scaling_factors=1;
		md.inversion.dxmin = 1e-6;
		%Additional parameters
		md.stressbalance.restol=1e-6;
		md.stressbalance.reltol=1e-5;
		md.stressbalance.abstol=NaN;

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
		md.settings.solver_residue_threshold = 1e-5;
		%Go solve
		md.cluster=cluster;
		md=solve(md,'sb');

		%Put results back into the model
		md.friction.C=md.results.StressbalanceSolution.FrictionC;
		md.initialization.vx=md.results.StressbalanceSolution.Vx;
		md.initialization.vy=md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end%}}}
	if perform(org, ['Transient_prep_MOLHO_', smb_model, suffix])% {{{

		md=loadmodel(org,['Inversion_MOLHO', friction_suffix]);
		% Step 1: set final_time to 2022
		md.timestepping.final_time = finalTime;

		% Step 2: extend the experiment to 2022, since Calving Front from Greene goes to 2022
		% get observed calving front
		mask = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [md.timestepping.start_time, md.timestepping.final_time]);
		% convert icemask to levelset distance
		distance = zeros(size(mask));
		for i = 1:size(mask,2)
			distance(1:end-1,i) = reinitializelevelset(md, mask(1:end-1,i));
		end
		distance(end,:) = mask(end,:);
		% transient spc
		md.levelset.spclevelset = distance;

		% Step 3: RACMO also goes to 2022, but not MAR(only to 2020), for MAR we will repeat the average of the last three years
		md.smb.mass_balance = [];
		% meltingrate
		timestamps = [md.timestepping.start_time, md.timestepping.final_time];
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
		md.frontalforcings.meltingrate(end,:) = timestamps;

		if strcmp(smb_model, 'RACMO')
         % load RACMO
         for i= md.timestepping.start_time: md.timestepping.final_time
            if (i>=1990) & (i<=2021)
               filename = ['/totten_1/ModelData/Greenland/RACMO23p2_2022_Greenland/smb/smb_rec.', num2str(i), '.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.nc'];

               %Get time
               T  = double(ncread(filename,'time')); % in days
               if leapyear(i)
                  time = i + (T)/366;
               else
                  time = i + (T)/365;
               end

               %Coordinates
               x=double(ncread(filename,'x'));
               y=double(ncread(filename,'y'));

               %Load SMB for this year
               disp(['Loading ' filename]);
               SMB=double(ncread(filename,'smb_rec'))/1000.0*12.0*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/month to mIE/yr
            elseif (i>= 1958) & (i<=1989)
               filename = ['/totten_1/ModelData/Greenland/RACMO23p2_2022_Greenland/smb/smb_rec.', num2str(i), '.BN_RACMO2.3p2_ERA40_ERAIn_FGRN055.1km.MM.nc'];
               %Get time
               T  = double(ncread(filename,'time')); % in month
               time = i + (T+0.5)/12;

               %Coordinates
               x=double(ncread(filename,'lon'));
               y=double(ncread(filename,'lat'));

               %Load SMB for this year
               disp(['Loading ' filename]);
               SMB=double(ncread(filename,'SMB_rec'))/1000.0*12.0*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/month to mIE/yr
            else
               disp(['Year ', num2str(i), ' is not coverd in RACMO data, skip for now!'])
               continue;
            end

            %Interpolate
            for m=1:numel(time)
               smb_mesh=InterpFromGridToMesh(x, y, SMB(:,:,m)', md.mesh.x, md.mesh.y, 0);
               md.smb.mass_balance = [md.smb.mass_balance, [smb_mesh;time(m)]];
            end
         end
		elseif strcmp(smb_model, 'MAR')
         for i= md.timestepping.start_time: md.timestepping.final_time
            if (i>=1950) & (i<=2019)
               filename = ['/totten_1/ModelData/Greenland/MARv3.11-ERA5/MARv3.11-monthly-ERA5-', num2str(i), '.nc'];
               disp(['Loading ' filename]);
            else
               disp(['Year ', num2str(i), ' is not coverd in the MAR data, skip for now!'])
               continue;
            end

            %Get time
            T  = double(ncread(filename,'time'));
            time = i + (T+0.5)/12;

            %Coordinates
            x=double(ncread(filename,'x'));
            y=double(ncread(filename,'y'));

            %Load SMB for this year
            SMB=double(ncread(filename,'SMB'))/1000.0*12.0*md.materials.rho_freshwater/md.materials.rho_ice; %from mmWE/month to mIE/yr

            %Interpolate
            for m=1:numel(time)
               smb_mesh=InterpFromGridToMesh(x, y, SMB(:,:,m)', md.mesh.x, md.mesh.y, 0);
               md.smb.mass_balance = [md.smb.mass_balance, [smb_mesh;time(m)]];
            end
         end

         if (md.timestepping.final_time>2020)
            % Take average of the last three years from the MAR, repeat to the finalTime
            last3pos = find((md.smb.mass_balance(end,:)>=2017) & (md.smb.mass_balance(end,:)<2020));
            monthlypos = reshape(last3pos, 12, 3);

            averageSMB = zeros(md.mesh.numberofvertices, 12);
            for i = 1:12
               averageSMB(:,i) = mean(md.smb.mass_balance(1:end-1, monthlypos(i,:)),2);
            end
            Next = ceil(md.timestepping.final_time - 2020);
            repeatAverageSMB = repmat(averageSMB,1,Next);
            repeatTime = 2020+linspace(0.5, Next*12-0.5,Next*12)/12;

            % put this to md.smb.mass_balance for 2020-2022
            md.smb.mass_balance = [md.smb.mass_balance, [repeatAverageSMB;repeatTime]];
         end
		else
			error(['Unknown SMB model: ', smb_model])
		end

		%Clean up
		savemodel(org,md);
	end%}}}
	if perform(org, ['Transient_MOLHO_', smb_model, suffix])% {{{

		md=loadmodel(org, ['Transient_prep_MOLHO_', smb_model, suffix]);

		% Set parameters
      md.inversion.iscontrol=0;
      md.settings.output_frequency = 1;
		md.timestepping.final_time = finalTime;

		md.cluster = cluster;
		md.stressbalance.restol = 1e-4;

		md.verbose.solution = 1;
		md.verbose.convergence = 1;
		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.settings.solver_residue_threshold = 1e-5;

		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM', 'CalvingMeltingrate'};

		md=solve(md,'Transient','runtimename',false);
		savemodel(org,md);

		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['cp ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end%}}}
