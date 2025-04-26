function fric=EstimateFric_Weertman(md) % {{{

   disp('Estimating friction coefficient - Weertman')
   % initial guess from driving stress
   min_velocity=0.1; % m/yr
   min_c=1;
   m=3;

   asurf=averaging(md,md.geometry.surface,20); % maybe executing 20 L2 projection is ok
   [sx,sy,s]=slope(md,asurf); % slope 's' comes on elements
   sslope=averaging(md,s,1); % average the slope once on the vertices, because 's' comes on elements, we need this data on vertices

   vel=md.inversion.vel_obs;
   vel=max(vel,min_velocity); % setting minimum velocity value
   vel=vel/md.constants.yts; %S.I.;

   driving_stress=md.materials.rho_ice*md.constants.g*md.geometry.thickness.*(sslope);
   c=sqrt(driving_stress./(vel.^(1/m)));
   c=max(c,min_c);

   fric=c; % initial guess

end%}}}
