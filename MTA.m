function Sol=MTA(MTAstruct)

%**************************************************************************
%                      METABOLIC TRANSITION ANALYSIS
%**************************************************************************
%
%This function implements the MTA simulating flux and concetration dynamics
%of a genome-scale model (GEM). The function takes a structure as input
%(MTAstruct) which contains the following required and optional fields:
%
%Required
%
%  GEM: genome-scale model
%
%  Vobj, Xobj,DataMinObj: vectors with the indexes of the reactions (Vobj)
%  or metabolites (Xobj) whose value should be maximized (minimized if
%  optional field ObjSense='min'). 
%  DataMinObj is logical, either 'T' or 'F', and indicates if the objective 
%  should be minimizing the bounds to the supplied experimental data 
%  (see optional fields Vt, Xt). At least one of the three fields must be 
%  non-empty. If DataMinObj = 'T', then at least one of the optional fields 
%  Vt or Xt must be supplied.
%
%Optional
%
%  Vt: a three column matrix, the first column displays the flux values of
%  a set of reactions, whose indeces in the GEM are specified in the second
%  column. The third column specifies the time step in which these flux
%  values should be fixed.For example, a row Vti in Vt may be: Vti=[1,4,1],
%  meaning that reaction 4 in the GEM should have a fixed flux value of 1
%  in time step 1 (which is the initial time step)
%
%  Xt: a three column matrix similar to Vt, in this case, the first column
%  displays the concentration values of a set of metabolites whose indices
%  in the GEM are specified in the second column. The third column
%  specifies the time step where the constraint should be applied. If Xt is
%  present, then the GEM must be pre-modified: a new metabolite (row in
%  GEM.S) must be included per each metabolite with data, this corresponds
%  to the total pool of the metabolite. This modification is added since
%  the data normally corresponds to the total pool measurement of a
%  metabolite, that is, regardless of its location in a subcellular
%  compartment. The new row, is created by summing all rows containing the
%  metabolite in different compartmens. Metabolite names, GEM.metNames must
%  be modified accordingly, as well as the indexes in Xt to much the newly
%  created pool metabolites
%
%  nsteps: a scalar indicating the number of time steps in the simulation
%
%  MaxX, MaxFlux, MaxRate,MaxDelta: scalars, maximum allowed values for
%  concentration, flux, rate and rate (or flux) increment per time step,
%  respectively. 
%
%  lbVt, ubVt, lbXt, ubXt, lbRt, ubRt: these are three column matrices
%  similar to Vt and Xt but specifying flux, concentration or concentration
%  rate bounds,respectively, for particular reactions/metabolites in
%  particular time points. First column indicates the bound value, second
%  column the index of the reaction or the metabolite in the GEM and third
%  column the time step where the bound should be applied.
%
%  TotalObj: logical, either 'T' of 'F', indicates if the objective
%  function should be optimized accross the time steps, 'T', (optimize the
%  sum) or only at the last time step, 'F'
%
%  SteadyState: logical, either 'T' or 'F', indicates if a steady state
%  (Sv=0) is impossed at the last time step or not, respectively
%
%  Coupling: character string, can be 'Fluxes' or 'Rates', indicating if
%  flux values are coupled (i.e. vi+1=vi+deltai) or rate values are coupled
%  (i.e. bi+1=bi+deltai, where Sv=b) between time steps. It can also be set
%  to 'F', in which case flux and rate values are uncoupled (i.e. each
%  value at time step i is independent from the preceding time steps). In
%  the latter case, the field "MaxRate" can be used to control the maximum
%  metabolite concentration change between consecutive time steps, and the
%  field MaxDelta is not used.
%
%  Mets2Plot, Rxns2Plot, Rates2Plot: numeric vectors containing the
%  indexes in the GEM of the metabolites, reactions, and rates (i.e.
%  metabolites in the GEM) that should be plotted. Default is empty unless
%  the alternative optima sampling method, AOSanalysis=T (see below), in
%  which case metabolites with a variablity value less than the 5%
%  percentile are plotted.
%
%  ObjSense: character string, either 'max' or 'min', to maximize or
%  minimize the objective function (Default is 'max' unless DataMinObj is
%  selected)
%
%  MaxEpsilon: a vector displaying the maximum (upper and lower) deviations
%  from the experimental data. If the experimental data are concentrations,
%  the length should be (nsteps+1)*NumberDataPoints. If fluxes then
%  nsteps*NumberDataPoints. If both then first concentrations then fluxes.
%
%  Tasks: an array of reaction indexes in the GEM. These reactions are
%  constrained to carry a non-zero flux (currently only works with
%  irreversible reactions)
%
%  Blocked: an array of reaction indexes in the GEM. These reactions are
%  constrained to carry 0 flux
%
%  Constant: an array of metabolite indexes in the GEM. These metabolites
%  are constrained to be at steady state throughout the simulation (i.e,
%  dx/dt = 0)
%
%  Emax: The maximum value of the variable epsilon in the AO sampling.
%  Default value is set to 1e6
%
%  Zoptimum: The value of the objective function corresponding to a
%  previously found MTA optimum. It is required to conduct the AO analysis
%
%  FVAanalysis: The type of AO analysis to be conducted. If 'T', then a
%  flux variability analysis (min and max values) is conducted.
%
%  AOSanalysis: The type of AO analysis to be conducted. If 'T', then an an
%  alternative optima sampling procedure is conducted.The AOS only searches
%  for closes solutions to a random concentration vector Xrand (no random
%  flux values are generated)
%
%  Mets2EvalFVA: a numerical array containing the indeces in the GEM of the
%  metabolites for which the FVA analysis is applied to obtain their AO min
%  and max values. Only in case FVAanalysis='T'
%
%  Rxns2EvalFVA: same as the previous field but containing the indeces in
%  the GEM of the reactions for which the FVA analysis is applied.
%
%  AOSnsamples: the number of sampled AO solutions, default is set to 1e2.
%  This field is only required if AOanalysis is set to 'AOS'
%
%  saveAOsamples, saveFVAranges: either 'T' or 'F' indicating if the
%  generated AO sampled trajectories or the concentration/flux ranges
%  should be saved to a mat file
%
%  Crossover: either 'T' or 'F', indicating if the crossover method should
%  be applied to the barrier algorithm in Gurobi. Including crossover may
%  increase accuracy at the cost of computational time
%
%**************************************************************************
%           Semidan, July 2016 (robaina@mpimp-golm.mpg.de)
%**************************************************************************

    %Argument evaluation
    RevRxns=find(MTAstruct.GEM.rev==1);
    if ~isfield(MTAstruct,'ConcFluxConstraint'),
        MTAstruct.ConcFluxConstraint='F';
    end   
    if MTAstruct.ConcFluxConstraint=='T',
        OrRxns=size(MTAstruct.GEM.S,2);
        MTAstruct.GEM.S=[MTAstruct.GEM.S,-MTAstruct.GEM.S(:,RevRxns)];
        MTAstruct.GEM.rev=zeros(size(MTAstruct.GEM.S,2),1);
    end
    Rxns=size(MTAstruct.GEM.S,2);
    Mets=size(MTAstruct.GEM.S,1);
    if ~isfield(MTAstruct,'DataMinObj'),
        MTAstruct.DataMinObj='F';
    end
    if ~isfield(MTAstruct,'GEM'),
        error('GEM is missing!')
    end
    if ~isfield(MTAstruct,'Vobj') && ~isfield(MTAstruct,'Xobj') && MTAstruct.DataMinObj=='F',
        error('Objective is missing!')
    end
    if MTAstruct.DataMinObj=='T' && ~isfield(MTAstruct,'Xt') && ~isfield(MTAstruct,'Vt'),
        error('Experimental data are missing!')
    end
    if ~isfield(MTAstruct,'nsteps'),
        MTAstruct.nsteps=10;
    end
    if ~isfield(MTAstruct,'MaxX'),
        MTAstruct.MaxX=10;
    end
    if ~isfield(MTAstruct,'MaxFlux'),
        MTAstruct.MaxFlux=1; 
    end
    if ~isfield(MTAstruct,'MaxRate'),
        MTAstruct.MaxRate=1; 
    end
    if ~isfield(MTAstruct,'Coupling'),
        MTAstruct.Coupling='Fluxes';
    end
    if ~isfield(MTAstruct,'MaxDelta') && MTAstruct.Coupling~='F',
        MTAstruct.MaxDelta=1e-3; 
    end
    if ~isfield(MTAstruct,'TotalObj'),
        MTAstruct.TotalObj='T';
    end
    if ~isfield(MTAstruct,'SteadyState'),
        MTAstruct.SteadyState='T';
    end
    if ~isfield(MTAstruct,'ObjSense'), 
        MTAstruct.ObjSense='max';
    end
    if MTAstruct.DataMinObj=='F' && ~isfield(MTAstruct,'MaxEpsilon'),
        D=.1;
    elseif MTAstruct.DataMinObj=='T' && ~isfield(MTAstruct,'MaxEpsilon'),
        D=1;
    end
    if isfield(MTAstruct,'Xt') && ~isfield(MTAstruct,'Vt'),
        if ~isfield(MTAstruct,'MaxEpsilon'),
            MTAstruct.MaxEpsilon=D*MTAstruct.Xt(:,1);
        end
    elseif ~isfield(MTAstruct,'Xt') && isfield(MTAstruct,'Vt'),
        if ~isfield(MTAstruct,'MaxEpsilon'),
            MTAstruct.MaxEpsilon=D*MTAstruct.Vt(:,1);
        end
    elseif isfield(MTAstruct,'Xt') && isfield(MTAstruct,'Vt'),
        if ~isfield(MTAstruct,'MaxEpsilon'),
            MTAstruct.MaxEpsilon=[D*MTAstruct.Xt(:,1);D*MTAstruct.Vt(:,1)];
        end
    end
    if strcmpi(MTAstruct.Coupling,'Fluxes'),
        DimDelta=Rxns;
    elseif strcmpi(MTAstruct.Coupling,'Rates'),
        DimDelta=Mets;
    elseif strcmpi(MTAstruct.Coupling,'F'),
        DimDelta=0;
    end
     if ~isfield(MTAstruct,'Emax'), 
        MTAstruct.Emax=1e6;
     end
     if ~isfield(MTAstruct,'FVAanalysis'),
         MTAstruct.FVAanalysis='F';
     end
     if ~isfield(MTAstruct,'AOSanalysis'),
         MTAstruct.AOSanalysis='F';
     end 
     if isfield(MTAstruct,'Zoptimum') && ~MTAstruct.AOSanalysis && ~MTAstruct.FVAanalysis,
         error('Please provide type of alternative optima analysis');
     end
     if ~(MTAstruct.AOSanalysis=='F' && MTAstruct.FVAanalysis=='F'),
         if ~isfield(MTAstruct,'Zoptimum'),
              error('Optimal objective value is missing!');
         end
         if MTAstruct.FVAanalysis=='T' && ~isfield(MTAstruct,'Mets2EvalFVA') && ~isfield(MTAstruct,'Rxns2EvalFVA'),
            error('Metabolites or Reactions to evaluate per FVA are missing!')
         end
         if MTAstruct.AOSanalysis && ~isfield(MTAstruct,'AOSnsamples'),
             MTAstruct.AOSnsamples=1e2;
         end
     end
     if ~isfield(MTAstruct,'Mets2Plot'),
        if MTAstruct.AOSanalysis=='F',
            MTAstruct.Mets2Plot=[];
        elseif MTAstruct.AOSanalysis=='T',
            MetPercentile=0.05;
        end
    end
    if ~isfield(MTAstruct,'Rxns2Plot'),
        MTAstruct.Rxns2Plot=[];
    end
    if ~isfield(MTAstruct,'Rates2Plot'),
        MTAstruct.Rates2Plot=[];
    end
    if ~isfield(MTAstruct,'saveAOsamples'),
        MTAstruct.saveAOsamples='F';
    end
    if ~isfield(MTAstruct,'saveDaos'),
        MTAstruct.saveDaos='F';
    end
    if ~isfield(MTAstruct,'saveFVAranges'),
        MTAstruct.saveFVAranges='F';
    end
    if ~isfield(MTAstruct,'Crossover'),
        MTAstruct.Crossover='T';
    end
    if ~isfield(MTAstruct,'FigDir'),
        MTAstruct.FigDir=[];
    end
    if ~isfield(MTAstruct,'AOSvarDist'),
        MTAstruct.AOSvarDist='F';
    end
    if ~isfield(MTAstruct,'ConcFluxConstraint'),
        MTAstruct.ConcFluxConstraint='F';
    end
    if ~isfield(MTAstruct,'ConcThreshold'),
        MTAstruct.ConcThreshold=1e3;
    end
         
    %construct Amat matrix 
    for k=1:MTAstruct.nsteps
      C1{k}=MTAstruct.GEM.S;
      C2{k}=-speye(Mets);
    end
    A1=sparse([blkdiag(C1{:}),blkdiag(C2{:})]);A1=[A1,sparse(size(A1,1),(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta)]; %Svi-bi=0    
    A31=sparse((MTAstruct.nsteps-1)*Mets,MTAstruct.nsteps*Mets);
    w=1;
    for k=1:MTAstruct.nsteps,
       A31(w:k*Mets,w:w+2*Mets-1)=[speye(Mets),-speye(Mets)];
       w=k*Mets+1;
    end
    A32=sparse([-blkdiag(C2{:}),A31]);
    A3=[sparse(MTAstruct.nsteps*Mets,MTAstruct.nsteps*Rxns),A32,sparse(MTAstruct.nsteps*Mets,(MTAstruct.nsteps-1)*DimDelta)]; %x1+b1-x2=0...
    if strcmpi(MTAstruct.Coupling,'F')==0,
        A21=sparse((MTAstruct.nsteps-1)*DimDelta,MTAstruct.nsteps*DimDelta);
        w=1;
        for k=1:MTAstruct.nsteps-1,
           C3{k}=-speye(DimDelta); 
           A21(w:k*DimDelta,w:w+2*DimDelta-1)=[speye(DimDelta),-speye(DimDelta)]; 
           w=k*DimDelta+1;
        end
        if strcmpi(MTAstruct.Coupling,'Fluxes'), 
            A2=[A21,sparse(size(A21,1),MTAstruct.nsteps*(2*Mets)+Mets),sparse(blkdiag(C3{:}))]; %v(i+1)-v(i)-delta(i)=0
        elseif strcmpi(MTAstruct.Coupling,'Rates'),
            A2=[sparse(size(A21,1),MTAstruct.nsteps*Rxns),A21,sparse(size(A21,1),(MTAstruct.nsteps+1)*Mets),sparse(blkdiag(C3{:}))]; %b(i+1)-b(i)-delta(i)=0
        end
        Amat=[A1;A2;A3]; 
    elseif MTAstruct.Coupling=='F'
       Amat=[A1;A3]; 
       A2=[];
    end
    
    %set bounds on reaction fluxes, metabolite rates and concentrations
    lbrates=-MTAstruct.MaxRate*ones(Mets,1);ubrates=MTAstruct.MaxRate*ones(Mets,1);
    lbdeltas=-MTAstruct.MaxDelta*ones(DimDelta,1);ubdeltas=MTAstruct.MaxDelta*ones(DimDelta,1);
    lbrxns=zeros(Rxns,1);lbrxns(MTAstruct.GEM.rev==1)=-MTAstruct.MaxFlux;
    ubrxns=MTAstruct.MaxFlux*ones(Rxns,1);
    
    if isfield(MTAstruct,'Tasks'), 
        lbrxns(MTAstruct.Tasks)=MTAstruct.MaxFlux*1e-6; 
    end
    if isfield(MTAstruct,'Blocked'),
        lbrxns(MTAstruct.Blocked)=0;
        ubrxns(MTAstruct.Blocked)=0;
    end
    if isfield(MTAstruct,'Constant'),
        lbrates(MTAstruct.Constant)=-1e-6;
        ubrates(MTAstruct.Constant)=1e-6;
    end
    
    %impose steady state at last time step?
    if MTAstruct.SteadyState=='T',
        FinalRatelb=zeros(Mets,1);
        FinalRateub=zeros(Mets,1);
    elseif MTAstruct.SteadyState=='F',
        FinalRatelb=lbrates;
        FinalRateub=ubrates;
    end
    
    totlbrxns=repmat(lbrxns,MTAstruct.nsteps,1);
    totubrxns=repmat(ubrxns,MTAstruct.nsteps,1);
    totlbmets=zeros(Mets*(MTAstruct.nsteps+1),1);
    totubmets=MTAstruct.MaxX*ones(Mets*(MTAstruct.nsteps+1),1);
    totlbrates=[repmat(lbrates,MTAstruct.nsteps-1,1);FinalRatelb];
    totubrates=[repmat(ubrates,MTAstruct.nsteps-1,1);FinalRateub];
   
    %include time-resolved bounds
    if isfield(MTAstruct,'lbVt'),
       for i=1:size(MTAstruct.lbVt,1),
           totlbrxns(MTAstruct.lbVt(i,2)+Rxns*(MTAstruct.lbVt(i,3)-1))=MTAstruct.lbVt(i,1);
       end
    end
    if isfield(MTAstruct,'ubVt'),
       for i=1:size(MTAstruct.ubVt,1),
           totubrxns(MTAstruct.ubVt(i,2)+Rxns*(MTAstruct.ubVt(i,3)-1))=MTAstruct.ubVt(i,1);
       end
    end
    if isfield(MTAstruct,'lbXt'),
       for i=1:size(MTAstruct.lbXt,1),
           totlbmets(MTAstruct.lbXt(i,2)+Mets*(MTAstruct.lbXt(i,3)-1))=MTAstruct.lbXt(i,1);
       end
    end
    if isfield(MTAstruct,'ubXt'),
       for i=1:size(MTAstruct.ubXt,1),
           totubmets(MTAstruct.ubXt(i,2)+Mets*(MTAstruct.ubXt(i,3)-1))=MTAstruct.ubXt(i,1);
       end
    end
    if isfield(MTAstruct,'lbRt'),
       for i=1:size(MTAstruct.lbRt,1),
           totlbrates(MTAstruct.lbRt(i,2)+Mets*(MTAstruct.lbRt(i,3)-1))=MTAstruct.lbRt(i,1);
       end
    end
    if isfield(MTAstruct,'ubRt'),
       for i=1:size(MTAstruct.ubRt,1),
           totubrates(MTAstruct.ubRt(i,2)+Mets*(MTAstruct.ubRt(i,3)-1))=MTAstruct.ubRt(i,1);
       end
    end
      
    lbvec=[totlbrxns;totlbrates;totlbmets;repmat(lbdeltas,MTAstruct.nsteps-1,1)];
    ubvec=[totubrxns;totubrates;totubmets;repmat(ubdeltas,MTAstruct.nsteps-1,1)];
    
    if ~isfield(MTAstruct,'Xt') && ~isfield(MTAstruct,'Vt'),
         rhsvec=zeros(size(Amat,1),1);
         sensevec=repmat('=',size(Amat,1),1);
    end
    
    if isfield(MTAstruct,'Vobj') || isfield(MTAstruct,'Xobj'),
         %construct objective vector
         cV=zeros(Rxns,1);
         cX=zeros(Mets,1);
         if isfield(MTAstruct,'Vobj'),
             cV(MTAstruct.Vobj)=1;
         end
         if isfield(MTAstruct,'Xobj'),
             cX(MTAstruct.Xobj)=1;
         end
         if MTAstruct.TotalObj=='T',
            ObjVector=[repmat(cV,MTAstruct.nsteps,1);zeros(MTAstruct.nsteps*Mets,1);repmat(cX,MTAstruct.nsteps+1,1);zeros((MTAstruct.nsteps-1)*DimDelta,1)];
         elseif MTAstruct.TotalObj=='F',
             ObjVector=[zeros((MTAstruct.nsteps-1)*Rxns,1);cV;zeros(MTAstruct.nsteps*Mets,1);zeros(MTAstruct.nsteps*Mets,1);cX;zeros((MTAstruct.nsteps-1)*DimDelta,1)];
         end
    end
    
    %Impose specific bounds throughout time steps, update solver structure
  
    %  metabolite concentrations
    if isfield(MTAstruct,'Xt') && ~isfield(MTAstruct,'Vt'),
       A41=sparse(size(MTAstruct.Xt,1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+size(MTAstruct.Xt,1));A42=A41;
       for i=1:size(MTAstruct.Xt,1),
           A41(i,[MTAstruct.nsteps*(Rxns+Mets)+MTAstruct.Xt(i,2)+Mets*(MTAstruct.Xt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+i])=[1,1];
           A42(i,[MTAstruct.nsteps*(Rxns+Mets)+MTAstruct.Xt(i,2)+Mets*(MTAstruct.Xt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+i])=[1,-1];
       end 
       Amat=[[Amat,sparse(size(Amat,1),size(MTAstruct.Xt,1))];A41;A42];
       
       if isfield(MTAstruct,'Vobj') || isfield(MTAstruct,'Xobj'),
           ObjVector=[ObjVector;zeros(size(MTAstruct.Xt,1),1)];
       elseif MTAstruct.DataMinObj=='T',
           ObjVector=[zeros(MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta,1);ones(size(MTAstruct.Xt,1),1)];   
       end
       lbvec=[lbvec;-MTAstruct.MaxEpsilon];ubvec=[ubvec;MTAstruct.MaxEpsilon]; 
       rhsvec=[zeros(size(A1,1)+size(A2,1)+size(A3,1),1);MTAstruct.Xt(:,1);MTAstruct.Xt(:,1)];
       sensevec=[repmat('=',size(A1,1)+size(A2,1)+size(A3,1),1);repmat('>',size(MTAstruct.Xt,1),1);repmat('<',size(MTAstruct.Xt,1),1)];          
    end
    
    %  reaction fluxes
    if isfield(MTAstruct,'Vt') && ~isfield(MTAstruct,'Xt'),
      A41=sparse(size(MTAstruct.Vt,1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+size(MTAstruct.Vt,1));A42=A41;
      for i=1:size(MTAstruct.Vt,1),
          A41(i,[MTAstruct.Vt(i,2)+Rxns*(MTAstruct.Vt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+i])=[1,1];
          A42(i,[MTAstruct.Vt(i,2)+Rxns*(MTAstruct.Vt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+i])=[1,-1];
      end
      Amat=[[Amat,sparse(size(Amat,1),size(MTAstruct.Vt,1))];A41;A42];
      
      if isfield(MTAstruct,'Vobj') || isfield(MTAstruct,'Xobj'),
           ObjVector=[ObjVector;zeros(size(MTAstruct.Vt,1),1)];
      elseif MTAstruct.DataMinObj=='T',
           ObjVector=[zeros(MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta,1);ones(size(MTAstruct.Vt,1),1)];   
      end
      lbvec=[lbvec;-MTAstruct.MaxEpsilon];ubvec=[ubvec;MTAstruct.MaxEpsilon]; 
      rhsvec=[zeros(size(A1,1)+size(A2,1)+size(A3,1),1);MTAstruct.Vt(:,1);MTAstruct.Vt(:,1)];
      sensevec=[repmat('=',size(A1,1)+size(A2,1)+size(A3,1),1);repmat('>',size(MTAstruct.Vt,1),1);repmat('<',size(MTAstruct.Vt,1),1)];    
    end
    
    %  metabolite concentrations and reaction fluxes requires creating
    %  deltas for rxns and mets
    if isfield(MTAstruct,'Xt') && isfield(MTAstruct,'Vt'),
       A41=sparse(size(MTAstruct.Xt,1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+size(MTAstruct.Xt,1)+size(MTAstruct.Vt,1));A42=A41;
       A43=sparse(size(MTAstruct.Vt,1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+size(MTAstruct.Xt,1)+size(MTAstruct.Vt,1));A44=A43;
       for i=1:size(MTAstruct.Xt,1),
           A41(i,[MTAstruct.nsteps*(Rxns+Mets)+MTAstruct.Xt(i,2)+Mets*(MTAstruct.Xt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+i])=[1,1];
           A42(i,[MTAstruct.nsteps*(Rxns+Mets)+MTAstruct.Xt(i,2)+Mets*(MTAstruct.Xt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+i])=[1,-1];
       end
       for i=1:size(MTAstruct.Vt,1),
           A43(i,[MTAstruct.Vt(i,2)+Rxns*(MTAstruct.Vt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+size(MTAstruct.Xt,1)+i])=[1,1];
           A44(i,[MTAstruct.Vt(i,2)+Rxns*(MTAstruct.Vt(i,3)-1),MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta+size(MTAstruct.Xt,1)+i])=[1,-1];
       end
       Amat=[[Amat,sparse(size(Amat,1),size(MTAstruct.Xt,1)+size(MTAstruct.Vt,1))];A41;A42;A43;A44];
       
       if isfield(MTAstruct,'Vobj') || isfield(MTAstruct,'Xobj'),
           ObjVector=[ObjVector;zeros(size(MTAstruct.Xt,1)+size(MTAstruct.Vt,1),1)];
       elseif MTAstruct.DataMinObj=='T',
           ObjVector=[zeros(MTAstruct.nsteps*(Rxns+Mets)+(MTAstruct.nsteps+1)*Mets+(MTAstruct.nsteps-1)*DimDelta,1);ones(size(MTAstruct.Xt,1)+size(MTAstruct.Vt,1),1)];   
       end
       lbvec=[lbvec;-MTAstruct.MaxEpsilon];ubvec=[ubvec;MTAstruct.MaxEpsilon]; 
       rhsvec=[zeros(size(A1,1)+size(A2,1)+size(A3,1),1);MTAstruct.Xt(:,1);MTAstruct.Xt(:,1);MTAstruct.Vt(:,1);MTAstruct.Vt(:,1)];
       sensevec=[repmat('=',size(A1,1)+size(A2,1)+size(A3,1),1);repmat('>',size(MTAstruct.Xt,1),1);repmat('<',size(MTAstruct.Xt,1),1);repmat('>',size(MTAstruct.Vt,1),1);repmat('<',size(MTAstruct.Vt,1),1)];     
    end
    
    if MTAstruct.ConcFluxConstraint=='T',
        %Impose null flux values when concentration of substrates is zero at
        %each time step
        A5=sparse(Mets,size(Amat,2));
        for j=0:(MTAstruct.nsteps-1),
           for i=1:Mets,
               FedRxns=find(MTAstruct.GEM.S(i,:)<0);
               A5(i+(j*Mets),[FedRxns+(j*Rxns),((Rxns+Mets)*MTAstruct.nsteps)+i+(j*Mets)])=[ones(1,length(FedRxns)),-MTAstruct.ConcThreshold];
           end
        end
          Amat=[Amat;A5];
          sensevec=[sensevec;repmat('<',size(A5,1),1)];
          rhsvec=[rhsvec;zeros(size(A5,1),1)];
    end
     
     if MTAstruct.FVAanalysis=='F' && MTAstruct.AOSanalysis=='F',
        %solve LP
        m.obj=ObjVector;
        m.A=Amat;
        m.modelsense=MTAstruct.ObjSense;
        m.lb=lbvec;
        m.ub=ubvec;
        m.rhs=rhsvec;
        m.sense=sensevec;
        params.outputflag=1;
        params.Threads=8;
        params.method=2; %0=primal simplex;1=dual simplex;2=barrier method
%         params.BarConvTol=0;
        if MTAstruct.Crossover=='F',
           params.crossover=0;
        end
        gur=gurobi(m,params);

        %get solution
        n=1;p=1;
        bmatrix=zeros(Mets,MTAstruct.nsteps);
        if MTAstruct.ConcFluxConstraint=='T',
            Vmatrix=zeros(OrRxns,MTAstruct.nsteps);
            for k=1:MTAstruct.nsteps,
                V=gur.x(n:k*Rxns);
                Vmatrix(:,k)=V(1:OrRxns);Vmatrix(RevRxns,k)=Vmatrix(RevRxns,k)-V((OrRxns+1):end);
                bmatrix(:,k)=gur.x(MTAstruct.nsteps*Rxns+p:MTAstruct.nsteps*Rxns+k*Mets);
                n=n+Rxns;
                p=p+Mets;
            end
        elseif MTAstruct.ConcFluxConstraint=='F',
            Vmatrix=zeros(Rxns,MTAstruct.nsteps);
            for k=1:MTAstruct.nsteps,
                Vmatrix(:,k)=gur.x(n:k*Rxns);
                bmatrix(:,k)=gur.x(MTAstruct.nsteps*Rxns+p:MTAstruct.nsteps*Rxns+k*Mets);
                n=n+Rxns;
                p=p+Mets;
            end
        end
        Xmatrix=zeros(Mets,MTAstruct.nsteps+1);
        p=1;
        for k=1:MTAstruct.nsteps+1,
            Xmatrix(:,k)=gur.x(MTAstruct.nsteps*(Rxns+Mets)+p:MTAstruct.nsteps*(Rxns+Mets)+k*Mets);
            p=p+Mets;
        end
        Xmatrix(abs(Xmatrix)<1e-8)=0;
        Vmatrix(abs(Vmatrix)<1e-8)=0;
        Sol.Fluxes=Vmatrix;
        Sol.Rates=bmatrix;
        Sol.Concentrations=Xmatrix;
        Sol.Zoptimum=gur.objval;

        %plot figures
        u=1;
        if ~isempty(MTAstruct.Mets2Plot) || ~isempty(MTAstruct.Rates2Plot),
            if isfield(MTAstruct.GEM,'metNames'),
                MetNames=MTAstruct.GEM.metNames;
            elseif ~isfield(MTAstruct.GEM,'metNames'),
                for i=1:size(MTAstruct.GEM.S,1),
                    MetNames{i}=['x',num2str(i)];
                end
            end
        end
        if ~isempty(MTAstruct.Mets2Plot),
            figure(u)
            plot(1:(MTAstruct.nsteps+1),Xmatrix(MTAstruct.Mets2Plot,:))
            if isfield(MTAstruct,'Xt') && ~isempty(intersect(MTAstruct.Mets2Plot,MTAstruct.Xt(:,2))),
                DataPoints=find(ismember(MTAstruct.Xt(:,2),MTAstruct.Mets2Plot));
                Xaxis=zeros(length(DataPoints),1);
                for k=1:length(DataPoints),
                    Xaxis(k)=MTAstruct.Xt(DataPoints(k),3);
                end
                hold on %Hold the Door! Hold Door! HolDoor! Hodoor! Hodor!
                plot(Xaxis,MTAstruct.Xt(DataPoints,1),'*k','MarkerSize',7)
            end
            ylabel('Concentration');
            xlabel('time steps')
            legend(MetNames(MTAstruct.Mets2Plot))
            u=u+1;
        end
        if ~isempty(MTAstruct.Rates2Plot),
            figure(u)
            plot(1:(MTAstruct.nsteps),bmatrix(MTAstruct.Mets2Plot,:))
            ylabel('dx/dt');
            xlabel('time steps')
            legend(MetNames(MTAstruct.Mets2Plot))
            u=u+1;
        end
        if ~isempty(MTAstruct.Rxns2Plot),
            if isfield(MTAstruct.GEM,'rxnNames'),
                RxnNames=MTAstruct.GEM.rxnNames;
            elseif ~isfield(MTAstruct.GEM,'rxnNames'),
                for j=1:size(MTAstruct.GEM.S,2),
                    RxnNames{j}=['v',num2str(j)];
                end
            end
            figure(u)
            plot(1:(MTAstruct.nsteps),Vmatrix(MTAstruct.Rxns2Plot,:))
            if isfield(MTAstruct,'Vt') && ~isempty(intersect(MTAstruct.Rxns2Plot,MTAstruct.Vt(:,2))),
                DataPoints=find(ismember(MTAstruct.Vt(:,2),MTAstruct.Rxns2Plot));
                Xaxis=zeros(length(DataPoints),1);
                for k=1:length(DataPoints),
                    Xaxis(k)=MTAstruct.Vt(DataPoints(k),3);
                end
                hold on 
                plot(Xaxis,MTAstruct.Vt(DataPoints,1),'*k','MarkerSize',7)
            end
            ylabel('Flux value');
            xlabel('time steps')
            legend(RxnNames(MTAstruct.Rxns2Plot))
        end
        
%**************************************************************************        
%                       MTA Alternative optima analysis  
%************************************************************************** 
  
     elseif MTAstruct.FVAanalysis=='T' || MTAstruct.AOSanalysis=='T',
        Amat=[Amat;sparse(ObjVector')]; %sum(vobj) = Za || vobj(tf) = Zb || sum(sum(deltas)) = Zc
        tic
        %AOS analysis
        if MTAstruct.AOSanalysis=='T',    
            oldAsize=size(Amat,2);
            B=sparse(Mets*(MTAstruct.nsteps+1),oldAsize);B(:,((Rxns+Mets)*MTAstruct.nsteps+1):((Rxns+Mets)*MTAstruct.nsteps+Mets*(MTAstruct.nsteps+1)))=speye(Mets*(MTAstruct.nsteps+1));
            B=[B,[-speye(Mets*(MTAstruct.nsteps+1)),speye(Mets*(MTAstruct.nsteps+1))]];
            Aaos=[[Amat,sparse(size(Amat,1),2*Mets*(MTAstruct.nsteps+1))];B]; %x - eplus + eminus = xrand
            lbaos=[lbvec;zeros(2*Mets*(MTAstruct.nsteps+1),1)];
            ubaos=[ubvec;MTAstruct.Emax*ones(2*Mets*(MTAstruct.nsteps+1),1)]; 
            ObjVector=[zeros(oldAsize,1);ones(2*Mets*(MTAstruct.nsteps+1),1)];       
            
            %prepare gurobi structure
            m.obj=ObjVector;
            m.A=Aaos;
            m.modelsense='min';
            m.lb=lbaos;
            m.ub=ubaos;
            m.sense=[sensevec;'=';repmat('=',Mets*(MTAstruct.nsteps+1),1)];
            params.outputflag=0;
%             params.BarConvTol=0;
            params.Threads=8;
            params.Crossover=1;
%             if MTAstruct.Crossover=='F',
%                params.crossover=0;
%             end
            params.method=2; %0=primal simplex;1=dual simplex;2=barrier method
            wbar = waitbar(0,'Sampling Alternative Optima...');
            
            %main Loop
            n=1;itecounter=1;
            while n<=MTAstruct.AOSnsamples && itecounter<2*MTAstruct.AOSnsamples,
               waitbar(n/MTAstruct.AOSnsamples)
               
               %generate random concentration vector
               Xrand=totubmets.*rand(Mets*(MTAstruct.nsteps+1),1);              
               m.rhs=[rhsvec;MTAstruct.Zoptimum;Xrand];
               Xmatrix=zeros(Mets,MTAstruct.nsteps+1);
               
                try
                   %solve LP
                   gur=gurobi(m,params);
                   
                   %get concentrations
                   p=1;
                   for k=1:MTAstruct.nsteps+1,
                      Xmatrix(:,k)=gur.x(MTAstruct.nsteps*(Rxns+Mets)+p:MTAstruct.nsteps*(Rxns+Mets)+k*Mets);
                      p=p+Mets;
                   end
                   Xmatrix(abs(Xmatrix)<1e-8)=0;
                   Xsamples.(['n',num2str(n)])=Xmatrix;  
                   
                   %quality control
                   Qc=Aaos*gur.x;
                   QC(n,1)=Qc(length(sensevec)+1);
                   n=n+1;
                   itecounter=itecounter+1;
                end
            end
            close(wbar);
            Sol.AOSQC=QC;
            
            if MTAstruct.AOSvarDist=='T',
                %quantify variability in alternative optimal concentration
                %trajectories(obtain total differences between trajectories)
                for k=1:Mets,
                    l=1;
                    for i=1:(n-2),
                        for j=(i+1):(n-1),
                            D(k,l)=sum((Xsamples.(['n',num2str(i)])(k,:)-Xsamples.(['n',num2str(j)])(k,:)).^2); 
                            l=l+1;
                        end
                    end 
                end
                if MTAstruct.saveDaos=='T',
                   Sol.Daos=D;
                end
                Sol.sumDaos=sum(D,2);
                Sol.meanDaos=mean(D,2);
                Sol.medianDaos=median(D,2);
            end
            
            %quantify variability in alternative optimal concentration
            %trajectories(obtain sum of partial entropies per each
            %time-step and across metabolites) and obtain average
            %trajectory across all samples for each metabolite
            meanX=zeros(Mets,MTAstruct.nsteps+1);sumPartEntro=zeros(Mets,1);
            for i=1:Mets,
                metsamples=[];
               for j=1:(n-1),
                  metsamples=[metsamples;Xsamples.(['n',num2str(j)])(i,:)];
               end
               meanX(i,:)=mean(metsamples);
               sumPartEntro(i,1)=sum(DistEntropy(metsamples',[],'T'));
            end
            
            if MTAstruct.saveAOsamples=='T',
                Sol.Xsamples=Xsamples;
            end
            Sol.meanX=meanX;
            Sol.sumPartEntro=sumPartEntro;
        end
        
        %FVA analysis
        if  MTAstruct.FVAanalysis=='T',
            %prepare gurobi structure
            m.A=Amat;
            m.lb=lbvec;
            m.ub=ubvec;
            m.sense=[sensevec;'='];
            m.rhs=[rhsvec;MTAstruct.Zoptimum];
            params.OutputFlag=0;
            params.method=2; %0=primal simplex;1=dual simplex;2=barrier method
%             params.TimeLimit=60;
%             if MTAstruct.Crossover=='F',
%                params.crossover=0;
%             end
            params.crossover=1;
            params.Threads=8;
%             params.BarConvTol=0;
                 
            %main loop: concentration ranges
            if isfield(MTAstruct,'Mets2EvalFVA'),
                Xminmatrix=zeros(length(MTAstruct.Mets2EvalFVA),MTAstruct.nsteps+1);
                Xmaxmatrix=Xminmatrix;
                for i=1:length(MTAstruct.Mets2EvalFVA),
                    wbar=waitbar(0,['Calculating bounds for ',MTAstruct.GEM.metNames{MTAstruct.Mets2EvalFVA(i)},'...']); 
                    for j=0:MTAstruct.nsteps,
                        waitbar((j+1)/(MTAstruct.nsteps+1))
                        %solve LPs
                        ObjVector=zeros(size(Amat,2),1);
                        ObjVector(Rxns*MTAstruct.nsteps+Mets*(MTAstruct.nsteps)+MTAstruct.Mets2EvalFVA(i)+j*Mets)=1;
                        m.obj=ObjVector;
                        m.modelsense='min';
                         try
                          gurmin=gurobi(m,params);
                         catch
                             gurmin.objval=nan;
                         end
                        m.modelsense='max';
                         try
                          gurmax=gurobi(m,params);
                         catch
                             gurmax.objval=nan;
                         end

                        %get min and max concentrations
                        try
                            Xminmatrix(i,j+1)=gurmin.objval;
                        catch
                            Xminmatrix(i,j+1)=nan;
                        end
                        try
                            Xmaxmatrix(i,j+1)=gurmax.objval;
                        catch
                            Xmaxmatrix(i,j+1)=nan;
                        end
                    end  
                    close(wbar);
                end
                Xminmatrix(abs(Xminmatrix)<1e-8)=0;
                Xmaxmatrix(abs(Xmaxmatrix)<1e-8)=0;
            end
            
            %main loop: flux ranges
            if isfield(MTAstruct,'Rxns2EvalFVA'),
                Vminmatrix=zeros(length(MTAstruct.Rxns2EvalFVA),MTAstruct.nsteps);
                Vmaxmatrix=Vminmatrix;
                wbar=waitbar(0,'Calculating AO flux range...'); 
                waitcount=1;
                for i=1:length(MTAstruct.Rxns2EvalFVA),
                    waitbar(waitcount/((MTAstruct.nsteps+1)*length(MTAstruct.Rxns2EvalFVA)))
                    for j=0:(MTAstruct.nsteps-1),
                        %solve LPs
                        waitcount=i*(j+1);
                        ObjVector=zeros(size(Amat,2),1);
                        ObjVector(MTAstruct.Rxns2EvalFVA(i)+j*Rxns)=1;
                        m.obj=ObjVector;
                        m.modelsense='min';
                        try
                          gurmin=gurobi(m,params);
                        catch
                            gurmin.objval=nan;
                        end
                        m.modelsense='max';
                        try
                          gurmax=gurobi(m,params);
                        catch
                            gurmax.objval=nan;
                        end

                        %get min and max concentrations (Needs to be
                        %adjusted to account splitted reversible reactions 
                        %when ConcFluxConstraint='T')
                         Vminmatrix(i,j+1)=gurmin.objval;
                         Vmaxmatrix(i,j+1)=gurmax.objval;
                    end                       
                end
                close(wbar);
                Vminmatrix(abs(Vminmatrix)<1e-8)=0;
                Vmaxmatrix(abs(Vmaxmatrix)<1e-8)=0;
            end
            
            if MTAstruct.saveFVAranges=='T',
                if isfield(MTAstruct,'Mets2EvalFVA'),
                   Sol.XminMat=Xminmatrix;
                   Sol.XmaxMat=Xmaxmatrix;
                end
                if isfield(MTAstruct,'Rxns2EvalFVA'),
                   Sol.VminMat=Vminmatrix;
                   Sol.VmaxMat=Vmaxmatrix;
                end
            end
        end 
        Sol.AOTime=toc;
        %plot figures
        if exist('MetPercentile','var');
            MTAstruct.Mets2Plot=find(Sol.meanDaos<quantile(Sol.meanDaos,MetPercentile));
%             MTAstruct.Mets2Plot=find(Sol.sumPartEntro<quantile(Sol.sumPartEntro,MetPercentile));
            Sol.Mets2Plot=[MTAstruct.GEM.metNames(MTAstruct.Mets2Plot),num2cell(MTAstruct.Mets2Plot)];
        end
        if ~isempty(MTAstruct.Mets2Plot),
            if isfield(MTAstruct.GEM,'metNames'),
                MetNames=MTAstruct.GEM.metNames;
            elseif ~isfield(MTAstruct.GEM,'metNames'),
                for i=1:size(MTAstruct.GEM.S,1),
                    MetNames{i}=['x',num2str(i)];
                end
            end
            
            for u=1:length(MTAstruct.Mets2Plot),
                fig(u)=figure('Color','w','Visible','off');
                if MTAstruct.AOSanalysis=='T',
                    for v=1:(n-1),
                        hold on
                        plot(1:(MTAstruct.nsteps+1),Xsamples.(['n',num2str(v)])(MTAstruct.Mets2Plot(u),:),'Color',[0.75,0.75,0.75])
                        %plot data points
                        if isfield(MTAstruct,'Xt') && ismember(MTAstruct.Mets2Plot(u),MTAstruct.Xt(:,2)),
                            DataPoints=find(ismember(MTAstruct.Xt(:,2),MTAstruct.Mets2Plot(u)));
                            Xaxis=zeros(length(DataPoints),1);
                            for k=1:length(DataPoints),
                                Xaxis(k)=MTAstruct.Xt(DataPoints(k),3);
                            end
                            hold on 
                            plot(Xaxis,MTAstruct.Xt(DataPoints,1),'*k','MarkerSize',7)
                        end
                    end
                    hold on
                    plot(1:(MTAstruct.nsteps+1),meanX(MTAstruct.Mets2Plot(u),:),'k','LineWidth',2)
                    ylabel('Concentration','FontSize',16);
                    xlabel('time steps','FontSize',16);
                    title(MetNames(MTAstruct.Mets2Plot(u)),'FontSize',16)
                    get(gca,'XTick');
                    set(gca,'FontSize',16)
                end
                if  MTAstruct.FVAanalysis=='T',
                    T=find(MTAstruct.Mets2EvalFVA==MTAstruct.Mets2Plot(u));
                    if ~isempty(T),
                        hold on
                        plot(1:(MTAstruct.nsteps+1),Xminmatrix(T,:),'--b','LineWidth',2)
                        hold on
                        plot(1:(MTAstruct.nsteps+1),Xmaxmatrix(T,:),'--r','LineWidth',2)
                        ylabel('Concentration','FontSize',16);
                        xlabel('time steps','FontSize',16);
                        get(gca,'XTick');
                        set(gca,'FontSize',16)
                        title(MetNames(MTAstruct.Mets2EvalFVA(T)),'FontSize',16)
                    end
                end
                print(fig(u),[MTAstruct.FigDir,'\',MTAstruct.GEM.metNames{MTAstruct.Mets2Plot(u)}],'-dtiff')
            end
        end
    end
end