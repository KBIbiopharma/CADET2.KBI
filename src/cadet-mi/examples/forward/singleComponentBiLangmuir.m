function singleComponentBiLangmuir()
%SINGLECOMPONENTBILANGMUIR Simple single component Bi-Langmuir simulation
%
% The second bound state is realized by adding a pseudo component to the
% system. The interstitial and mobile phase of these pseudo components are
% not used (i.e. only transport equations are solved). The solid phase of
% the pseudo components serves as second bound state.
%
% The layout of the components is as follows:
%  [Comp1a, ... , CompNa,  Comp1b, ..., CompNb]
%   |<- "Real Comp" ->|   |<- "Pseudo Comp" ->|
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    model = ModelGRM();
    
    % General
    model.nComponents = 2; % nComponents has to be twice the amount of "real components"

    % Initial conditions
	model.initialMobileConcentration = [0.0 0.0];  % Second (artificial) component has to be 0.0
    model.initialSolidConcentration = [0.0 0.0];   % Second (artificial) component has to be 0.0
    
    % Adsorption
    model.kineticBindingModel = true;
    model.bindingModel = MultiComponentBiLangmuirBinding();
    model.bindingParameters.MCBL_KA1         = [1.14];
    model.bindingParameters.MCBL_KD1         = [0.002];
    model.bindingParameters.MCBL_QMAX1       = [4.88];
    model.bindingParameters.MCBL_KA2         = [0.3];
    model.bindingParameters.MCBL_KD2         = [0.004];
    model.bindingParameters.MCBL_QMAX2       = [1.12];
    
    % Transport
    model.dispersionColumn          = 5.75e-8;
    model.filmDiffusion             = [6.9e-6 6.9e-6];      % Same for second (artificial) component
    model.diffusionParticle         = [6.07e-11 6.07e-11];  % Same for second (artificial) component
    model.diffusionParticleSurface  = [0.0 0.0];            % Same for second (artificial) component
    model.interstitialVelocity      = 5.75e-4;

    % Geometry
    model.columnLength        = 0.014;
    model.particleRadius      = 4.5e-5;
    model.porosityColumn      = 0.37;
    model.porosityParticle    = 0.75;
    
    % Inlet
    model.nInletSections = 1;
    model.sectionTimes = [0.0 14000];
    model.sectionContinuity = [];
    
    model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
    model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
    model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
    model.sectionCubic          = zeros(model.nComponents, model.nInletSections);

    % Sec 1
    model.sectionConstant(1,1)  = 7.14e-3;  % component 1

    % Discretization
    disc = DiscretizationGRM();
    disc.nCellsColumn = 16;
    disc.nCellsParticle = 4;
    
    % Simulator
    sim = Simulator(model, disc);
    sim.solutionTimes = linspace(0, 14000, 2001);
   
    % Run
    result = sim.simulate();
            
    plot(result.solution.time, result.solution.outlet(:,1));
    legend('Comp 1');
    grid on;
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson¹,
%                         Andreas Puettmann¹, Sebastian Schnittert¹,
%                         Samuel Leweke¹
%                                      
%    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
