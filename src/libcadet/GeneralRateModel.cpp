// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2014: Eric von Lieres¹, Joel Andersson,
//                         Andreas Puettmann¹, Sebastian Schnittert¹,
//                         Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <iomanip>
#include <limits>
#include <typeinfo>

#include "GeneralRateModel.hpp"
#include "TimeIntegrator.hpp"
#include "JacobianData.hpp"
#include "SchurSolver.hpp"
#include "WenoScheme.hpp"
#include "ParticleDiscretization.hpp"
#include "AdsorptionModel.hpp"

namespace cadet {

//int GeneralRateModel::count = 0;
//int GeneralRateModel::countS = 0;

GeneralRateModel::GeneralRateModel(SimulatorPImpl& sim) :
    ChromatographyModel(sim, GENERAL_RATE_MODEL),
    _psim(sim),
    _jac (sim.getSchurSolver().getJacobianData()),
    _ws  (sim.getWenoScheme()),
    _c_in(_cc.ncomp(), 0.0)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double inf = std::numeric_limits<double>::infinity();

    // Scalar parameters
    addParam(Parameter<active> (COL_LENGTH,     e2s(COL_LENGTH),     -1, 0.0, 0.0, 0.0, CADET_STRICT, inf, CADET_STRICT)); ///todo re-check loose / strict
    addParam(Parameter<active> (COL_POROSITY,   e2s(COL_POROSITY),   -1, 0.0, 0.0, 0.0, CADET_LOOSE,  1.0, CADET_LOOSE));
    addParam(Parameter<active> (COL_DISPERSION, e2s(COL_DISPERSION), -1, 0.0, 0.0, 0.0, CADET_LOOSE,  inf, CADET_STRICT));
    addParam(Parameter<active> (VELOCITY,       e2s(VELOCITY),       -1, 0.0, 0.0, 0.0, CADET_STRICT, inf, CADET_STRICT));

    addParam(Parameter<active> (PAR_RADIUS,     e2s(PAR_RADIUS),     -1, 0.0, 0.0, 0.0, CADET_STRICT, inf, CADET_STRICT));
    addParam(Parameter<active> (PAR_POROSITY,   e2s(PAR_POROSITY),   -1, 0.0, 0.0, 0.0, CADET_STRICT, 1.0, CADET_STRICT));

    // Vectorial parameters
    for (int comp = 0; comp < _cc.ncomp(); ++comp)
    {
        addParam(Parameter<active> (FILM_DIFFUSION,    e2s(FILM_DIFFUSION),    comp, 0.0, 0.0, 0.0, CADET_LOOSE, inf, CADET_STRICT));
        addParam(Parameter<active> (PAR_DIFFUSION,     e2s(PAR_DIFFUSION),     comp, 0.0, 0.0, 0.0, CADET_LOOSE, inf, CADET_STRICT));
        addParam(Parameter<active> (PAR_SURFDIFFUSION, e2s(PAR_SURFDIFFUSION), comp, 0.0, 0.0, 0.0, CADET_LOOSE, inf, CADET_STRICT));
    }

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


GeneralRateModel::~GeneralRateModel()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

///todo Checkl destructors and other stuff... rule of three, for all classes!

///todo why are all parameters active? shouldn't it be possible to use double params for computations?

int GeneralRateModel::residualDae(double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res, void* userData)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    _timerResDae.start();

    //============================================
    double* y    = NV_DATA_S(NV_y);
    double* yDot = NV_DATA_S(NV_yDot);
    double* res  = NV_DATA_S(NV_res);
    //============================================

    _timerResDaePar.start();

#ifndef VERIFY_ANALYTICAL_JAC
    // Residual evaluation including analytical jacobian computation
    if (_psim.getTimeIntegrator().useAnalyticJacobian())
#endif
        residualColumnParticle<double, double, double, true> (t, y, yDot, res);


#ifndef VERIFY_ANALYTICAL_JAC
    // Residual evaluation including AD jacobian computation
    else
#endif
    {
        // reinitialize actives
        for (int i = 0; i < _cc.neq(); ++i)
        {
            // Copy content of state vector (NV) to AD state vector, init directional derivatives with zero
            _ti.getYAd(i).setValue(y[i]);
            _ti.getResAd(i).setValue(0.0);

            for (int dir = 0; dir < _ti.getJacAdDirs(); ++dir)
                _ti.getResAd(i).setADValue(dir, 0.0);
        }

        residualColumnParticle<active, active, double, false> (t, _ti.getYAd(), yDot, _ti.getResAd());

        // copy residual values from AD residual vector to residual vector (NV)
        for (int i = 0; i < _cc.neq(); ++i)
            res[i] = _ti.getResAd(i).getValue();

#ifdef VERIFY_ANALYTICAL_JAC
        // Comparison of analytical an AD jacobian implementation
        _psim.getSchurSolver().getJacobianData().compareWithAd(_ti.getResAd(), _ti.getDiagDir());
#endif

        // Copy Jacobian entries from resAd to band Jacobian data structures...
        _psim.getSchurSolver().getJacobianData().setFromAd(_ti.getResAd(), _ti.getDiagDir());
    }

    _timerResDaePar.stop();

    // now take care of the boundaries (without any sensitivity computation)...
    residualBoundaries<double, double> (y, yDot, res);

    // mark jacobian factorization for next call of schurSolve
    _psim.getTimeIntegrator().setFactorizeJac(true);


    _timerResDae.stop();

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}




int GeneralRateModel::residualSens(int ns, double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res,
        N_Vector* NV_yS, N_Vector* NV_ySDot, N_Vector* NV_resS,
        void* userData, N_Vector NV_tmp1, N_Vector NV_tmp2, N_Vector NV_tmp3)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    _timerResSens.start();

    NV_tmp1 = _ti.getNvTemp1(); // stores result of (dF / dy) * s
    NV_tmp2 = _ti.getNvTemp2(); // stores result of (dF / dyDot) * sDot

    double* tmp1 = NV_DATA_S(NV_tmp1);
    double* tmp2 = NV_DATA_S(NV_tmp2);

    double* resS;

    double* y    = NV_DATA_S(NV_y);
    double* yDot = NV_DATA_S(NV_yDot);

    _timerResSensPar.start();

    residualColumnParticle<double, active, active, false> (t, y, yDot, _ti.getResAd());

    _timerResSensPar.stop();

    //=============================================================

    residualBoundaries<active, active> (y, yDot, _ti.getResAd());

    for (int param = 0; param < _ti.getNSensParams(); param++)
    {
        //=============================================================================
        // Directional derivative (dF / dy) * s
        dFdy_times_s(NV_yS[param], NV_tmp1);

        // Directional derivative (dF / dyDot) * sDot
        dFdyDot_times_sDot(NV_ySDot[param], NV_tmp2);

        resS = NV_DATA_S(NV_resS[param]);

        _timerResSensPar.start();

        // Complete sens residual is the sum:
        #pragma omp parallel for
        for (int i = 0; i < _cc.neq(); i++)
            resS[i] = tmp1[i] + tmp2[i] + _ti.getResAd(i).getADValue(_ti.getJacAdDirs() + param);

        _timerResSensPar.stop();
        //=============================================================================
    }

    // do not factorize Jacobians at next call of schurSolve
    _psim.getTimeIntegrator().setFactorizeJac(false);

    _timerResSens.stop();
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}



void GeneralRateModel::calcIC(const double t)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Compute column residual without time derivative (->nullptr) and store this in ydot.
    // This means ydot will contain the sum of the convective and the dispersive flux term.
    residualColumn<double, double, double, false>(t, _ti.getY(), nullptr, _ti.getYDot());

    double invBetaC              = 1.0 / getValue<double>(COL_POROSITY) - 1.0;
    double surfaceToVolumeRatio  = 3.0 / getValue<double>(PAR_RADIUS);

    // ydot shall fulfill ydot = -column_residual - film_flux
    // ==> negate and subtract boundary contribution
    for (int eqc = 0; eqc < _cc.neq_col(); ++eqc)
    {
        _ti.getYDot()[eqc] *= -1.0;
        _ti.getYDot()[eqc] -= invBetaC * surfaceToVolumeRatio * _cc.offsetBnd(_ti.getNvY())[eqc]; // <- boundary contribution!
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::calcICSens(const double t)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Quick return if we have no sensitivities activated
    if (_ti.getNSensParams() < 1) return;

    // call DAE residual to compute Jacobian dF/dy
    residualDae(t, _ti.getNvY(), _ti.getNvYDot(), _ti.getNvTemp1(), nullptr);

    // call residuals for sensitivities
    residualColumnParticle<double, active, active, false>(t, _ti.getY(), _ti.getYDot(), _ti.getResAd());
    residualBoundaries<active, active>(_ti.getY(), _ti.getYDot(), _ti.getResAd());

    for (int param = 0; param < _ti.getNSensParams(); ++param)
    {
        // compute Jacobian times sensitivity -> Js
        dFdy_times_s(_ti.getNvYS(param), _ti.getNvTemp1());

        // For the column, the matrix in front of sdot, dF/dydot, is simply the identity.
        // ==>   sdot = -J*s - dF/dp    if parameter has direct influence
        // ==>   sdot = -J*s            if parameter has no direct influence
        for (int eqc = 0; eqc < _cc.neq_col(); eqc++)
            _ti.getYSDot(param)[eqc] = - _ti.getTemp1()[eqc] - _ti.getResAd(eqc).getADValue(_ti.getJacAdDirs() + param);
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::specialSetup()
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    _dc_indp = std::vector<std::vector<double> >(_cc.ncomp(), std::vector<double>(_ti.getMaxSensInletParams(), 0.0));

    assembleOffdiagJac();

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}



// this residual function now only handles column and particles
// boundaries are treated differently elsewhere
template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
int GeneralRateModel::residualColumnParticle(const double t, const StateType* y, const double* yDot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    #pragma omp parallel for
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        if (pblk == -1)
            residualColumn<StateType, ResidType, ParamType, wantJac> (t, y, yDot, res);
        else
        {
            const StateType* par_y    = y    + _cc.neq_col() + pblk * _cc.neq_par(); ///todo use offset inline functions
            const double*    par_yDot = yDot + _cc.neq_col() + pblk * _cc.neq_par();
            ResidType*       par_res  = res  + _cc.neq_col() + pblk * _cc.neq_par();

            residualParticle<StateType, ResidType, ParamType, wantJac> (t, pblk, par_y, par_yDot, par_res);
        }
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}


template <>
void GeneralRateModel::setInletParamDerivatives<double>(std::vector<double>& concInlet)
{
    // If ParamType is a double, no derivatives are needed
}

template <>
void GeneralRateModel::setInletParamDerivatives<active>(std::vector<active>& concInlet)
{
    if (_ti.getNSensInletParams() > 0)
    {
        int localParamIndex = 0;
        for (int param = 0; param < _ti.getMaxSensInletParams(); ++param)
        {
            if (_ti.getInletParamIsSensitive(param)) // Found a sensitive inlet parameter
            {
                for (int comp = 0; comp < _cc.ncomp(); ++comp)
                {
                    concInlet.at(comp).setADValue(_ti.getJacAdDirs() + _ti.getNSensModelParams() + localParamIndex,
                            _dc_indp.at(comp).at(param));
                }
                localParamIndex++;
            }
        }
    }
}

//       t     input       current time
//       y     input       [array] pointer to first column element of state vector
//      yp     input       [array] pointer to first column element of derivative state vector
//       p     input       [scalar] pointer to parameter data structure containing all model parameters (sensisitve as well as non-sensitive)
//     res    output       [array] calculated residual values
//  csdata                 ChromsimData structure
//
template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
int GeneralRateModel::residualColumn(const double t, const StateType* y, const double* yDot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    ParamType u   = getValue<ParamType>(VELOCITY);
    ParamType h   = getValue<ParamType>(COL_LENGTH) / double(_cc.ncol());
    ParamType h2  = h * h;
    ParamType d_c = getValue<ParamType>(COL_DISPERSION);

    // WENO help variables
    const int& WK   = _cc.max_wk();           // Max. WENO-order
    StateType* work = new StateType[3 * WK];  // required by "wenoReconstruct"
    StateType* v    = new StateType[2 * WK];  // Stencil space
    double*    Dvm  = new double[2 * WK - 1]; // Derivatives of vm

    StateType* w = v + WK;  // Stencil pointer bend for readability - shortcut for v[WK + x]

    double* jac;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Clear content of the inlet concentrations derivative matrix
    for (int i = 0; i < _cc.ncomp(); ++i)
        _dc_indp.at(i).assign(_ti.getMaxSensInletParams(), 0.0);

    // Evaluate the user-specified function for the inlet concentration
    _ti.inletConcentration(t, _ti.getSection(t), _c_in, _ti.getInletParamIsSensitive(), _dc_indp);

    // Copy content of double vector "_c_in" into a vector of ParamTypes
    std::vector<ParamType> concInlet(_c_in.begin(),_c_in.end());

    // In case ParamType is active, copy derivatives to the respective active type, else do nothing
    setInletParamDerivatives<ParamType>(concInlet);
    /////////////////////////////////////////////////////////////////////////////////////////////////


    for (int comp = 0; comp < _cc.ncomp(); ++comp)
    {
        // Add time derivative
        if (yDot != nullptr)
            for (int col = 0; col < _cc.ncol(); ++col)
                _cc.colC<ResidType> (res, col, comp) = _cc.colC<double> (yDot, col, comp);
        else
            for (int col = 0; col < _cc.ncol(); ++col)
                _cc.colC<ResidType> (res, col, comp) = 0.0;

        // Stencil:
        w[2] = _cc.colC<StateType> (y, 2, comp);
        w[1] = _cc.colC<StateType> (y, 1, comp);
        w[0] = _cc.colC<StateType> (y, 0, comp);

        StateType vm = 0.0; // reconstructed value

        for (int i = 0; i < 2 * WK - 1; ++i)
            Dvm[i] = 0.0;

        int wk = 0;  // Current WENO-order

        for (int col = 0; col < _cc.ncol(); ++col)
        {
            // Jacobian entries
            if (wantJac)
            {
                jac = _jac.getJacC(comp,col);
                // Initialize Jacobian to zero
                for (int i = -WK; i < WK; ++i)  // -3 <= i <= 2
                    jac[WK + i] = 0.0;
            }

            // Add dispersion
            if (col < _cc.ncol() - 1) // right side
            {
                _cc.colC<ResidType> (res, col, comp) -= d_c / h2 * (w[ 1] - w[0]);
                // Jacobian entries
                if (wantJac)
                {
                    jac[WK]     += ParamType(d_c / h2);
                    jac[WK + 1] -= ParamType(d_c / h2);
                }
            }

            if (col > 0) // left side
            {
                _cc.colC<ResidType> (res, col, comp) -= d_c / h2 * (w[-1] - w[0]);
                // Jacobian entries
                if (wantJac)
                {
                    jac[WK]     += ParamType(d_c / h2);
                    jac[WK - 1] -= ParamType(d_c / h2);
                }
            }

            //-----------------------------------------------------------------
            // Add convection through this cell's left face
            if (col == 0)
                // for the first cell, the concentration at its left face
                // is determined by the inflow concentration
                _cc.colC<ResidType> (res, col, comp) -= u / h * concInlet.at(comp);
            else
                // ... for all other cells use the reconstructed value ...
                // Remember that vm still contains the reconstructed value
                // of the previous cell's *right* face,
                // which is identical to this cell's left face!
                _cc.colC<ResidType> (res, col, comp) -= u / h * vm;
            // Jacobian entries
            if (wantJac)
            {
                for (int i = 0; i < 2 * wk - 1; ++i)
                    jac[WK - wk + i] -= ParamType(u / h * Dvm[i]);
            }
            //-----------------------------------------------------------------


            // Boundaries
            int bnd = 0;
            switch (_ws.getBoundaryModel())
            {
            case 0: // Lower WENO order
                // This very statement selects the max. weno order for the current column cell
                // wk = min(maxWKleft, maxWKright)
                wk = min(min(col + 1, _ws.getWenoOrder()), min(_cc.ncol() - col, _ws.getWenoOrder()));
                break;

            case 1: // Zero weights
                wk = _ws.getWenoOrder();
                if (col < wk - 1)
                    bnd = -(wk - 1 - col);
                else if (col > _cc.ncol() - wk)
                    bnd = _cc.ncol() - col;
                break;

            case 2: // Zero weights for p != 0
                if (col == 0)
                    wk = 1;
                else
                {
                    wk = _ws.getWenoOrder();
                    if (col < wk - 1)
                        bnd = -(wk - 1 - col);
                    else if (col > _cc.ncol() - wk)
                        bnd = _cc.ncol() - col;
                }
                break;

            case 3: // Large ghost points
                wk = _ws.getWenoOrder();
                if (col == 0)
                {
                    w[-1] = 1e20;
                    w[-2] = 1e50;
                }
                else if (col == _cc.ncol() - 2)
                    w[2] = 1e20;
                else if (col == _cc.ncol() - 1)
                    w[2] = 1e50;
                break;

            default:
                std::ostringstream ss;
                ss << "GeneralRateModel::residualColumn(): Wrong boundary model specified - accepted values [0-3]: " << wk;
                throw CadetException(ss.str());
                break;
            }

            // Reconstruct concentration on this cell's right face
            _ws.wenoReconstruct<StateType, wantJac>(wk, comp, bnd, w - wk + 1, &vm, Dvm, work);

            // Right side
            _cc.colC<ResidType> (res, col, comp) += u / h * vm;
            // Jacobian entries
            if (wantJac)
                for (int i = 0; i < 2 * wk - 1; ++i)
                    jac[WK - wk + i + 1] += ParamType(u / h * Dvm[i]);

            // Update stencil
            w[-3] = w[-2];
            w[-2] = w[-1];
            w[-1] = w[ 0];
            w[ 0] = w[ 1];
            w[ 1] = w[ 2];
            w[ 2] = _cc.colC<StateType> (y, col + 3, comp);
        }
    }

    delete [] work;
    delete [] Dvm;
    delete [] v;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}




template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
int GeneralRateModel::residualParticle(const double t, const int pblk, const StateType* y, const double* ydot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Simulation parameters
    ParamType radius            = getValue<ParamType> (PAR_RADIUS);
    ParamType inv_beta_p        = 1.0 / getValue<ParamType> (PAR_POROSITY) - 1.0;
    std::vector<ParamType> d_p  = getValueForAllComp<ParamType> (PAR_DIFFUSION);
    std::vector<ParamType> d_s  = getValueForAllComp<ParamType> (PAR_SURFDIFFUSION);

    double z = 1.0 / double(_cc.npblk()) * (0.5 + pblk) ; // Current z coordinate in the column - needed in externally dependent adsorption kinetic
    ParamType dr;

    ParamType pTypeInfo; // This is only to determine the residual function that is called

    // go to the first diagonal element of Jacobian (ie. skip undefined values)
    double* jac = _jac.getJacP(pblk) + _jac.ku_par();

    // Loop over particle cells
    for (int par = 0; par < _cc.npar(); ++par)
    {
        // Geometry
        ParamType outer_area_per_volume = _pd.getOuterSurfAreaPerVolume(par) / radius;
        ParamType inner_area_per_volume = _pd.getInnerSurfAreaPerVolume(par) / radius;

        // Mobile phase
        for (int comp = 0; comp < _cc.ncomp(); ++comp, ++res, ++y, ++ydot, jac += _jac.ld_jp())
        {
            // Add time derivatives
            *res = *ydot + inv_beta_p * ydot[_cc.ncomp()];

            // set dres / dc_i and dres / dq_i = 0
            if (wantJac)
            {
                jac[0]           = 0.0;
                jac[_cc.ncomp()] = 0.0;
            }

            // Add flow through outer surface
            if (par != 0)
            {
                // difference between two cell-centers
                dr = (_pd.getParCellCoords().at(par - 1) - _pd.getParCellCoords().at(par)) * radius;

                // Gradient approximation
                ResidType grad_c = (y[-2 * _cc.ncomp()] - y[0]) / dr;
                ResidType grad_q = (y[-_cc.ncomp()] - y[_cc.ncomp()]) / dr;

                // Molecular diffusion contribution
                *res -= outer_area_per_volume * d_p.at(comp) * grad_c;

                // Surface diffusion contribution
                *res -= outer_area_per_volume * d_s.at(comp) * inv_beta_p * grad_q;

                if (wantJac)
                {
                        // Jacobian entrys w.r.t. this cell's concentrations
                        jac[0]           += ParamType(outer_area_per_volume * d_p.at(comp) / dr);                      // dres / dc_p,i^(p,j)
                        jac[_cc.ncomp()] += ParamType(outer_area_per_volume * inv_beta_p * d_s.at(comp) / dr);         // dres / dq_i^(p,j)

                        // Jacobian entrys w.r.t. neighboring cell's concentrations
                        jac[-2 * _cc.ncomp()] = ParamType(-outer_area_per_volume * d_p.at(comp) / dr);                 // dres / dc_p,i^(p,j-1)
                        jac[-_cc.ncomp()]     = ParamType(-outer_area_per_volume * inv_beta_p * d_s.at(comp) / dr);    // dres / dq_i^(p,j-1)
                }
//                else
//                    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
            }

            // Add flow through inner surface
            if (par != _cc.npar() - 1)
            {
                // difference between two cell-centers
                dr = (_pd.getParCellCoords().at(par) - _pd.getParCellCoords().at(par + 1)) * radius;

                // Gradient approximation
                ResidType grad_c = (y[0] - y[2 * _cc.ncomp()]) / dr;
                ResidType grad_q = (y[_cc.ncomp()] - y[3 * _cc.ncomp()]) / dr;

                // Molecular diffusion contribution
                *res += inner_area_per_volume * d_p.at(comp) * grad_c;

                // Surface diffusion contribution
                *res += inner_area_per_volume * d_s.at(comp) * inv_beta_p * grad_q;

                if (wantJac)
                {
                        // Jacobian entrys w.r.t. this cell's concentrations
                        jac[0]           += ParamType(inner_area_per_volume * d_p.at(comp) / dr);                    // dres / dc_p,i^(p,j)
                        jac[_cc.ncomp()] += ParamType(inner_area_per_volume * inv_beta_p * d_s.at(comp) / dr);       // dres / dq_i^(p,j)

                        // Jacobian entrys w.r.t. neighboring cell's concentrations
                        jac[2 * _cc.ncomp()] = ParamType(-inner_area_per_volume * d_p.at(comp) / dr);                // dres / dc_p,i^(p,j+1)
                        jac[3 * _cc.ncomp()] = ParamType(-inner_area_per_volume * inv_beta_p * d_s.at(comp) / dr);   // dres / dq_i^(p,j+1)
                }
//                still a bad hack with the operator double() !!! ///todo remove this bad hack
//                else
//                    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
            }
        }

        // Bound phase
        for (int comp = 0; comp < _cc.ncomp(); ++comp, ++res, ++y, ++ydot, jac += _jac.ld_jp())
        {
            // === call the generic isotherm residual function =============
            _am.evaluateResidual(t, z, comp, y, res, &pTypeInfo);
            if (wantJac)
            {
                if (typeid(y) == typeid(const double*))
                    _am.setJacobian(t, z, comp, (const double*)(y), jac);
                else
                    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
            }

            // === Add time derivative for adsorption models ===============
            if (_am.isDifferential(comp)) *res += *ydot;
        }
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}



template <typename ResidType, typename ParamType>
int GeneralRateModel::residualBoundaries(const double* y, const double* yDot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    ParamType invBetaC          = 1.0 / getValue<ParamType>(COL_POROSITY) - 1.0;
    ParamType epsP              = getValue<ParamType>(PAR_POROSITY);
    ParamType radius            = getValue<ParamType>(PAR_RADIUS);
    std::vector<ParamType> kf   = getValueForAllComp<ParamType>(FILM_DIFFUSION);
    std::vector<ParamType> dp   = getValueForAllComp<ParamType>(PAR_DIFFUSION);

    ParamType surfaceToVolumeRatio = 3.0 / radius;
    ParamType outerAreaPerVolume   = _pd.getOuterSurfAreaPerVolume(0) / radius;

    ParamType jacCB_val = invBetaC * surfaceToVolumeRatio;
    ParamType jacPB_val = -outerAreaPerVolume / epsP;

    ParamType* kf_FV = new ParamType[_cc.ncomp()];  // kf for finite volumes

    const double relOuterShellHalfRadius = 0.5 * _pd.getCellSize(0);
    for (int comp = 0; comp < _cc.ncomp(); ++comp)
        kf_FV[comp] = 1.0 / (radius * relOuterShellHalfRadius / epsP / dp.at(comp) + 1.0 / kf.at(comp));

    int eq;  // Index for the current equation we work on

    ResidType* res_col   = res;
    ResidType* res_par   = res + _cc.neq_col();
    ResidType* res_bound = res + _cc.neq_col() + _cc.npblk() * _cc.neq_par();

    const double* y_col   = y;
    const double* y_par   = y + _cc.neq_col();
    const double* y_bound = y + _cc.neq_col() + _cc.npblk() * _cc.neq_par();

    //========================================================================
    // J_b part
    //========================================================================
    for (int comp = 0; comp < _cc.ncomp(); ++comp)
        for (int bnd = 0; bnd < _cc.nbnd(); ++bnd)
        {
            eq = bnd + comp * _cc.nbnd();
            res_bound[eq] = y_bound[eq];
        }
    //========================================================================


    //========================================================================
    // J_c,b part
    //========================================================================
    for (eq = 0; eq < _cc.neq_col(); ++eq)
        res_col[eq] += jacCB_val * y_bound[eq];
    //========================================================================


    //========================================================================
    // J_b,c part
    //========================================================================
    for (int bnd = 0; bnd < _cc.nbnd(); ++bnd)
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = bnd + comp * _cc.nbnd();
            res_bound[eq] += -kf_FV[comp] * y_col[eq];
        }
    //========================================================================


    //========================================================================
    // J_p,b part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = pblk + comp * _cc.npblk();
            res_par[pblk * _cc.neq_par() + comp] += jacPB_val * y_bound[eq];
        }
    //========================================================================


    //========================================================================
    // J_b,p part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); pblk++)
        for (int comp = 0; comp < _cc.ncomp(); comp++)
        {
            eq = pblk + comp * _cc.npblk();
            res_bound[eq] += kf_FV[comp] * y_par[comp + pblk * _cc.neq_par()];
        }
    //========================================================================

    delete [] kf_FV;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}




void GeneralRateModel::dFdy_times_s(N_Vector NV_s, N_Vector NV_ret)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    N_VScale(1.0, NV_s, NV_ret);

    lapackInt_t     n;
    lapackInt_t     kl;
    lapackInt_t     ku;
    lapackInt_t     lda;

    double      alpha = 1.0;
    double      beta = 0.0;
    lapackInt_t inc = 1;

    char trans[] = "Trans";

    _timerResSensPar.start();

    #pragma omp parallel for private(n, kl, ku, lda)
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        //==========================================================================
        // Column data
        //==========================================================================
        if (pblk == -1)
        {
            n   = _cc.neq_col();
            kl  = _jac.kl_col();
            ku  = _jac.ku_col();
            lda = _jac.ld_jc();

            DGBMV(trans, &n, &n, &kl, &ku, &alpha, _jac.getJacC(),
                    &lda, _cc.offsetCol(NV_s), &inc, &beta,
                    _cc.offsetCol(NV_ret), &inc);

            _jac.sparseMV(_jac.getJacCB(), _jac.numel_jcb(), 1.0, _cc.offsetBnd(NV_s), _cc.offsetCol(NV_ret));
        }
        //==========================================================================
        // Particle data
        //==========================================================================
        else
        {
            n   = _cc.neq_par();
            kl  = _jac.kl_par();
            ku  = _jac.ku_par();
            lda = _jac.ld_jp();

            DGBMV(trans, &n, &n, &kl, &ku, &alpha, _jac.getJacP(pblk),
                    &lda, _cc.offsetPar(NV_s, pblk), &inc, &beta,
                    _cc.offsetPar(NV_ret, pblk), &inc);

            _jac.sparseMV(_jac.getJacPB(pblk), _jac.numel_jpb(), 1.0, _cc.offsetBnd(NV_s), _cc.offsetPar(NV_ret, pblk));
        }
        //==========================================================================
    }
    _timerResSensPar.stop();

    //==========================================================================
    // Boundary data
    //==========================================================================
    _jac.sparseMV(_jac.getJacBC(), _jac.numel_jbc(), 1.0, _cc.offsetCol(NV_s), _cc.offsetBnd(NV_ret));

    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
        _jac.sparseMV(_jac.getJacBP(pblk), _jac.numel_jbp(), 1.0, _cc.offsetPar(NV_s, pblk), _cc.offsetBnd(NV_ret));
    //==========================================================================

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}



void GeneralRateModel::dFdyDot_times_sDot(N_Vector NV_sDot, N_Vector NV_ret)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double invBetaP = 1.0 / getValue<double> (PAR_POROSITY) - 1.0;

    _timerResSensPar.start();

    #pragma omp parallel for
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        if (pblk == -1) // column
        {
            double* sDot = _cc.offsetCol(NV_sDot);
            double* ret  = _cc.offsetCol(NV_ret);

            for (int i = 0; i < _cc.neq_col(); ++i) // loop over column equations
                ret[i] = sDot[i];
        }
        else // particle
        {
            double* sDot = _cc.offsetPar(NV_sDot, pblk);
            double* ret  = _cc.offsetPar(NV_ret, pblk);

            for (int par = 0; par < _cc.npar(); ++par)
                for (int comp = 0; comp < _cc.ncomp(); ++comp)
                    _cc.parC<double> (ret, par, comp) = _cc.parC<double> (sDot, par, comp)
                    + invBetaP * _cc.parQ<double> (sDot, par, comp);

            for (int par = 0; par < _cc.npar(); ++par)
                for (int comp = 0; comp < _cc.ncomp(); ++comp)
                    if (_am.isDifferential(comp)) ///todo might be faster like this: parQ = _am.isDifferential(comp) * parQ ...
                        _cc.parQ<double> (ret, par, comp) = _cc.parQ<double> (sDot, par, comp);
                    else
                        _cc.parQ<double> (ret, par, comp) = 0.0;
        }
    }
    _timerResSensPar.stop();

    double* dFdyDot = _cc.offsetBnd(NV_ret);
    for (int eqb = 0; eqb < _cc.neq_bnd(); ++eqb) // loop over boundary equations
        dFdyDot[eqb] = 0.0;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::assembleOffdiagJac() throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    SparseMatrixElement* jcb;
    SparseMatrixElement* jbc;
    SparseMatrixElement* jpb;
    SparseMatrixElement* jbp;

    double invBetaC          = 1.0 / getValue<double>(COL_POROSITY) - 1.0;
    double epsP              = getValue<double>(PAR_POROSITY);
    double radius            = getValue<double>(PAR_RADIUS);
    std::vector<double> kf   = getValueForAllComp<double>(FILM_DIFFUSION);
    std::vector<double> dp   = getValueForAllComp<double>(PAR_DIFFUSION);


    double surfaceToVolumeRatio = 3.0 / radius;
    double outerAreaPerVolume   = _pd.getOuterSurfAreaPerVolume(0) / radius;

    log::emit<Debug2>() << CURRENT_FUNCTION << ": radius = " << radius << log::endl;
    log::emit<Debug2>() << CURRENT_FUNCTION << ": _pd.getOuterSurfAreaPerVolume(0) = " << _pd.getOuterSurfAreaPerVolume(0) << log::endl;

    log::emit<Debug2>() << CURRENT_FUNCTION << ": surfaceToVolumeRatio = " << surfaceToVolumeRatio << log::endl;
    log::emit<Debug2>() << CURRENT_FUNCTION << ": outerAreaPerVolume = " << outerAreaPerVolume << log::endl;

    double jacCB_val = invBetaC * surfaceToVolumeRatio;
    double jacPB_val = -outerAreaPerVolume / epsP;

    log::emit<Debug2>() << CURRENT_FUNCTION << ": jacCB_val = " << jacCB_val << log::endl;
    log::emit<Debug2>() << CURRENT_FUNCTION << ": jacPB_val = " << jacPB_val << log::endl;

    double* kf_FV = new double[_cc.ncomp()];  // kf for finite volumes ///todo compute and store kf_FV in the specialSetup routine!!!

    for (int comp = 0; comp < _cc.ncomp(); ++comp)
    {
        kf_FV[comp] = 1.0 / (0.5 * radius / _cc.npar() / epsP / dp.at(comp) + 1.0 / kf.at(comp));
        log::emit<Debug2>() << CURRENT_FUNCTION << ": kf_FV[" << comp << "] = " << kf_FV[comp] << log::endl;
    }

    int eq;  // Index for the current equation we work on

    //========================================================================
    // J_c,b part
    //========================================================================
    jcb = _jac.getJacCB();
    for (eq = 0; eq < _cc.neq_col(); ++eq)
    {
        jcb->setElement(eq, eq, jacCB_val);
        jcb++;
    }
    //========================================================================


    //========================================================================
    // J_b,c part
    //========================================================================
    jbc = _jac.getJacBC();
    for (int col = 0; col < _cc.ncol(); ++col)
    {
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = col + comp * _cc.ncol();
            jbc->setElement(eq, eq, -kf_FV[comp]);
            jbc++;
        }
    }
    //========================================================================



    //========================================================================
    // J_p,b part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
    {
        jpb = _jac.getJacPB(pblk);
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = pblk + comp * _cc.ncol();
            jpb->setElement(comp, eq, jacPB_val);
            jpb++;
        }
    }
    //========================================================================


    //========================================================================
    // J_b,p part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
    {
        jbp = _jac.getJacBP(pblk);
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = pblk + comp * _cc.ncol();
            jbp->setElement(eq, comp, kf_FV[comp]);
            jbp++;
        }
    }
    //========================================================================

    delete [] kf_FV;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


} // namespace cadet