// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson¹,
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

#ifndef ADSORPTIONMODEL_MPM_HPP_
#define ADSORPTIONMODEL_MPM_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Mobile Phyase Modulators adsorption model
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_MPM : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_MPM(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, MOBILE_PHASE_MODULATORS)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        _kA.reserve(_cc.ncomp());
        _kD.reserve(_cc.ncomp());
        _qMax.reserve(_cc.ncomp());
        _beta.reserve(_cc.ncomp());
        _gamma.reserve(_cc.ncomp());

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _kA.push_back(Parameter<active> (MPM_KA,    e2s(MPM_KA),    comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kA[comp]);
            _kD.push_back(Parameter<active> (MPM_KD,    e2s(MPM_KD),    comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kD[comp]);
            _qMax.push_back(Parameter<active> (MPM_QMAX,  e2s(MPM_QMAX),  comp, -1, 0.0, 0.0, 0.0, true,  inf, true));
            addParam(_qMax[comp]);
            _beta.push_back(Parameter<active> (MPM_BETA,  e2s(MPM_BETA),  comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_beta[comp]);
            _gamma.push_back(Parameter<active> (MPM_GAMMA, e2s(MPM_GAMMA), comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_gamma[comp]);
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_MPM()
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Public members

    virtual void setIsKinetic(bool isKinetic)
    {
        _isKinetic = isKinetic;
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
            _isDifferential.at(comp) = isKinetic;
    }

    virtual void evaluateResidual(const double t, const double z, const int comp, const active * q, active * res, const active * p) const
        { evaluateResidual<active, active, active>(comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const active * q, active * res, const double * p) const
        { evaluateResidual<active, active, double>(comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, active * res, const active * p) const
        { evaluateResidual<double, active, active>(comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, double * res, const double * p) const
        { evaluateResidual<double, double, double>(comp, q, res); }


    virtual void setJacobian(const double t, const double z, const int comp, const double* q, double* jac) const throw (CadetException)
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double              ka    = _kA[comp].getValue<double>();
        const double              kd    = _kD[comp].getValue<double>();
        const double              beta  = _beta[comp].getValue<double>();
        const double              gamma = _gamma[comp].getValue<double>();

        // Only protein components
        if (comp > 0)
        {
            // Liquid phase concentration
            const double* c = q - _cc.ncomp();

            // Liquid phase salt concentration
            const double c0 = c[-comp];

            const double ka_mpm = ka * exp(gamma * c0);
            const double kd_mpm = kd * pow(c0, beta);

            double qsum = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                qsum -= q[-comp + j] / _qMax[j].getValue<double>();
            }

            // Jacobian
            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[-comp + j] = ka_mpm * (*c) * _qMax[comp].getValue<double>() / _qMax[j].getValue<double>();              // dresi/dqj

            jac[0]            += kd_mpm;                                                // dresi/dqi
            jac[-_cc.ncomp()]  = -ka_mpm * _qMax[comp].getValue<double>() * qsum;                        // dresi/dci

            jac[-_cc.ncomp() - comp] = -ka_mpm * (*c) * _qMax[comp].getValue<double>() * qsum * gamma
                    + kd * beta * *q * pow(c0, beta - 1);                               // dresi/dc0
        }
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

private:
    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const ParamType              ka    = _kA[comp].getValue<ParamType>();
        const ParamType              kd    = _kD[comp].getValue<ParamType>();
        const ParamType              beta  = _beta[comp].getValue<ParamType>();
        const ParamType              gamma = _gamma[comp].getValue<ParamType>();

        // Liquid phase concentration
        const StateType* c = q - _cc.ncomp();

        if (comp == 0)
        {
            // Salt component
            *res = (ResidType) 0.0;
        }
        else
        {
            // Protein components
            
            // Liquid phase salt concentration
            const StateType c0 = c[-comp];
            
            ResidType ka_mpm = ka * exp(gamma * c0);
            ResidType kd_mpm = kd * pow(c0, beta);

            ResidType qsum = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                qsum -= q[-comp + j] / _qMax[j].getValue<ParamType>();
            }

            // Residual
            *res = - (ka_mpm * (*c) * _qMax[comp].getValue<ParamType>() * qsum - kd_mpm * *q);
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    std::vector<Parameter<active>>  _kA;
    std::vector<Parameter<active>>  _kD;
    std::vector<Parameter<active>>  _qMax;
    std::vector<Parameter<active>>  _beta;
    std::vector<Parameter<active>>  _gamma;    
};

} // namespace cadet

#endif // ADSORPTIONMODEL_MPM_HPP_
