// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright (c) 2008-2012: Eric von Lieres¹, Joel Andersson¹,
//                           Andreas Puettmann¹, Sebastian Schnittert¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0
//  which accompanies this distribution, and is available at
//  http://www.gnu.org/licenses/gpl.html
//  ---------------------------------------------------------------------------
//  Author    Sebastian Schnittert <schnittert@gmail.com>
//  Version   $Id: AdsorptionModel_THM.hpp 519 2012-10-08 12:55:41Z schnittert $
// =============================================================================

#ifndef ADSORPTIONMODEL_EXT_SMAPH_HPP_
#define ADSORPTIONMODEL_EXT_SMAPH_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Steric Mass Action adsorption model with externally dependent parameters
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_EXTSMAPH : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_EXTSMAPH(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, EXTERNAL_STERIC_MASS_ACTION_PH), 
        _lambda(EXTSMAPH_LAMBDA,     e2s(EXTSMAPH_LAMBDA),     -1, -1, 0.0, 0.0, -std::numeric_limits<double>::infinity(), true, std::numeric_limits<double>::infinity(), true)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        _kA.reserve(_cc.ncomp());
        _kAE.reserve(_cc.ncomp());
        _kAEE.reserve(_cc.ncomp());

        _kD.reserve(_cc.ncomp());
        _kDE.reserve(_cc.ncomp());
        _kDEE.reserve(_cc.ncomp());

        _nu.reserve(_cc.ncomp());
        _nuP.reserve(_cc.ncomp());
        _nuPP.reserve(_cc.ncomp());

        _sigma.reserve(_cc.ncomp());
        _sigmaP.reserve(_cc.ncomp());
        _sigmaPP.reserve(_cc.ncomp());

        addParam(_lambda);

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _kA.push_back(Parameter<active> (EXTSMAPH_KA,         e2s(EXTSMAPH_KA),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kA[comp]);
            _kAE.push_back(Parameter<active> (EXTSMAPH_KA_E,       e2s(EXTSMAPH_KA_E),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kAE[comp]);
            _kAEE.push_back(Parameter<active> (EXTSMAPH_KA_EE,      e2s(EXTSMAPH_KA_EE),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kAEE[comp]);

            _kD.push_back(Parameter<active> (EXTSMAPH_KD,         e2s(EXTSMAPH_KD),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kD[comp]);
            _kDE.push_back(Parameter<active> (EXTSMAPH_KD_E,       e2s(EXTSMAPH_KD_E),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDE[comp]);
            _kDEE.push_back(Parameter<active> (EXTSMAPH_KD_EE,      e2s(EXTSMAPH_KD_EE),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDEE[comp]);

            _nu.push_back(Parameter<active> (EXTSMAPH_NU,         e2s(EXTSMAPH_NU),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_nu[comp]);
            _nuP.push_back(Parameter<active> (EXTSMAPH_NU_P,       e2s(EXTSMAPH_NU_P),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_nuP[comp]);
            _nuPP.push_back(Parameter<active> (EXTSMAPH_NU_PP,      e2s(EXTSMAPH_NU_PP),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_nuPP[comp]);

            _sigma.push_back(Parameter<active> (EXTSMAPH_SIGMA,      e2s(EXTSMAPH_SIGMA),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_sigma[comp]);
            _sigmaP.push_back(Parameter<active> (EXTSMAPH_SIGMA_P,    e2s(EXTSMAPH_SIGMA_P),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_sigmaP[comp]);
            _sigmaPP.push_back(Parameter<active> (EXTSMAPH_SIGMA_PP,   e2s(EXTSMAPH_SIGMA_PP),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_sigmaPP[comp]);
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_EXTSMAPH()
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Public members

    virtual void setIsKinetic(bool isKinetic)
    {
        _isKinetic = isKinetic;
        for (int comp = 1; comp < _cc.ncomp(); ++comp)  // start only at comp 1, since salt-eq. is always non-differential
            _isDifferential.at(comp) = isKinetic;
    }

    virtual void evaluateResidual(const double t, const double z, const int comp, const active * q, active * res, const active * p) const
        { evaluateResidual<active, active, active>(t, z, comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const active * q, active * res, const double * p) const
        { evaluateResidual<active, active, double>(t, z, comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, active * res, const active * p) const
        { evaluateResidual<double, active, active>(t, z, comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, double * res, const double * p) const
        { evaluateResidual<double, double, double>(t, z, comp, q, res); }

    virtual void setJacobian(const double t, const double z, const int comp, const double* q, double* jac) const throw (CadetException)
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double              ka        = _kA[comp].getValue<double>();
        const double              ka_E      = _kAE[comp].getValue<double>();
        const double              ka_EE     = _kAEE[comp].getValue<double>();
        const double              kd        = _kD[comp].getValue<double>();
        const double              kd_E      = _kDE[comp].getValue<double>();
        const double              kd_EE     = _kDEE[comp].getValue<double>();

        // Temperature
        double temp;
        _externalBase->externalProfile(z, t, &temp);

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        if (comp == 0)  // Salt component
        {
            jac[0] = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                jac[j] = _nu[j].getValue<double>() + temp * (_nuP[j].getValue<double>() + temp * (_nuPP[j].getValue<double>() ));
            }
        }
        else  // Protein component
        {
            // Salt concentrations in liquid and solid phase
            const double c0 = c[-comp];
            const double q0 = q[-comp];

            double q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalSigma = _sigma[j].getValue<double>() + temp * (_sigmaP[j].getValue<double>() + temp * (_sigmaPP[j].getValue<double>() ));
                q0_bar -= finalSigma * q[-comp + j];
            }

            const double finalNu = _nu[comp].getValue<double>() + temp * (_nuP[comp].getValue<double>() + temp * (_nuPP[comp].getValue<double>() ));
            const double c0_pow_nu     = pow(c0, finalNu);
            const double q0_bar_pow_nu = pow(q0_bar, finalNu);

            //const double finalKa = ka + temp * (ka_T + temp * (ka_TT));
            //const double finalKd = kd + temp * (kd_T + temp * (kd_TT));
            const double finalKa = pow(double(10.0), ka + temp * (ka_E + temp * (ka_EE)) );
            const double finalKd = pow(double(10.0), kd + temp * (kd_E + temp * (kd_EE)) );

            // Jacobian
            jac[-_cc.ncomp() - comp] = finalKd * (*q) * finalNu * c0_pow_nu / c0;                      // dres_i / dc0
            jac[-_cc.ncomp()] = -finalKa * q0_bar_pow_nu;                                              // dres_i / dci
            jac[-comp] = -finalKa * (*c) * finalNu * q0_bar_pow_nu / q0_bar;                           // dres_i / dq0

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalSigma = _sigma[j].getValue<double>() + temp * (_sigmaP[j].getValue<double>() + temp * (_sigmaPP[j].getValue<double>() ));
                jac[-comp + j] = -finalKa * (*c) * finalNu * q0_bar_pow_nu / q0_bar * (-finalSigma);  // dres_i / dqj
            }

            jac[0] += finalKd * c0_pow_nu;                                                               // dres_i / dqi
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const double t, const double z, const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const ParamType              ka        = _kA[comp].getValue<ParamType>();
        const ParamType              ka_E      = _kAE[comp].getValue<ParamType>();
        const ParamType              ka_EE     = _kAEE[comp].getValue<ParamType>();
        const ParamType              kd        = _kD[comp].getValue<ParamType>();
        const ParamType              kd_E      = _kDE[comp].getValue<ParamType>();
        const ParamType              kd_EE     = _kDEE[comp].getValue<ParamType>();
        const ParamType              lambda    = _lambda.getValue<ParamType>();

        // Temperature
        double temp;
        _externalBase->externalProfile(z, t, &temp);

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        if (comp == 0)
        {
            // Salt component
            *res = *q - ( lambda );

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                ResidType finalNu = _nu[j].getValue<double>() + temp * (_nuP[j].getValue<double>() + temp * (_nuPP[j].getValue<double>() ));
                *res += finalNu * q[j];
            }
        }
        else
        {
            // Protein 
            // Salt concentrations in liquid and solid phase
            const StateType c0 = c[-comp];
            const StateType q0 = q[-comp];

            ResidType q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                ResidType finalSigma = _sigma[j].getValue<double>() + temp * (_sigmaP[j].getValue<double>() + temp * (_sigmaPP[j].getValue<double>() ));
                q0_bar -= finalSigma * q[-comp + j];
            }

            ResidType finalNu = _nu[comp].getValue<double>() + temp * (_nuP[comp].getValue<double>() + temp * (_nuPP[comp].getValue<double>() ));
            ResidType c0_pow_nu = pow(c0, finalNu);
            ResidType q0_bar_pow_nu = pow(q0_bar, finalNu);

            //ResidType finalKa = ka + temp * (ka_T + temp * (ka_TT));
            //ResidType finalKd = kd + temp * (kd_T + temp * (kd_TT));
            ResidType finalKa = pow(ResidType(10.0), ka + temp * (ka_E + temp * (ka_EE)) );
            ResidType finalKd = pow(ResidType(10.0), kd + temp * (kd_E + temp * (kd_EE)) );

            // Residual
            *res = finalKd * (*q) * c0_pow_nu - finalKa * (*c) * q0_bar_pow_nu;
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    std::vector<Parameter<active>>  _kA;
    std::vector<Parameter<active>>  _kAE;
    std::vector<Parameter<active>>  _kAEE;

    std::vector<Parameter<active>>  _kD;
    std::vector<Parameter<active>>  _kDE;
    std::vector<Parameter<active>>  _kDEE;

    std::vector<Parameter<active>>  _nu;
    std::vector<Parameter<active>>  _nuP;
    std::vector<Parameter<active>>  _nuPP;

    std::vector<Parameter<active>>  _sigma;
    std::vector<Parameter<active>>  _sigmaP;
    std::vector<Parameter<active>>  _sigmaPP;

    Parameter<active> _lambda;
};

} // namespace cadet

#endif // ADSORPTIONMODEL_EXT_SMAPH_HPP_
