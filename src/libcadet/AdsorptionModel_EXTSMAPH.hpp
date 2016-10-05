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
        _kAT.reserve(_cc.ncomp());
        _kATT.reserve(_cc.ncomp());

        _kD.reserve(_cc.ncomp());
        _kDT.reserve(_cc.ncomp());
        _kDTT.reserve(_cc.ncomp());

        _nu.reserve(_cc.ncomp());
        _nuT.reserve(_cc.ncomp());
        _nuTT.reserve(_cc.ncomp());

        _sigma.reserve(_cc.ncomp());
        _sigmaT.reserve(_cc.ncomp());
        _sigmaTT.reserve(_cc.ncomp());

        addParam(_lambda);

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _kA.push_back(Parameter<active> (EXTSMAPH_KA,         e2s(EXTSMAPH_KA),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kA[comp]);
            _kAT.push_back(Parameter<active> (EXTSMAPH_KA_T,       e2s(EXTSMAPH_KA_T),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kAT[comp]);
            _kATT.push_back(Parameter<active> (EXTSMAPH_KA_TT,      e2s(EXTSMAPH_KA_TT),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kATT[comp]);

            _kD.push_back(Parameter<active> (EXTSMAPH_KD,         e2s(EXTSMAPH_KD),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kD[comp]);
            _kDT.push_back(Parameter<active> (EXTSMAPH_KD_T,       e2s(EXTSMAPH_KD_T),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDT[comp]);
            _kDTT.push_back(Parameter<active> (EXTSMAPH_KD_TT,      e2s(EXTSMAPH_KD_TT),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDTT[comp]);

            _nu.push_back(Parameter<active> (EXTSMAPH_NU,         e2s(EXTSMAPH_NU),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_nu[comp]);
            _nuT.push_back(Parameter<active> (EXTSMAPH_NU_T,       e2s(EXTSMAPH_NU_T),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_nuT[comp]);
            _nuTT.push_back(Parameter<active> (EXTSMAPH_NU_TT,      e2s(EXTSMAPH_NU_TT),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_nuTT[comp]);

            _sigma.push_back(Parameter<active> (EXTSMAPH_SIGMA,      e2s(EXTSMAPH_SIGMA),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_sigma[comp]);
            _sigmaT.push_back(Parameter<active> (EXTSMAPH_SIGMA_T,    e2s(EXTSMAPH_SIGMA_T),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_sigmaT[comp]);
            _sigmaTT.push_back(Parameter<active> (EXTSMAPH_SIGMA_TT,   e2s(EXTSMAPH_SIGMA_TT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_sigmaTT[comp]);
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
        const double              ka_T      = _kAT[comp].getValue<double>();
        const double              ka_TT     = _kATT[comp].getValue<double>();
        const double              kd        = _kD[comp].getValue<double>();
        const double              kd_T      = _kDT[comp].getValue<double>();
        const double              kd_TT     = _kDTT[comp].getValue<double>();

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
                jac[j] = _nu[j].getValue<double>() + temp * (_nuT[j].getValue<double>() + temp * (_nuTT[j].getValue<double>() ));
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
                const double finalSigma = _sigma[j].getValue<double>() + temp * (_sigmaT[j].getValue<double>() + temp * (_sigmaTT[j].getValue<double>() ));
                q0_bar -= finalSigma * q[-comp + j];
            }

            const double finalNu = _nu[comp].getValue<double>() + temp * (_nuT[comp].getValue<double>() + temp * (_nuTT[comp].getValue<double>() ));
            const double c0_pow_nu     = pow(c0, finalNu);
            const double q0_bar_pow_nu = pow(q0_bar, finalNu);

            const double finalKa = ka + temp * (ka_T + temp * (ka_TT));
            const double finalKd = kd + temp * (kd_T + temp * (kd_TT));

            // Jacobian
            jac[-_cc.ncomp() - comp] = finalKd * (*q) * finalNu * c0_pow_nu / c0;                      // dres_i / dc0
            jac[-_cc.ncomp()] = -finalKa * q0_bar_pow_nu;                                              // dres_i / dci
            jac[-comp] = -finalKa * (*c) * finalNu * q0_bar_pow_nu / q0_bar;                           // dres_i / dq0

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalSigma = _sigma[j].getValue<double>() + temp * (_sigmaT[j].getValue<double>() + temp * (_sigmaTT[j].getValue<double>() ));
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
        const ParamType              ka_T      = _kAT[comp].getValue<ParamType>();
        const ParamType              ka_TT     = _kATT[comp].getValue<ParamType>();
        const ParamType              kd        = _kD[comp].getValue<ParamType>();
        const ParamType              kd_T      = _kDT[comp].getValue<ParamType>();
        const ParamType              kd_TT     = _kDTT[comp].getValue<ParamType>();
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
                ResidType finalNu = _nu[j].getValue<double>() + temp * (_nuT[j].getValue<double>() + temp * (_nuTT[j].getValue<double>() ));
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
                ResidType finalSigma = _sigma[j].getValue<double>() + temp * (_sigmaT[j].getValue<double>() + temp * (_sigmaTT[j].getValue<double>() ));
                q0_bar -= finalSigma * q[-comp + j];
            }

            ResidType finalNu = _nu[comp].getValue<double>() + temp * (_nuT[comp].getValue<double>() + temp * (_nuTT[comp].getValue<double>() ));
            ResidType c0_pow_nu = pow(c0, finalNu);
            ResidType q0_bar_pow_nu = pow(q0_bar, finalNu);

            //ResidType finalKa = ka + temp * (ka_T + temp * (ka_TT));
            //ResidType finalKd = kd + temp * (kd_T + temp * (kd_TT));
            ResidType finalKa = pow(ResidType(10.0), ka + temp * (ka_T + temp * (ka_TT)) );
            ResidType finalKd = pow(ResidType(10.0), kd + temp * (kd_T + temp * (kd_TT)) );

            // Residual
            *res = finalKd * (*q) * c0_pow_nu - finalKa * (*c) * q0_bar_pow_nu;
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    std::vector<Parameter<active>>  _kA;
    std::vector<Parameter<active>>  _kAT;
    std::vector<Parameter<active>>  _kATT;

    std::vector<Parameter<active>>  _kD;
    std::vector<Parameter<active>>  _kDT;
    std::vector<Parameter<active>>  _kDTT;

    std::vector<Parameter<active>>  _nu;
    std::vector<Parameter<active>>  _nuT;
    std::vector<Parameter<active>>  _nuTT;

    std::vector<Parameter<active>>  _sigma;
    std::vector<Parameter<active>>  _sigmaT;
    std::vector<Parameter<active>>  _sigmaTT;

    Parameter<active> _lambda;
};

} // namespace cadet

#endif // ADSORPTIONMODEL_EXT_SMAPH_HPP_
