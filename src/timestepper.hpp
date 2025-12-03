#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void doStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void doStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

    class ImprovedEuler : public TimeStepper
  {
    Vector<> m_vecf;      // f-Wert
    Vector<> m_ytilde;    // Zwischenwert ỹ
  public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
      : TimeStepper(rhs),
        m_vecf(rhs->dimF()),
        m_ytilde(rhs->dimX())
    { }

    void doStep(double tau, VectorView<double> y) override
    {
      // 1) f(y_n)
      this->m_rhs->evaluate(y, m_vecf);

      // 2) ỹ = y_n + (tau/2) * f(y_n)
      m_ytilde = y;
      m_ytilde += 0.5 * tau * m_vecf;

      // 3) f(ỹ)
      this->m_rhs->evaluate(m_ytilde, m_vecf);

      // 4) y_{n+1} = y_n + tau * f(ỹ)
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void doStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };

    class CrankNicolson : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
    std::shared_ptr<ConstantFunction> m_fold;   // f(y_n)

  public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
      : TimeStepper(rhs),
        m_tau(std::make_shared<Parameter>(0.0))
    {
      // y_n und f(y_n) als konstante Funktionen
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      m_fold = std::make_shared<ConstantFunction>(rhs->dimF());

      // y_new
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());

      // Funktionsgleichung:
      // G(y_new) = y_new - y_old - 0.5*tau*( f(y_old) + f(y_new) )
      //
      // 0.5*tau * f(y_new)
      //auto half_tau = 0.5 * m_tau;
      auto term_new = 0.5 * (m_tau * m_rhs);

      // 0.5*tau * f(y_old)  -> hier benutzen wir m_fold
      auto term_old = 0.5 * (m_tau * m_fold);

      m_equ = ynew - m_yold - term_new - term_old;
    }

    void doStep(double tau, VectorView<double> y) override
    {
      // y_n speichern
      m_yold->set(y);

      // f(y_n) in m_fold speichern
      Vector<> tmp_f(m_rhs->dimF());
      m_rhs->evaluate(y, tmp_f);
      m_fold->set(tmp_f);

      m_tau->set(tau);

      // Newton löst G(y_{n+1}) = 0, Ausgangswert = y_n
      NewtonSolver(m_equ, y);
    }
  };



  

}


#endif
