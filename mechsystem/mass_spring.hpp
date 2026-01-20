#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

class DistanceConstraint
{
public:
  double length;
  std::array<Connector,2> connectors;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  std::vector<DistanceConstraint> m_constraints;
  Vec<D> m_gravity=0.0;
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  size_t addConstraint (DistanceConstraint c)
  {
    m_constraints.push_back (c);
    return m_constraints.size()-1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }
  auto & constraints() { return m_constraints; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D*mss.masses().size() + mss.constraints().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size() + mss.constraints().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;

    size_t nmasses = mss.masses().size();
    auto xmat = x.asMatrix(nmasses, D);
    auto fmat = f.asMatrix(nmasses, D);

    // Extract Lagrange multipliers from state vector
    VectorView<double> lambda(mss.constraints().size(), x.data() + D*nmasses);

    // Gravity and spring forces
    for (size_t i = 0; i < nmasses; i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    for (auto spring : mss.springs())
      {
        auto [c1,c2] = spring.connectors;
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        double force = spring.stiffness * (norm(p1-p2)-spring.length);
        Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);
        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force*dir12;
      }

    // Add constraint forces: F_i += lambda_j * grad_x g_j
    for (size_t j = 0; j < mss.constraints().size(); j++)
      {
        auto constraint = mss.constraints()[j];
        auto [c1, c2] = constraint.connectors;
        
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        // Gradient of constraint g = |p1 - p2|^2 - l^2
        // dg/dp1 = 2(p1 - p2), dg/dp2 = -2(p1 - p2)
        Vec<D> grad_g = 2.0 * (p1 - p2);
        
        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += lambda(j) * grad_g;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= lambda(j) * grad_g;
      }

    // Divide by mass for acceleration
    for (size_t i = 0; i < nmasses; i++)
      fmat.row(i) *= 1.0/mss.masses()[i].mass;

    // Constraint equations: g(x) = |p1-p2|^2 - l^2 = 0
    VectorView<double> g(mss.constraints().size(), f.data() + D*nmasses);
    for (size_t j = 0; j < mss.constraints().size(); j++)
      {
        auto constraint = mss.constraints()[j];
        auto [c1, c2] = constraint.connectors;
        
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        double dist = norm(p1 - p2);
        g(j) = dist*dist - constraint.length * constraint.length;
      }
  }
  
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    
    size_t nmasses = mss.masses().size();
    auto xmat = x.asMatrix(nmasses, D);
    VectorView<double> lambda(mss.constraints().size(), x.data() + D*nmasses);
    
    // Build Jacobian block by block
    // Block 1: d(acceleration)/d(positions)
    for (auto spring : mss.springs())
      {
        auto [c1, c2] = spring.connectors;
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);
        
        Vec<D> d = p1 - p2;
        double r = norm(d);
        Vec<D> dir = (1.0/r) * d;
        
        // Spring force: F = -k(r - l) * dir
        // dF/d(p1) = -k[dir ⊗ dir + (r-l)/r * (I - dir ⊗ dir)]
        // dF/d(p2) = -dF/d(p1)
        double coeff1 = spring.stiffness / r;
        double coeff2 = spring.stiffness * (r - spring.length) / (r * r);
        
        // Hessian contribution for mass equations
        if (c1.type == Connector::MASS)
          {
            size_t i1 = c1.nr;
            double mass1 = mss.masses()[i1].mass;
            // dF1/d(p1)
            for (int a = 0; a < D; a++)
              for (int b = 0; b < D; b++)
                {
                  double val = coeff1 * dir(a) * dir(b) - coeff2 * (a == b ? 1.0 : 0.0);
                  df(i1*D + a, i1*D + b) += val / mass1;
                }
            // dF1/d(p2)
            if (c2.type == Connector::MASS)
              {
                size_t i2 = c2.nr;
                double mass2 = mss.masses()[i2].mass;
                for (int a = 0; a < D; a++)
                  for (int b = 0; b < D; b++)
                    {
                      double val = -coeff1 * dir(a) * dir(b) + coeff2 * (a == b ? 1.0 : 0.0);
                      df(i1*D + a, i2*D + b) += val / mass1;
                    }
              }
          }
        
        if (c2.type == Connector::MASS)
          {
            size_t i2 = c2.nr;
            double mass2 = mss.masses()[i2].mass;
            // dF2/d(p2)
            for (int a = 0; a < D; a++)
              for (int b = 0; b < D; b++)
                {
                  double val = coeff1 * dir(a) * dir(b) - coeff2 * (a == b ? 1.0 : 0.0);
                  df(i2*D + a, i2*D + b) += val / mass2;
                }
            // dF2/d(p1)
            if (c1.type == Connector::MASS)
              {
                size_t i1 = c1.nr;
                double mass1 = mss.masses()[i1].mass;
                for (int a = 0; a < D; a++)
                  for (int b = 0; b < D; b++)
                    {
                      double val = -coeff1 * dir(a) * dir(b) + coeff2 * (a == b ? 1.0 : 0.0);
                      df(i2*D + a, i1*D + b) += val / mass2;
                    }
              }
          }
      }
    
    // Block 2: Constraint force Hessian contributions: d(lambda_j * grad_g_j)/d(x)
    for (size_t j = 0; j < mss.constraints().size(); j++)
      {
        auto constraint = mss.constraints()[j];
        auto [c1, c2] = constraint.connectors;
        
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);
        
        // d^2 g / d(p1)^2 = 2*I
        if (c1.type == Connector::MASS)
          {
            size_t i1 = c1.nr;
            double mass1 = mss.masses()[i1].mass;
            for (int a = 0; a < D; a++)
              for (int b = 0; b < D; b++)
                df(i1*D + a, i1*D + b) += 2.0 * lambda(j) * (a == b ? 1.0 : 0.0) / mass1;
          }
        
        // d^2 g / d(p1)d(p2) = -2*I
        if (c1.type == Connector::MASS && c2.type == Connector::MASS)
          {
            size_t i1 = c1.nr;
            size_t i2 = c2.nr;
            double mass1 = mss.masses()[i1].mass;
            for (int a = 0; a < D; a++)
              for (int b = 0; b < D; b++)
                df(i1*D + a, i2*D + b) -= 2.0 * lambda(j) * (a == b ? 1.0 : 0.0) / mass1;
          }
        
        // d^2 g / d(p2)^2 = 2*I
        if (c2.type == Connector::MASS)
          {
            size_t i2 = c2.nr;
            double mass2 = mss.masses()[i2].mass;
            for (int a = 0; a < D; a++)
              for (int b = 0; b < D; b++)
                df(i2*D + a, i2*D + b) += 2.0 * lambda(j) * (a == b ? 1.0 : 0.0) / mass2;
          }
        
        // d^2 g / d(p2)d(p1) = -2*I
        if (c1.type == Connector::MASS && c2.type == Connector::MASS)
          {
            size_t i1 = c1.nr;
            size_t i2 = c2.nr;
            double mass2 = mss.masses()[i2].mass;
            for (int a = 0; a < D; a++)
              for (int b = 0; b < D; b++)
                df(i2*D + a, i1*D + b) -= 2.0 * lambda(j) * (a == b ? 1.0 : 0.0) / mass2;
          }
      }
    
    // Block 3: d(acceleration)/d(lambda) = grad_g transposed
    for (size_t j = 0; j < mss.constraints().size(); j++)
      {
        auto constraint = mss.constraints()[j];
        auto [c1, c2] = constraint.connectors;
        
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);
        
        Vec<D> grad_g = 2.0 * (p1 - p2);
        
        if (c1.type == Connector::MASS)
          {
            size_t i1 = c1.nr;
            double mass1 = mss.masses()[i1].mass;
            for (int a = 0; a < D; a++)
              df(i1*D + a, D*nmasses + j) = grad_g(a) / mass1;
          }
        if (c2.type == Connector::MASS)
          {
            size_t i2 = c2.nr;
            double mass2 = mss.masses()[i2].mass;
            for (int a = 0; a < D; a++)
              df(i2*D + a, D*nmasses + j) = -grad_g(a) / mass2;
          }
      }
    
    // Block 4: d(constraint)/d(positions) = jacobian of constraint
    for (size_t j = 0; j < mss.constraints().size(); j++)
      {
        auto constraint = mss.constraints()[j];
        auto [c1, c2] = constraint.connectors;
        
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);
        
        Vec<D> grad_g = 2.0 * (p1 - p2);
        
        if (c1.type == Connector::MASS)
          {
            size_t i1 = c1.nr;
            for (int a = 0; a < D; a++)
              df(D*nmasses + j, i1*D + a) = grad_g(a);
          }
        if (c2.type == Connector::MASS)
          {
            size_t i2 = c2.nr;
            for (int a = 0; a < D; a++)
              df(D*nmasses + j, i2*D + a) = -grad_g(a);
          }
      }
    
    // Block 5: d(constraint)/d(lambda) = 0 (already initialized to zero)
  }
  
};

#endif
