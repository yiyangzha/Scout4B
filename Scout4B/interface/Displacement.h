#ifndef Scout4B_Displacement_h
#define Scout4B_Displacement_h

#include "Scout4B/Scout4B/interface/KinFitUtils.h"
#include "Scout4B/Scout4B/interface/KinematicFitResult.h"

namespace Scout4B
{
  /* The class is used to compute displacement information for a given
     decay tree stored in KinematicFitResult with respect to a given
     production vertex. Each instance of the class has a name to mark
     what vertex type is used. It also keeps a copy of the vertex and
     the index number of the primary vertex if applicable. */

  class Displacement
  {
  public:
    Displacement(const char *name, const KinematicFitResult &fit,
                 const reco::Vertex &vertex, int pv_index) : fit_(fit), prodVertex_(vertex), pvIndex_(pv_index), name_(name)
    {
      compute_displacement();
    }

    Displacement(const char *name) : name_(name) {}
    double decayLength() const { return decayLength_; }
    double decayLengthErr() const { return decayLengthErr_; }
    double decayLengthSig() const { return decayLengthErr_ > 0 ? decayLength_ / decayLengthErr_ : 0; }
    double distaceOfClosestApproach() const { return distaceOfClosestApproach_; }
    double distaceOfClosestApproachErr() const { return distaceOfClosestApproachErr_; }
    double distaceOfClosestApproachSig() const { return distaceOfClosestApproachSig_; }
    double longitudinalImpactParameter() const { return longitudinalImpactParameter_; }
    double longitudinalImpactParameterErr() const { return longitudinalImpactParameterErr_; }
    double longitudinalImpactParameterSig() const { return longitudinalImpactParameterSig_; }
    double decayTime() const { return decayTime_; }
    double decayTimeError() const { return decayTimeError_; }
    double decayTimeXY() const { return decayTimeXY_; }
    double decayTimeXYError() const { return decayTimeXYError_; }
    double alpha() const { return alpha_; }
    double alphaErr() const { return alphaErr_; }
    double alphaXY() const { return alphaXY_; }
    double alphaXYErr() const { return alphaXYErr_; }
    int pvIndex() const { return pvIndex_; }
    const std::string &name() const { return name_; }
    const reco::Vertex &prodVertex() const { return prodVertex_; }

  private:
    void compute_displacement();
    double decayLength_ = -1.0;
    double decayLengthErr_ = 0.0;
    double distaceOfClosestApproach_ = -1.0;
    double distaceOfClosestApproachErr_ = 0.0;
    double distaceOfClosestApproachSig_ = 0.0;
    double longitudinalImpactParameter_ = 0.0;
    double longitudinalImpactParameterErr_ = 0.0;
    double longitudinalImpactParameterSig_ = 0.0;
    double decayTime_ = -999.;
    double decayTimeError_ = -999.;
    double decayTimeXY_ = -999.;
    double decayTimeXYError_ = -999.;
    double alpha_ = -999.;
    double alphaErr_ = -999.;
    double alphaXY_ = -999.;
    double alphaXYErr_ = -999.;

    const KinematicFitResult fit_;
    const reco::Vertex prodVertex_;
    int pvIndex_ = -1;
    std::string name_;
  };

  // WARNING: Displacements inherits from std::vector<T>, which does
  // not have a virtual destructor. Do not delete objects of this
  // class through a base class pointer to avoid potential memory
  // leaks or undefined behavior.

  class Displacements : public std::vector<Displacement>
  {
  public:
    const Displacement &get(const std::string &name) const
    {
      auto it = std::find_if(this->begin(), this->end(), [&](const Displacement &obj)
                             { return obj.name() == name; });
      if (it == this->end())
      {
        throw std::runtime_error("Cannot find displacement information with name '" + name + "'");
      }
      return *it;
    }
  };

}

#endif
