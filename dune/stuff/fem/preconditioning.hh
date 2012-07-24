#ifndef DUNE_STUFF_PRECONDITIONING_HH
#define DUNE_STUFF_PRECONDITIONING_HH

namespace Dune {

  namespace Stuff {

    namespace Fem {
// ! allow any class fullfilling the Operator concept to be used as a preconditioner
template< class Operator, template< class T, class F > class Solver, class RangeDiscreteFunctionType >
class OperatorBasedPreconditioner
{
  Operator& operator_;
  const typename RangeDiscreteFunctionType::FunctionSpaceType & range_space_;
  const bool right_preconditioning_;
  typedef Solver< RangeDiscreteFunctionType, Operator >
  SolverType;
  SolverType solver_;

public:
  OperatorBasedPreconditioner(Operator& op,
                              const typename RangeDiscreteFunctionType::FunctionSpaceType& range_space,
                              const bool right_preconditioning = false,
                              double solver_accuracy = 1e-7)
    : operator_(op)
      , range_space_(range_space)
      , right_preconditioning_(right_preconditioning)
      , solver_(operator_,
                solver_accuracy /*rel limit*/,
                solver_accuracy /*abs limit*/,
                1 /*not working iteration limit*/,
                false /*verbose*/)
  {}

  template< class VecType >
  void precondition(const VecType* tmp, VecType* dest) const {
    multOEM(tmp, dest);
  }

  template< class VECtype >
  void multOEM(const VECtype* x, VECtype* ret) const {
    operator_.multOEM(x, ret);
  }

  bool rightPrecondition() const {
    return right_preconditioning_;
  }
};

} // namespace Fem

} // namespace Stuff

} // namespace Dune

#endif // DUNE_STUFF_PRECONDITIONING_HH

/** Copyright (c) 2012, Rene Milk
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above copyright notice, this
   *    list of conditions and the following disclaimer.
   * 2. Redistributions in binary form must reproduce the above copyright notice,
   *    this list of conditions and the following disclaimer in the documentation
   *    and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   *
   * The views and conclusions contained in the software and documentation are those
   * of the authors and should not be interpreted as representing official policies,
   * either expressed or implied, of the FreeBSD Project.
   **/