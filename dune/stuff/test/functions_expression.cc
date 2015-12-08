// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"
#include "functions.hh"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/expression.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif //HAVE_ALUGRID

struct ExpressionFunctionTest : public DS::FunctionTest< TESTFUNCTIONTYPE >
{
    virtual void check() const
    {
        Dune::Stuff::Common::Configuration config = FunctionType::default_config();
        const std::unique_ptr<const FunctionType> function(FunctionType::create(config));
        config["expression"] = "[2*x[0] 3*x[0] 4*x[0]; 1 sin(x[0]) 0; cos(x[0]) x[0] 0]";
        config["gradient"]   = "[2 0 0; 0 0 0; -sin(x[0]) 0 0]";
        config["gradient.1"] = "[3 0 0; cos(x[0]) 0 0; 1 0 0]";
        config["gradient.2"] = "[4 0 0; 0 0 0; 0 0 0]";
        const std::unique_ptr<const FunctionType> function2(FunctionType::create(config));
        const std::unique_ptr<const FunctionType> function3(new FunctionType(
                                                                           "x", "sin(x[0])", 3, FunctionType::static_id(), {"cos(x[0])", "0", "0"}));
    }
};

TEST_F(ExpressionFunctionTest, provides_required_methods) { this->check(); }

#else // HAVE_DUNE_GRID

TEST(DISABLED_ExpressionFunctionTest, provides_required_methods) {}

#endif // HAVE_DUNE_GRID
