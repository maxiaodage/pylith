// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/DirichletBIE.hh
 *
 * @brief C++ implementation of Dirichlet using Spectral Boundary integral
 * implementation defined in [Xiao et al, 2018, IJNAMG].
 */

#if !defined(pylith_bc_dirichletbie_hh)
#define pylith_bc_dirichletbie_hh

// Include directives ---------------------------------------------------
#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/ConstraintPointwise.hh" // ISA ConstraintPointwise

#include "pylith/topology/topologyfwd.hh" // USES Field

// Dirichlet ----------------------------------------------------
/// @brief Dirichlet (prescribed values at degrees of freedom) boundary
/// conditions with points on a boundary.
class pylith::bc::DirichletBIE :
    public pylith::bc::BoundaryCondition,
    public pylith::feassemble::ConstraintPointwise {

    friend class DirichletAuxiliaryFactory; // factory for auxiliary fields
    friend class TestDirichletBIE;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    DirichletBIE(void);

    /// Destructor.
    ~DirichletBIE(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Initialize boundary condition.
     *
     * @param[in] solution Solution field.
     */
    void initialize(const pylith::topology::Field& solution);

    /** Set constrained values in solution field.
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    void setSolution(pylith::topology::Field* solution,
                     const double t);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void _auxFieldSetup(const pylith::topology::Field& solution) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _boundaryMesh;   ///< Boundary mesh.
    pylith::topology::FieldBase::Description _description; ///< Description for constrained field.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    DirichletBIE(const DirichletBIE&); ///< Not implemented.
    const DirichletBIE& operator=(const DirichletBIE&); ///< Not implemented.

}; // class DirichletBIE

#endif // pylith_bc_dirichletbie_hh


// End of file
