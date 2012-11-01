// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticitylgdeformexplicitapp.

#include "ElasticityExplicitLgDeformData2DLinear.hh"

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_spaceDim = 2;

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_cellDim = 2;

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_numVertices = 3;

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_numBasis = 3;

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_matType = "ElasticPlaneStrain";

const char* pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_matDBFilename = "data/elasticplanestrain.spatialdb";

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_matLabel = "elastic strain 2-D";

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_dt =   1.00000000e-02;

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_dtStableExplicit =   1.50923086e-04;

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_gravityVec[] = {
  0.00000000e+00, -1.00000000e+08,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_vertices[] = {
  2.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  5.00000000e-01,
 -1.00000000e+00, -2.00000000e-01,
};

const int pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_cells[] = {
0,1,2,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_quadPts[] = {
 -3.33333333e-01, -3.33333333e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_quadWts[] = {
  2.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_basis[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_basisDerivRef[] = {
 -5.00000000e-01, -5.00000000e-01,
  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_fieldTIncr[] = {
  1.30000000e+00, -9.00000000e-01,
  1.40000000e+00,  1.50000000e+00,
  5.00000000e-01, -9.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_fieldT[] = {
  1.60000000e+00, -8.00000000e-01,
  9.00000000e-01,  7.00000000e-01,
 -2.00000000e-01, -1.10000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_fieldTmdt[] = {
  8.00000000e-01,  1.00000000e-01,
  5.00000000e-01,  3.00000000e-01,
 -1.00000000e-01, -6.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_valsResidual[] = {
 -5.23631714e+11,  5.83306915e+11,
  1.65371368e+11, -5.78689962e+11,
  3.58250720e+11, -4.62016160e+09,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::_valsJacobian[] = {
  4.58333333e+06,  4.58333333e+06,
  4.58333333e+06,  4.58333333e+06,
  4.58333333e+06,  4.58333333e+06,
};

pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::ElasticityExplicitLgDeformData2DLinear(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  dtStableExplicit = _dtStableExplicit;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityExplicitLgDeformData2DLinear::~ElasticityExplicitLgDeformData2DLinear(void)
{}


// End of file
