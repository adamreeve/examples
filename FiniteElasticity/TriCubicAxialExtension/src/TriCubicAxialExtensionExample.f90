!> \file
!> $Id: TriCubicAxialExtensionExample.f90 20 2007-05-28 20:22:52Z cpb $
!> \author Kumar Mithraratne
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): 
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example FiniteElasticity/TriCubicAxialExtension/src/TriCubicAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TriCubicAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TriCubicAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM TRICUBICAXIALEXTENSIONEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber,MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex

  TYPE(CMISSBasisType) :: Basis1,Basis2
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType) :: Elements,LinearElements

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

  !parameters
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !local variables
  INTEGER(CMISSIntg) :: NumberOfSpatialCoordinates
  INTEGER(CMISSIntg) :: TotalNumberOfNodes
  INTEGER(CMISSIntg) :: NumberOfMeshDimensions,NumberOfMeshComponents
  INTEGER(CMISSIntg) :: TotalNumberOfElements
  INTEGER(CMISSIntg) :: FieldGeometryNumberOfVariables,FieldGeometryNumberOfComponents
  INTEGER(CMISSIntg) :: FieldFibreNumberOfVariables,FieldFibreNumberOfComponents
  INTEGER(CMISSIntg) :: FieldMaterialNumberOfVariables,FieldMaterialNumberOfComponents
  INTEGER(CMISSIntg) :: FieldDependentNumberOfVariables,FieldDependentNumberOfComponents

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Set diagnostics on
  CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=1
  NumberGlobalYElements=1
  NumberGlobalZElements=1
  NumberOfDomains=1

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  NumberOfSpatialCoordinates=3
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis functions (1 - tri-cubic Hermite ; 2 - tri-linear Lagrange)
  CALL CMISSBasisTypeInitialise(Basis1,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis1,Err)
  CALL CMISSBasisTypeSet(Basis1,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis1,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(Basis1,(/CMISSBasisCubicHermiteInterpolation, &
    & CMISSBasisCubicHermiteInterpolation,CMISSBasisCubicHermiteInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1, &
    & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisCreateFinish(Basis1,Err)

  CALL CMISSBasisTypeInitialise(Basis2,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,Basis2,Err)
  CALL CMISSBasisTypeSet(Basis2,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis2,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(Basis2,(/CMISSBasisLinearLagrangeInterpolation, &
    & CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis2, &
    & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisCreateFinish(Basis2,Err)

  !Create a mesh
  NumberOfMeshDimensions=3
  NumberOfMeshComponents=2
  TotalNumberOfElements=1
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL CMISSMeshNumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMeshNumberOfElementsSet(Mesh,TotalNumberOfElements,Err)

  !define nodes for the mesh
  TotalNumberOfNodes=8
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)

  !Cubic Hermite - Geometry
  CALL CMISSMeshElementsTypeInitialise(Elements,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,MeshComponentNumber,Basis1,Elements,Err)
  CALL CMISSMeshElementsNodesSet(Elements,1,(/1,2,3,4,5,6,7,8/),Err)
  CALL CMISSMeshElementsCreateFinish(Elements,Err)

  !Linear Lagrange - Hydrostatic pressure, Fibre angles and Material properties
  CALL CMISSMeshElementsTypeInitialise(LinearElements,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,LinearMeshComponentNumber,Basis2,LinearElements,Err)
  CALL CMISSMeshElementsNodesSet(LinearElements,1,(/1,2,3,4,5,6,7,8/),Err)
  CALL CMISSMeshElementsCreateFinish(LinearElements,Err)

  CALL CMISSMeshCreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  FieldGeometryNumberOfVariables=1
  FieldGeometryNumberOfComponents=3
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,MeshComponentNumber,Err)
  CALL CMISSFieldScalingTypeSetAndLock(GeometricField,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Set node parameters
  !node 1
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,1,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,1,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,1,3,1.0_CMISSDP,Err)

  !node 2
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,2,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,2,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,2,3,1.0_CMISSDP,Err)

  !node 3
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,3,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,3,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,3,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,3,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,3,3,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,3,3,1.0_CMISSDP,Err)

  !node 4
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,4,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,4,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,4,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,4,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,4,3,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,4,3,1.0_CMISSDP,Err)

  !node 5
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,5,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,5,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,5,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,5,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,5,3,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,5,3,1.0_CMISSDP,Err)

  !node 6
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,6,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,6,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,6,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,6,3,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,6,3,1.0_CMISSDP,Err)

  !node 7
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,7,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,7,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,7,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,7,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,7,3,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,7,3,1.0_CMISSDP,Err)

  !node 8
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,8,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,8,1,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,8,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,8,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,8,3,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,8,3,1.0_CMISSDP,Err)

  !Create a fibre field and attach it to the geometric field
  FieldFibreNumberOfVariables=1
  FieldFibreNumberOfComponents=3
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField,CMISSFieldUVariableType,FieldFibreNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,1,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,2,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,3,LinearMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Create a material field and attach it to the geometric field
  FieldMaterialNumberOfVariables=1
  FieldMaterialNumberOfComponents=2
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,1,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,2,LinearMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 3.0 respectively.
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,4.0_CMISSDP,Err)

  !Create a dependent field with two variables and four components
  FieldDependentNumberOfVariables=2
  FieldDependentNumberOfComponents=4
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSFieldCreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,FieldDependentNumberOfComponents, &
    & Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelUDelNVariableType,FieldDependentNumberOfComponents, &
    & Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,1,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,2,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,3,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,4,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,1,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,2,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,3,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,4,LinearMeshComponentNumber,Err)
  CALL CMISSFieldScalingTypeSetAndLock(DependentField,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldCreateFinish(DependentField,Err)

  !Create the EquationsSet
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,EquationsSet,Err)
  CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetMooneyRivlinSubtype,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-8.0_CMISSDP,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,1,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,1,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,1,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,1,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,1,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,1,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,1,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,1,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,1,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,1,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,1,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,1,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,1,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,1,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,1,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,1,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,1,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,1,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,1,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,1,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,1,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,1,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,1,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,1,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,2,1,CMISSBoundaryConditionFixed,1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,2,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,2,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,2,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,2,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,2,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,2,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,2,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,2,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,2,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,2,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,2,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,2,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,2,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,2,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,2,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,2,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,2,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,2,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,2,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,2,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,2,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,2,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,2,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,3,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,3,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,3,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,3,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,3,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,3,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,3,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,3,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,3,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,3,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,3,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,3,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,3,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,3,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,3,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,3,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,3,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,3,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,3,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,3,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,3,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,3,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,3,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,3,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,4,1,CMISSBoundaryConditionFixed,1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,4,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,4,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,4,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,4,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,4,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,4,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,4,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,4,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,4,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,4,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,4,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,4,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,4,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,4,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,4,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,4,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,4,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,4,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,4,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,4,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,4,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,4,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,4,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,5,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,5,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,5,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,5,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,5,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,5,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,5,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,5,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,5,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,5,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,5,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,5,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,5,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,5,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,5,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,5,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,5,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,5,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,5,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,5,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,5,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,5,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,5,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,5,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,6,1,CMISSBoundaryConditionFixed,1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,6,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,6,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,6,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,6,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,6,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,6,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,6,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,6,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,6,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,6,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,6,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,6,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,6,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,6,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,6,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,6,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,6,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,6,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,6,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,6,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,6,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,6,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,6,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,7,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,7,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,7,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,7,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,7,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,7,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,7,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,7,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,7,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,7,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,7,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,7,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,7,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,7,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,7,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,7,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,7,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,7,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,7,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,7,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,7,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,7,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,7,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,7,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,8,1,CMISSBoundaryConditionFixed,1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,8,1,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,8,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,8,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,8,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,8,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,8,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,8,1,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,8,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,8,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,8,2,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,8,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,8,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,8,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,8,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,8,2,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,8,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,2,8,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,3,8,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,4,8,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,5,8,3,CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,6,8,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,7,8,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,8,8,3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)

  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemFiniteElasticityType, &
    & CMISSProblemNoSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"TriCubicAxialExtension","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"TriCubicAxialExtension","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  !CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program completed."

  STOP

END PROGRAM TRICUBICAXIALEXTENSIONEXAMPLE

