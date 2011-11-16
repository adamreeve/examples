!> \file
!> \author Chris Bradley
!> \brief This is an example program which solves a weakly coupled finite-elasticity-Membrane-finite-elasticity-Membrane problem in two regions using openCMISS calls.
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
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):Jessica Jor, Thiranja Prasad Babarenda Gamage
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

!> \example InterfaceExamples/Coupled-FiniteElasticityMembrane-FiniteElasticityMembrane/src/Coupled-FiniteElasticityMembrane-FiniteElasticityMembrane.f90
!! Example program which solves a weakly coupled finite-elasticity-Membrane-finite-elasticity-Membrane problem in two regions using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/Coupled-FiniteElasticity-FiniteElasticity/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/Coupled-FiniteElasticity-FiniteElasticity/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM COUPLEDFINITEELASTICITYFINITEELASTICITY

  USE OPENCMISS
  USE MPI
  
#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
 
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: Region1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Region2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMappingBasisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMesh1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMesh2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeneratedMeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Mesh1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Mesh2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Decomposition1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Decomposition2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceDecompositionUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeometricFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FibreField1UserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FibreField2UserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: MaterialField1UserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2UserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentField1UserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentField2UserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField1UserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField2UserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: LagrangeFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: PenaltyFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: FibreField1NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FibreField1NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: FibreField2NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FibreField2NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: MaterialField1NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: DependentField1NumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: DependentField2NumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialField1NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: DependentField1NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: DependentField2NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: InterfaceUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: InterfaceConditionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CoupledProblemUserNumber=1

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI,NUMBER_OF_NODE_XI
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  INTEGER(CMISSIntg) :: EquationsSet1Index,EquationsSet2Index
  INTEGER(CMISSIntg) :: InterfaceConditionIndex
  INTEGER(CMISSIntg) :: Mesh1Index,Mesh2Index
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: y_element_idx,z_element_idx,mesh_local_y_node,mesh_local_z_node
  REAL(CMISSDP) :: XI2(2),XI3(3), HEIGHT,WIDTH,LENGTH
  LOGICAL :: UNCOUPLED,FORCE_BC,INCOMPRESSIBLE,PENALTY

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis1,Basis2,InterfaceBasis,InterfaceMappingBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem1,CoordinateSystem2,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition1,Decomposition2,InterfaceDecomposition
  TYPE(CMISSEquationsType) :: Equations1,Equations2
  TYPE(CMISSEquationsSetType) :: EquationsSet1,EquationsSet2
  TYPE(CMISSFieldType) :: GeometricField1,GeometricField2,InterfaceGeometricField,DependentField1,DependentField2,LagrangeField, &
    & EquationsSetField1,EquationsSetField2,PenaltyField
  TYPE(CMISSFieldType) :: FibreField1,FibreField2,MaterialField1,MaterialField2
  TYPE(CMISSFieldsType) :: Fields1,Fields2,Fields3,InterfaceFields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh1,GeneratedMesh2,InterfaceGeneratedMesh
  TYPE(CMISSMeshType) :: Mesh1,Mesh2,InterfaceMesh
  TYPE(CMISSInterfaceMeshConnectivityType) :: InterfaceMeshConnectivity
  TYPE(CMISSInterfaceType) :: Interface
  TYPE(CMISSInterfaceConditionType) :: InterfaceCondition
  TYPE(CMISSInterfaceEquationsType) :: InterfaceEquations
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: CoupledProblem
  TYPE(CMISSRegionType) :: Region1,Region2,WorldRegion
  TYPE(CMISSSolverType) :: NonLinearSolver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: NonLinearSolverEquations
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: Err
  
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

  PENALTY=.TRUE.
  FORCE_BC=.TRUE. 
  UNCOUPLED=.FALSE. 
  IF(UNCOUPLED.EQV..TRUE.) THEN
    NUMBER_GLOBAL_X_ELEMENTS=2
    WIDTH=2.0_CMISSDP
  ELSE
    NUMBER_GLOBAL_X_ELEMENTS=1
    WIDTH=1.0_CMISSDP
  ENDIF
  HEIGHT=1.0_CMISSDP
  LENGTH=1.0_CMISSDP
  NUMBER_GLOBAL_Y_ELEMENTS=1
  NUMBER_GLOBAL_Z_ELEMENTS=0
  INTERPOLATION_TYPE = CMISSBasisLinearLagrangeInterpolation

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Set error handling mode
  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)
 
  !Set diganostics for testing
  !CALL CMISSDiagnosticsSetOn(CMISSAllDiagType,[1,2,3,4,5],"Diagnostics",["SOLVER_MAPPING_CALCULATE         ", &
  !  & "SOLVER_MATRIX_STRUCTURE_CALCULATE"],Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  
  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(1) << == '
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem1,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystem1UserNumber,CoordinateSystem1,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem1,3,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem1,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem1,Err)

  !Start the creation of a new RC coordinate system for the second region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(2) << == '
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem2,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystem2UserNumber,CoordinateSystem2,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem2,3,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem2,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem2,Err) 
  
  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION(1) << == '
  CALL CMISSRegionTypeInitialise(Region1,Err)
  CALL CMISSRegionCreateStart(Region1UserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegionLabelSet(Region1,"Region1",Err)
  !Set the regions coordinate system to the 1st RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region1,CoordinateSystem1,Err)
  !Finish the creation of the first region
  CALL CMISSRegionCreateFinish(Region1,Err)

  !Start the creation of the second region
  PRINT *, ' == >> CREATING REGION(2) << == '
  CALL CMISSRegionTypeInitialise(Region2,Err)
  CALL CMISSRegionCreateStart(Region2UserNumber,WorldRegion,Region2,Err)
  CALL CMISSRegionLabelSet(Region2,"Region2",Err)
  !Set the regions coordinate system to the 2nd RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region2,CoordinateSystem2,Err)
  !Finish the creation of the second region
  CALL CMISSRegionCreateFinish(Region2,Err) 

  SELECT CASE(INTERPOLATION_TYPE)
  CASE(CMISSBasisLinearLagrangeInterpolation)
    NUMBER_OF_GAUSS_XI=2
  CASE(CMISSBasisQuadraticLagrangeInterpolation)
    NUMBER_OF_GAUSS_XI=3
  CASE(CMISSBasisCubicLagrangeInterpolation,CMISSBasisCubicHermiteInterpolation)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT

  !Start the creation of a bI/tri-linear-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(1) << == '
  CALL CMISSBasisTypeInitialise(Basis1,Err)
  CALL CMISSBasisCreateStart(Basis1UserNumber,Basis1,Err)
  !Set the basis to be a tri-XXX basis
  CALL CMISSBasisNumberOfXiSet(Basis1,2,Err)
  CALL CMISSBasisInterpolationXiSet(Basis1,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis1,Err) 
   
  !Start the creation of a bI/tri-linear-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(2) << == '
  CALL CMISSBasisTypeInitialise(Basis2,Err)
  CALL CMISSBasisCreateStart(Basis2UserNumber,Basis2,Err)
  !Set the basis to be a tri-XXX basis
  CALL CMISSBasisNumberOfXiSet(Basis2,2,Err)
  CALL CMISSBasisInterpolationXiSet(Basis2,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis2,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis2,Err) 
  
  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING GENERATED MESH(1) << == '
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh1,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshCreateStart(GeneratedMesh1UserNumber,Region1,GeneratedMesh1,Err)
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh1,CMISSGeneratedMeshRegularMeshType,Err)
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh1,Basis1,Err)
  !Define the mesh on the first region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh1,[WIDTH,HEIGHT,0.0_CMISSDP],Err) 
  CALL CMISSGeneratedMeshOriginSet(GeneratedMesh1,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err) 
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh1,[NUMBER_GLOBAL_X_ELEMENTS, &
    & NUMBER_GLOBAL_Y_ELEMENTS],Err)
  !Finish the creation of a generated mesh in the first region 
  CALL CMISSMeshTypeInitialise(Mesh1,Err) 
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh1,Mesh1UserNumber,Mesh1,Err)

  !Start the creation of a generated mesh in the second region
  PRINT *, ' == >> CREATING GENERATED MESH(2) << == '
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh2,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshCreateStart(GeneratedMesh2UserNumber,Region2,GeneratedMesh2,Err)
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh2,CMISSGeneratedMeshRegularMeshType,Err)
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh2,Basis2,Err)
  !Define the mesh on the first region
  IF(UNCOUPLED.EQV..TRUE.) THEN
    CALL CMISSGeneratedMeshOriginSet(GeneratedMesh2,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err) 
  ELSE
    CALL CMISSGeneratedMeshOriginSet(GeneratedMesh2,[WIDTH,0.0_CMISSDP,0.0_CMISSDP],Err) 
  ENDIF
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh2,[WIDTH,HEIGHT,0.0_CMISSDP],Err) 
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh2,[NUMBER_GLOBAL_X_ELEMENTS, &
    & NUMBER_GLOBAL_Y_ELEMENTS],Err)
  !Finish the creation of a generated mesh in the first region 
  CALL CMISSMeshTypeInitialise(Mesh2,Err) 
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh2,Mesh2UserNumber,Mesh2,Err)

  !Create a decomposition for mesh1
  PRINT *, ' == >> CREATING MESH(1) DECOMPOSITION << == '
  CALL CMISSDecompositionTypeInitialise(Decomposition1,Err)
  CALL CMISSDecompositionCreateStart(Decomposition1UserNumber,Mesh1,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition1,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition1,Err)

  !Create a decomposition for mesh2
  PRINT *, ' == >> CREATING MESH(2) DECOMPOSITION << == '
  CALL CMISSDecompositionTypeInitialise(Decomposition2,Err)
  CALL CMISSDecompositionCreateStart(Decomposition2UserNumber,Mesh2,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition2,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition2,Err)
  
  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH(1) GEOMETRIC FIELD << == '
  CALL CMISSFieldTypeInitialise(GeometricField1,Err)
  CALL CMISSFieldCreateStart(GeometricField1UserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField1,Decomposition1,Err)
  CALL CMISSFieldTypeSet(GeometricField1,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricField1,GeometricField1NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField1,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField1,CMISSFieldUVariableType,GeometricField1NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(GeometricField1,Err)

  !Start to create a default (geometric) field on the second region
  PRINT *, ' == >> CREATING MESH(2) GEOMETRIC FIELD << == '
  CALL CMISSFieldTypeInitialise(GeometricField2,Err)
  CALL CMISSFieldCreateStart(GeometricField2UserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField2,Decomposition2,Err)
  CALL CMISSFieldTypeSet(GeometricField2,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricField2,GeometricField2NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField2,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField2,CMISSFieldUVariableType,GeometricField2NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField2,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField2,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField2,CMISSFieldUVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(GeometricField2,Err)

  !Update the geometric field parameters for the first field
  PRINT *, ' == >> INITIALIZE GEOMETRIC FIELD(1) FROM GENERATED MESH(1) << == '
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField1,GeneratedMesh1,Err)

  !Update the geometric field parameters for the second field
  PRINT *, ' == >> INITIALIZE GEOMETRIC FIELD(2) FROM GENERATED MESH(2) << == '
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField2,GeneratedMesh2,Err)

  !######################################### FINITE ELASTICITY ##########################################

  !Start to create a default (fibre) field on the first region
  PRINT *, ' == >> CREATING MESH(1) FIBRE FIELD << == '
  CALL CMISSFieldTypeInitialise(FibreField1,Err)
  CALL CMISSFieldCreateStart(FibreField1UserNumber,Region1,FibreField1,Err)
  CALL CMISSFieldTypeSet(FibreField1,CMISSFieldFibreType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(FibreField1,Decomposition1,Err)
  CALL CMISSFieldGeometricFieldSet(FibreField1,GeometricField1,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField1,FibreField1NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(FibreField1,CMISSFieldUVariableType,"Fibre",Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField1,CMISSFieldUVariableType,FibreField1NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(FibreField1,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField1,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField1,CMISSFieldUVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(FibreField1,Err)

  !Start to create a default (fibre) field on the second region
  PRINT *, ' == >> CREATING MESH(2) FIBRE FIELD << == '
  CALL CMISSFieldTypeInitialise(FibreField2,Err)
  CALL CMISSFieldCreateStart(FibreField2UserNumber,Region2,FibreField2,Err)
  CALL CMISSFieldTypeSet(FibreField2,CMISSFieldFibreType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(FibreField2,Decomposition2,Err)
  CALL CMISSFieldGeometricFieldSet(FibreField2,GeometricField2,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField2,FibreField2NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(FibreField2,CMISSFieldUVariableType,"Fibre",Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField2,CMISSFieldUVariableType,FibreField2NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(FibreField2,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField2,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField2,CMISSFieldUVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(FibreField2,Err)

  !Start to create a material field on the first region
  PRINT *, ' == >> CREATING MESH(1) MATERIAL FIELD << == '
  CALL CMISSFieldTypeInitialise(MaterialField1,Err)
  CALL CMISSFieldCreateStart(MaterialField1UserNumber,Region1,MaterialField1,Err)
  CALL CMISSFieldTypeSet(MaterialField1,CMISSFieldMaterialType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(MaterialField1,Decomposition1,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialField1,GeometricField1,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField1,MaterialField1NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(MaterialField1,CMISSFieldUVariableType,"Material",Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField1,CMISSFieldUVariableType,MaterialField1NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(MaterialField1,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialField1,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialField1,CMISSFieldUVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(MaterialField1,Err)

  !Start to create a material field on the second region
  PRINT *, ' == >> CREATING MESH(2) MATERIAL FIELD << == '
  CALL CMISSFieldTypeInitialise(MaterialField2,Err)
  CALL CMISSFieldCreateStart(MaterialField2UserNumber,Region2,MaterialField2,Err)
  CALL CMISSFieldTypeSet(MaterialField2,CMISSFieldMaterialType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(MaterialField2,Decomposition2,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialField2,GeometricField2,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField2,MaterialField2NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(MaterialField2,CMISSFieldUVariableType,"Material",Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField2,CMISSFieldUVariableType,MaterialField2NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(MaterialField2,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialField2,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialField2,CMISSFieldUVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(MaterialField2,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 1.0 respectively on the first region
  PRINT *, ' == >> INITIALIZE MATERIAL FIELD(1) << == '
  CALL CMISSFieldComponentValuesInitialise(MaterialField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)
  !Set membrane thickness (only specifiy this for 3D space)
  CALL CMISSFieldComponentValuesInitialise(MaterialField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,0.1_CMISSDP,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 1.0 respectively on the second region
  PRINT *, ' == >> INITIALIZE MATERIAL FIELD(2) << == '
  CALL CMISSFieldComponentValuesInitialise(MaterialField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)
  !Set membrane thickness (only specifiy this for 3D space)
  CALL CMISSFieldComponentValuesInitialise(MaterialField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,0.1_CMISSDP,Err)

  !Start to create a dependent field on the first region
  PRINT *, ' == >> CREATING MESH(1) DEPENDENT FIELD << == '
  CALL CMISSFieldTypeInitialise(DependentField1,Err)
  CALL CMISSFieldCreateStart(DependentField1UserNumber,Region1,DependentField1,Err)
  CALL CMISSFieldTypeSet(DependentField1,CMISSFieldGeneralType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(DependentField1,Decomposition1,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField1,GeometricField1,Err)
  CALL CMISSFieldDependentTypeSet(DependentField1,CMISSFieldDependentType,Err) 
  CALL CMISSFieldNumberOfVariablesSet(DependentField1,DependentField1NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(DependentField1,CMISSFieldUVariableType,"Dependent1",Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField1,CMISSFieldUVariableType,DependentField1NumberOfComponents,Err)  
  CALL CMISSFieldNumberOfComponentsSet(DependentField1,CMISSFieldDelUDelNVariableType,DependentField1NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldUVariableType,3,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldDelUDelNVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldDelUDelNVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldDelUDelNVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(DependentField1,Err)

  !Start to create a dependent field on the second region
  PRINT *, ' == >> CREATING MESH(2) DEPENDENT FIELD << == '
  CALL CMISSFieldTypeInitialise(DependentField2,Err)
  CALL CMISSFieldCreateStart(DependentField2UserNumber,Region2,DependentField2,Err)
  CALL CMISSFieldTypeSet(DependentField2,CMISSFieldGeneralType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(DependentField2,Decomposition2,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField2,GeometricField2,Err)
  CALL CMISSFieldDependentTypeSet(DependentField2,CMISSFieldDependentType,Err) 
  CALL CMISSFieldNumberOfVariablesSet(DependentField2,DependentField2NumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(DependentField2,CMISSFieldUVariableType,"Dependent2",Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField2,CMISSFieldUVariableType,DependentField2NumberOfComponents,Err)  
  CALL CMISSFieldNumberOfComponentsSet(DependentField2,CMISSFieldDelUDelNVariableType,DependentField2NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldUVariableType,3,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldDelUDelNVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldDelUDelNVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldDelUDelNVariableType,3,1,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(DependentField2,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  PRINT *, ' == >> INITIALIZE DEPENDENT FIELD(1) << == '
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  PRINT *, ' == >> INITIALIZE DEPENDENT FIELD(2) << == '
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)

  !Create the equations set for the first region
  PRINT *, ' == >> CREATING EQUATION SET(1) << == '
  CALL CMISSFieldTypeInitialise(EquationsSetField1,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSet1,Err)
  !Set the equations set to be a Finite Elasticity problem
  CALL CMISSEquationsSetCreateStart(EquationsSet1UserNumber,Region1,FibreField1,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetMembraneSubtype,EquationsSetField1UserNumber,&
    & EquationsSetField1,EquationsSet1,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet1,Err)

  !Create the equations set for the second region
  PRINT *, ' == >> CREATING EQUATION SET(2) << == '
  CALL CMISSFieldTypeInitialise(EquationsSetField2,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSet2,Err)
  !Set the equations set to be a Finite Elasticity problem
  CALL CMISSEquationsSetCreateStart(EquationsSet2UserNumber,Region2,FibreField2,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetMembraneSubtype,EquationsSetField2UserNumber,&
    & EquationsSetField2,EquationsSet2,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet2,Err)

  !Create the equations set dependent field variables for the first equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD EQUATION SET(1) << == '
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet1,DependentField1UserNumber,DependentField1,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet1,Err)

  !Create the equations set dependent field variables for the second equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD EQUATION SET(2) << == '
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet2,DependentField2UserNumber,DependentField2,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet2,Err)

  !Create the equations set dependent field variables for the first equations set
  PRINT *, ' == >> CREATING MATERIAL FIELD EQUATION SET(1) << == '
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet1,MaterialField1UserNumber,MaterialField1,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet1,Err)

  !Create the equations set dependent field variables for the second equations set
  PRINT *, ' == >> CREATING MATERIAL FIELD EQUATION SET(2) << == '
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet2,MaterialField2UserNumber,MaterialField2,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet2,Err)

  !Create the equations set equations for the first equations set
  PRINT *, ' == >> CREATING EQUATIONS(1) << == '
  CALL CMISSEquationsTypeInitialise(Equations1,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet1,Equations1,Err)
  !Set the equations matrices sparsity type
  !CALL CMISSEquationsSparsityTypeSet(Equations1,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations1,CMISSEquationsFullMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations1,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations1,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations1,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations1,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet1,Err)

  !Create the equations set equations for the second equations set
  PRINT *, ' == >> CREATING EQUATIONS(2) << == '
  CALL CMISSEquationsTypeInitialise(Equations2,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet2,Equations2,Err)
  !Set the equations matrices sparsity type
  !CALL CMISSEquationsSparsityTypeSet(Equations2,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations2,CMISSEquationsFullMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations2,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations2,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations2,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations2,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet2,Err)

!  !######################################### START CREATING INTERFACE ##########################################

  !Create an interface between the two meshes
  PRINT *, ' == >> CREATING INTERFACE << == '
  CALL CMISSInterfaceTypeInitialise(Interface,Err)
  CALL CMISSInterfaceCreateStart(InterfaceUserNumber,WorldRegion,Interface,Err)
  CALL CMISSInterfaceLabelSet(Interface,"Interface",Err)
  !Add in the two meshes
  CALL CMISSInterfaceMeshAdd(Interface,Mesh1,Mesh1Index,Err)
  CALL CMISSInterfaceMeshAdd(Interface,Mesh2,Mesh2Index,Err)
  !Finish creating the interface
  CALL CMISSInterfaceCreateFinish(Interface,Err)

  !Start the creation of a (bi)-linear-Lagrange basis
  PRINT *, ' == >> CREATING INTERFACE BASIS << == '
  CALL CMISSBasisTypeInitialise(InterfaceBasis,Err)
  CALL CMISSBasisCreateStart(InterfaceBasisUserNumber,InterfaceBasis,Err)
  !Set the basis to be a linear Lagrange basis
  CALL CMISSBasisNumberOfXiSet(InterfaceBasis,1,Err)
  CALL CMISSBasisInterpolationXiSet(InterfaceBasis,[CMISSBasisLinearLagrangeInterpolation],Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(InterfaceBasis,Err)

  !Start the creation of a bi-linear-Lagrange basis
  PRINT *, ' == >> CREATING INTERFACE MAPPING BASIS << == '
  CALL CMISSBasisTypeInitialise(InterfaceMappingBasis,Err)
  CALL CMISSBasisCreateStart(InterfaceMappingBasisUserNumber,InterfaceMappingBasis,Err)
  !Set the basis to be a linear Lagrange basis
  CALL CMISSBasisNumberOfXiSet(InterfaceMappingBasis,1,Err)
  CALL CMISSBasisInterpolationXiSet(InterfaceMappingBasis,[CMISSBasisLinearLagrangeInterpolation],Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(InterfaceMappingBasis,Err)
  
  !Start the creation of a generated mesh for the interface
  PRINT *, ' == >> CREATING INTERFACE GENERATED MESH << == '
  CALL CMISSGeneratedMeshTypeInitialise(InterfaceGeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(InterfaceGeneratedMeshUserNumber,Interface,InterfaceGeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(InterfaceGeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(InterfaceGeneratedMesh,InterfaceBasis,Err)
  CALL CMISSGeneratedMeshOriginSet(InterfaceGeneratedMesh,[WIDTH,0.0_CMISSDP,0.0_CMISSDP],Err)
  CALL CMISSGeneratedMeshExtentSet(InterfaceGeneratedMesh,[0.0_CMISSDP,HEIGHT,0.0_CMISSDP],Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(InterfaceGeneratedMesh,[NUMBER_GLOBAL_Y_ELEMENTS],Err)
  !Finish the creation of a generated mesh in interface
  CALL CMISSMeshTypeInitialise(InterfaceMesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(InterfaceGeneratedMesh,InterfaceMeshUserNumber,InterfaceMesh,Err)

  !Couple the interface meshes
  PRINT *, ' == >> CREATING INTERFACE MESHES CONNECTIVITY << == '
  CALL CMISSInterfaceMeshConnectivityTypeInitialise(InterfaceMeshConnectivity,Err)
  CALL CMISSInterfaceMeshConnectivityCreateStart(Interface,InterfaceMesh,InterfaceMeshConnectivity,Err)
  CALL CMISSInterfaceMeshConnectivitySetBasis(InterfaceMeshConnectivity,InterfaceMappingBasis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,4)
    NUMBER_OF_NODE_XI=2
  CASE(2)
    NUMBER_OF_NODE_XI=3
  CASE(3)
    NUMBER_OF_NODE_XI=4
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
    DO y_element_idx=1,NUMBER_GLOBAL_Y_ELEMENTS
      !Map the interface element to the elements in mesh 1
      CALL CMISSInterfaceMeshConnectivityElementNumberSet(InterfaceMeshConnectivity,y_element_idx,Mesh1Index, &
        y_element_idx*NUMBER_GLOBAL_X_ELEMENTS,Err)
      XI2 = [ 1.0_CMISSDP, 0.0_CMISSDP ]
      CALL CMISSInterfaceMeshConnectivityElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh1Index, &
        & y_element_idx*NUMBER_GLOBAL_X_ELEMENTS,1,1,XI2,Err)
      XI2 = [ 1.0_CMISSDP, 1.0_CMISSDP ]
      CALL CMISSInterfaceMeshConnectivityElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh1Index, &
        & y_element_idx*NUMBER_GLOBAL_X_ELEMENTS,2,1,XI2,Err)
      !Map the interface element to the elements in mesh 2
      CALL CMISSInterfaceMeshConnectivityElementNumberSet(InterfaceMeshConnectivity,y_element_idx,Mesh2Index, &
        & 1+(y_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS,Err)
      DO mesh_local_y_node = 1,NUMBER_OF_NODE_XI-1
        XI2 = [ 0.0_CMISSDP, REAL(mesh_local_y_node-1,CMISSDP)/REAL(NUMBER_OF_NODE_XI-1,CMISSDP) ]
        CALL CMISSInterfaceMeshConnectivityElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh2Index, &
          & 1+(y_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS,1,1,XI2,Err)
        XI2 = [ 0.0_CMISSDP, REAL(mesh_local_y_node,CMISSDP)/REAL(NUMBER_OF_NODE_XI-1,CMISSDP) ]
        CALL CMISSInterfaceMeshConnectivityElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh2Index, &
          & 1+(y_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS,2,1,XI2,Err)
      ENDDO !mesh_local_y_node
    ENDDO !y_element_idx
  CALL CMISSInterfaceMeshConnectivityCreateFinish(InterfaceMeshConnectivity,Err)

  !Create a decomposition for the interface mesh
  PRINT *, ' == >> CREATING INTERFACE DECOMPOSITION << == '
  CALL CMISSDecompositionTypeInitialise(InterfaceDecomposition,Err)
  CALL CMISSDecompositionCreateStart(InterfaceDecompositionUserNumber,InterfaceMesh,InterfaceDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(InterfaceDecomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(InterfaceDecomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(InterfaceDecomposition,Err)

  !Start to create a default (geometric) field on the Interface
  PRINT *, ' == >> CREATING INTERFACE GEOMETRIC FIELD << == '
  CALL CMISSFieldTypeInitialise(InterfaceGeometricField,Err)
  CALL CMISSFieldCreateStart(InterfaceGeometricFieldUserNumber,Interface,InterfaceGeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(InterfaceGeometricField,InterfaceDecomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(InterfaceGeometricField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(InterfaceGeometricField,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(InterfaceGeometricField,CMISSFieldUVariableType,3,1,Err)
   !Finish creating the first field
  CALL CMISSFieldCreateFinish(InterfaceGeometricField,Err)

  !Update the geometric field parameters for the interface field
  PRINT *, ' == >> INITIALIZE INTERFACE GEOMETRIC FIELD FROM INTERFACE GENERATED MESH << == '
  CALL CMISSGeneratedMeshGeometricParametersCalculate(InterfaceGeometricField,InterfaceGeneratedMesh,Err)

  !Create an interface condition between the two meshes
  PRINT *, ' == >> CREATING INTERFACE CONDITIONS << == '
  CALL CMISSInterfaceConditionTypeInitialise(InterfaceCondition,Err)
  CALL CMISSInterfaceConditionCreateStart(InterfaceConditionUserNumber,Interface,InterfaceGeometricField, &
    & InterfaceCondition,Err)
  !Specify the method for the interface condition
  IF(PENALTY) THEN
    CALL CMISSInterfaceConditionMethodSet(InterfaceCondition,CMISSInterfaceConditionPenaltyMethod,Err)
  ELSE
    CALL CMISSInterfaceConditionMethodSet(InterfaceCondition,CMISSInterfaceConditionLagrangeMultipliers,Err)
  ENDIF
  !Specify the type of interface condition operator
  CALL CMISSInterfaceConditionOperatorSet(InterfaceCondition,CMISSInterfaceConditionFieldContinuityOperator,Err)
  !Add in the dependent variables from the equations sets
  CALL CMISSInterfaceConditionDependentVariableAdd(InterfaceCondition,Mesh1Index,EquationsSet1, &
    & CMISSFieldUVariableType,Err)
  CALL CMISSInterfaceConditionDependentVariableAdd(InterfaceCondition,Mesh2Index,EquationsSet2, &
    & CMISSFieldUVariableType,Err)
  !Finish creating the interface condition
  CALL CMISSInterfaceConditionCreateFinish(InterfaceCondition,Err)

  !Create the Lagrange multipliers field
  PRINT *, ' == >> CREATING INTERFACE LAGRANGE FIELD << == '
  CALL CMISSFieldTypeInitialise(LagrangeField,Err)
  CALL CMISSInterfaceConditionLagrangeFieldCreateStart(InterfaceCondition,LagrangeFieldUserNumber,LagrangeField,Err)
  !Finish the Lagrange multipliers field
  CALL CMISSInterfaceConditionLagrangeFieldCreateFinish(InterfaceCondition,Err)

  IF(PENALTY) THEN
    !Create the Penalty field
    PRINT *, ' == >> CREATING INTERFACE PENALTY FIELD << == '
    CALL CMISSFieldTypeInitialise(PenaltyField,Err)
    CALL CMISSInterfaceConditionPenaltyFieldCreateStart(InterfaceCondition,PenaltyFieldUserNumber,PenaltyField,Err)
    !Finish the Penalty field
    CALL CMISSInterfaceConditionPenaltyFieldCreateFinish(InterfaceCondition,Err)
    !Set the penalty field coefficients
    CALL CMISSFieldComponentValuesInitialise(PenaltyField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,100.0_CMISSDP,Err)
    CALL CMISSFieldComponentValuesInitialise(PenaltyField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,100.0_CMISSDP,Err)
    CALL CMISSFieldComponentValuesInitialise(PenaltyField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,100.0_CMISSDP,Err)
  ENDIF

  !Create the interface condition equations
  PRINT *, ' == >> CREATING INTERFACE EQUATIONS << == '
  CALL CMISSInterfaceEquationsTypeInitialise(InterfaceEquations,Err)
  CALL CMISSInterfaceConditionEquationsCreateStart(InterfaceCondition,InterfaceEquations,Err)
  !Set the interface equations sparsity
  !CALL CMISSInterfaceEquationsSparsitySet(InterfaceEquations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSInterfaceEquationsSparsitySet(InterfaceEquations,CMISSEquationsFullMatrices,Err)
  !Set the interface equations output
  CALL CMISSInterfaceEquationsOutputTypeSet(InterfaceEquations,CMISSEquationsMatrixOutput,Err)
  !Finish creating the interface equations
  CALL CMISSInterfaceConditionEquationsCreateFinish(InterfaceCondition,Err)

  !######################################### FINISH CREATING INTERFACE ##########################################

  !Start the creation of a coupled problem.
  PRINT *, ' == >> CREATING PROBLEM << == '
  CALL CMISSProblemTypeInitialise(CoupledProblem,Err)
  CALL CMISSProblemCreateStart(CoupledProblemUserNumber,CoupledProblem,Err)
  CALL CMISSProblemSpecificationSet(CoupledProblem,CMISSProblemElasticityClass, &
    & CMISSProblemFiniteElasticityType,CMISSProblemNoSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(CoupledProblem,Err)

  !Start the creation of the problem control loop for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM CONTROL LOOP << == '
  CALL CMISSProblemControlLoopCreateStart(CoupledProblem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(CoupledProblem,Err)
 
  !Start the creation of the problem solver for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVERS << == '
  CALL CMISSSolverTypeInitialise(NonLinearSolver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(CoupledProblem,Err)
  CALL CMISSProblemSolverGet(CoupledProblem,CMISSControlLoopNode,1,NonLinearSolver,Err)
  !CALL CMISSSolverOutputTypeSet(NonLinearSolver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(NonLinearSolver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(NonLinearSolver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(NonLinearSolver,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(NonLinearSolver,CMISSSolverSolverMatrixOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonLinearSolver,CMISSSolverNewtonJacobianFDCalculated,Err)
  !CALL CMISSSolverNewtonMaximumFunctionEvaluationsSet(NonLinearSolver,1,Err)
  CALL CMISSSolverNewtonMaximumIterationsSet(NonLinearSolver,1000,Err)
  !CALL CMISSSolverLibraryTypeSet(NonLinearSolver,CMISSSolverMUMPSLibrary,Err)
  CALL CMISSSolverNewtonLinearSolverGet(NonLinearSolver,LinearSolver,Err)
  CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverLAPACKLibrary,Err)
  !CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearIterativeSolveType,Err)
  !CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolver,1,Err)
  !CALL CMISSSolverOutputTypeSet(LinearSolver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(LinearSolver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(LinearSolver,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(LinearSolver,CMISSSolverSolverMatrixOutput,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(CoupledProblem,Err)

  !Start the creation of the problem solver equations for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVER EQUATIONS << == '
  CALL CMISSSolverTypeInitialise(NonLinearSolver,Err)
  CALL CMISSSolverEquationsTypeInitialise(NonLinearSolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(CoupledProblem,Err)
  !Get the solve equations
  CALL CMISSProblemSolverGet(CoupledProblem,CMISSControlLoopNode,1,NonLinearSolver,Err)
  CALL CMISSSolverSolverEquationsGet(NonLinearSolver,NonLinearSolverEquations,Err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquationsSparsityTypeSet(NonLinearSolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(NonLinearSolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the first equations set
  CALL CMISSSolverEquationsEquationsSetAdd(NonLinearSolverEquations,EquationsSet1,EquationsSet1Index,Err)
  !Add in the second equations set
  CALL CMISSSolverEquationsEquationsSetAdd(NonLinearSolverEquations,EquationsSet2,EquationsSet2Index,Err)
  IF(UNCOUPLED.EQV..FALSE.) THEN
    !Add in the interface condition
    CALL CMISSSolverEquationsInterfaceConditionAdd(NonLinearSolverEquations,InterfaceCondition,InterfaceConditionIndex,Err)
  ENDIF
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(CoupledProblem,Err)

  !Export the fields
  PRINT *, ' == >> EXPORTING FIELDS << == '
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSFieldsTypeInitialise(Fields1,Err)
    CALL CMISSFieldsTypeCreate(Region1,Fields1,Err)
    CALL CMISSFieldIONodesExport(Fields1,"3DCoupled-FiniteElasticity-FiniteElasticity_1_Initial","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields1,"3DCoupled-FiniteElasticity-FiniteElasticity_1_Initial","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields1,Err)
    CALL CMISSFieldsTypeInitialise(Fields2,Err)
    CALL CMISSFieldsTypeCreate(Region2,Fields2,Err)
    CALL CMISSFieldIONodesExport(Fields2,"3DCoupled-FiniteElasticity-FiniteElasticity_2_Initial","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields2,"3DCoupled-FiniteElasticity-FiniteElasticity_2_Initial","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields2,Err)
    IF(UNCOUPLED.EQV..FALSE.) THEN
      CALL CMISSFieldsTypeInitialise(Fields3,Err)
      CALL CMISSFieldsTypeCreate(Interface,Fields3,Err)
      CALL CMISSFieldIONodesExport(Fields3,"3DCoupled-FiniteElasticity-FiniteElasticity_Interface_Initial","FORTRAN",Err)
      CALL CMISSFieldIOElementsExport(Fields3,"3DCoupled-FiniteElasticity-FiniteElasticity_Interface_Initial","FORTRAN",Err)
      CALL CMISSFieldsTypeFinalise(Fields3,Err)
    ENDIF
  ENDIF

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(NonLinearSolverEquations,BoundaryConditions,Err)

  IF(UNCOUPLED.EQV..FALSE.) THEN
    !============================================================================================================================
    ! COUPLED
    !============================================================================================================================
    !Start the creation of the solver equations boundary conditions for the first region
    PRINT *, ' == >> CREATING BOUNDARY CONDITIONS(1) << == '
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,1,1, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,3,1, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,1,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,2,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,1,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,2,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,3,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,4,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    !Start the creation of the solver equations boundary conditions for the second region
    PRINT *, ' == >> CREATING BOUNDARY CONDITIONS(2) << == '
    IF(FORCE_BC.EQV..TRUE.) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldDelUDelNVariableType,1,1,2,1, &
        & CMISSBoundaryConditionFixed, &
        & -0.2_CMISSDP,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldDelUDelNVariableType,1,1,4,1, &
        & CMISSBoundaryConditionFixed, &
        & -0.2_CMISSDP,Err)
    ELSE
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,2,1, &
        & CMISSBoundaryConditionFixed, &
        & 2.1_CMISSDP,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,4,1, &
        & CMISSBoundaryConditionFixed, &
        & 2.1_CMISSDP,Err)
    ENDIF

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,1,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,2,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,1,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,2,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,LagrangeField,CMISSFieldUVariableType,1,1,1,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,LagrangeField,CMISSFieldUVariableType,1,1,2,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,LagrangeField,CMISSFieldUVariableType,1,1,1,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
  ELSE
    !============================================================================================================================
    ! UNCOUPLED
    !============================================================================================================================
    !Start the creation of the solver equations boundary conditions for the first region
    PRINT *, ' == >> CREATING BOUNDARY CONDITIONS(1) << == '
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,1,1, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,4,1, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,1,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,2,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,3,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,1,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,2,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,3,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,4,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,5,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,6,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    IF(FORCE_BC.EQV..TRUE.) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldDelUDelNVariableType,1,1,3,1, &
        & CMISSBoundaryConditionFixed, &
        & -0.2_CMISSDP,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldDelUDelNVariableType,1,1,6,1, &
        & CMISSBoundaryConditionFixed, &
        & -0.2_CMISSDP,Err)
    ELSE
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,3,1, &
        & CMISSBoundaryConditionFixed, &
        & 2.1_CMISSDP,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField1,CMISSFieldUVariableType,1,1,6,1, &
        & CMISSBoundaryConditionFixed, &
        & 2.1_CMISSDP,Err)
    ENDIF
    !Start the creation of the solver equations boundary conditions for the second region
    PRINT *, ' == >> CREATING BOUNDARY CONDITIONS(2) << == '
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,1,1, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,4,1, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,1,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,2,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,3,2, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,1,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,2,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,3,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,4,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,5,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,6,3, &
      & CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)

    IF(FORCE_BC.EQV..TRUE.) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldDelUDelNVariableType,1,1,3,1, &
        & CMISSBoundaryConditionFixed, &
        & -0.2_CMISSDP,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldDelUDelNVariableType,1,1,6,1, &
        & CMISSBoundaryConditionFixed, &
        & -0.2_CMISSDP,Err)
    ELSE
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,3,1, &
        & CMISSBoundaryConditionFixed, &
        & 2.1_CMISSDP,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField2,CMISSFieldUVariableType,1,1,6,1, &
        & CMISSBoundaryConditionFixed, &
        & 2.1_CMISSDP,Err)
    ENDIF
  ENDIF

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(NonLinearSolverEquations,Err)

  !Solve the problem
  PRINT *, ' == >> SOLVING PROBLEM << == '
  CALL CMISSProblemSolve(CoupledProblem,Err)

  !Export the fields
  PRINT *, ' == >> EXPORTING FIELDS << == '
  CALL CMISSFieldsTypeInitialise(Fields1,Err)
  CALL CMISSFieldsTypeCreate(Region1,Fields1,Err)
  CALL CMISSFieldIONodesExport(Fields1,"3DCoupled-FiniteElasticityMembrane-FiniteElasticityMembrane_1","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields1,"3DCoupled-FiniteElasticityMembrane-FiniteElasticityMembrane_1","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields1,Err)
  CALL CMISSFieldsTypeInitialise(Fields2,Err)
  CALL CMISSFieldsTypeCreate(Region2,Fields2,Err)
  CALL CMISSFieldIONodesExport(Fields2,"3DCoupled-FiniteElasticityMembrane-FiniteElasticityMembrane_2","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields2,"3DCoupled-FiniteElasticityMembrane-FiniteElasticityMembrane_2","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields2,Err)
  IF(UNCOUPLED.EQV..FALSE.) THEN
    CALL CMISSFieldsTypeInitialise(Fields3,Err)
    CALL CMISSFieldsTypeCreate(Interface,Fields3,Err)
    CALL CMISSFieldIONodesExport(Fields3,"3DCoupled-FiniteElasticityMembrane-FiniteElasticityMembrane_Interface","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields3,"3DCoupled-FiniteElasticityMembrane-FiniteElasticityMembrane_Interface","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields3,Err)
  ENDIF
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM COUPLEDFINITEELASTICITYFINITEELASTICITY
