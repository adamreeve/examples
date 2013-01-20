!> \file
!> \author Adam Reeve
!> \brief This is an example program to solve a finite elasticity problem with
!> stress boundary conditions using openCMISS calls.
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

!> \example FiniteElasticity/UniAxialStress/src/UniAxialStressExample.f90
!! Example program to solve a finite elasticity problem with stress boundary conditions using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialStress/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialStress/build-gnu'>Linux GNU Build</a>
!<

PROGRAM UniAxialStressExample

  USE OPENCMISS
  USE MPI

  IMPLICIT NONE

  !problem setup parameters
  REAL(CMISSDP), PARAMETER :: stress=10.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: height=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: width=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: length=1.0_CMISSDP
  INTEGER(CMISSIntg), PARAMETER :: interpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: pressureInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  LOGICAL, PARAMETER :: usePressureBasis=.FALSE.
  INTEGER(CMISSIntg), PARAMETER :: numberOfGaussXi=2
  INTEGER(CMISSIntg) :: numberGlobalXElements=8
  INTEGER(CMISSIntg) :: numberGlobalYElements=8
  INTEGER(CMISSIntg) :: numberGlobalZElements=8
  INTEGER(CMISSIntg) :: numberLoadSteps=2

  !CMISS object user numbers
  INTEGER(CMISSIntg), PARAMETER :: coordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: regionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: basisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: pressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: generatedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: meshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: decompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: fieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: fieldMaterialUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: fieldDependentUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: fieldEquationsSetUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: equationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: problemUserNumber=1

  !Program variables
  INTEGER(CMISSIntg) :: equationsSetIndex
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,numberOfDomains,computationalNodeNumber
  INTEGER(CMISSIntg) :: nodeNumber,nodeDomain,nodeIdx,numberOfNodes
  INTEGER(CMISSIntg) :: componentIdx
  REAL(CMISSDP) :: nodePos(3)
  REAL(CMISSDP) :: eps=1.0e-10
  REAL(CMISSDP) :: force, cnForce, totalForce, resultStress
  INTEGER(CMISSIntg) :: mpiError

  !CMISS variables
  TYPE(CMISSBasisType) :: basis,pressureBasis
  TYPE(CMISSBoundaryConditionsType) :: boundaryConditions
  TYPE(CMISSControlLoopType) :: controlLoop
  TYPE(CMISSCoordinateSystemType) :: coordinateSystem,worldCoordinateSystem
  TYPE(CMISSDecompositionType) :: decomposition
  TYPE(CMISSEquationsType) :: equations
  TYPE(CMISSEquationsSetType) :: equationsSet
  TYPE(CMISSFieldType) :: geometricField,materialField,dependentField,equationsSetField
  TYPE(CMISSFieldsType) :: fields
  TYPE(CMISSGeneratedMeshType) :: generatedMesh
  TYPE(CMISSMeshType) :: mesh
  TYPE(CMISSNodesType) :: nodes
  TYPE(CMISSProblemType) :: problem
  TYPE(CMISSRegionType) :: region,worldRegion
  TYPE(CMISSSolverType) :: solver,linearSolver
  TYPE(CMISSSolverEquationsType) :: solverEquations

  INTEGER(CMISSIntg) :: err

  !Initialise OpenCMISS
  CALL CMISSInitialise(worldCoordinateSystem,worldRegion,err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(numberOfComputationalNodes,err)
  CALL CMISSComputationalNodeNumberGet(computationalNodeNumber,err)

  numberOfDomains=numberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystem_Initialise(coordinateSystem,err)
  CALL CMISSCoordinateSystem_CreateStart(coordinateSystemUserNumber,coordinateSystem,err)
  CALL CMISSCoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(region,err)
  CALL CMISSRegion_CreateStart(regionUserNumber,worldRegion,region,err)
  CALL CMISSRegion_LabelSet(region,"example_region",err)
  CALL CMISSRegion_CoordinateSystemSet(region,coordinateSystem,err)
  CALL CMISSRegion_CreateFinish(region,err)

  !Define geometric basis
  CALL CMISSBasis_Initialise(basis,err)
  CALL CMISSBasis_CreateStart(basisUserNumber,basis,err)
  SELECT CASE(InterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(basis,CMISS_BASIS_SIMPLEX_TYPE,err)
  END SELECT
  IF(numberGlobalZElements==0) THEN
    CALL CMISSBasis_NumberOfXiSet(basis,2,err)
    CALL CMISSBasis_InterpolationXiSet(basis,[InterpolationType,InterpolationType],err)
    IF(NumberOfGaussXi>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(basis,[NumberOfGaussXi,NumberOfGaussXi],err)
    END IF
  ELSE
    CALL CMISSBasis_NumberOfXiSet(basis,3,err)
    CALL CMISSBasis_InterpolationXiSet(basis,[InterpolationType,InterpolationType,InterpolationType],err)
    IF(NumberOfGaussXi>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(basis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],err)
    END IF
  END IF
  CALL CMISSBasis_CreateFinish(basis,err)

  !Define pressure basis
  IF(UsePressureBasis) THEN
    CALL CMISSBasis_Initialise(pressureBasis,err)
    CALL CMISSBasis_CreateStart(pressureBasisUserNumber,pressureBasis,err)
    SELECT CASE(pressureInterpolationType)
    CASE(1,2,3,4)
      CALL CMISSBasis_TypeSet(pressureBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
    CASE(7,8,9)
      CALL CMISSBasis_TypeSet(pressureBasis,CMISS_BASIS_SIMPLEX_TYPE,err)
    END SELECT
    IF(numberGlobalZElements==0) THEN
      CALL CMISSBasis_NumberOfXiSet(pressureBasis,2,err)
      CALL CMISSBasis_InterpolationXiSet(pressureBasis,[pressureInterpolationType,pressureInterpolationType],err)
      IF(NumberOfGaussXi>0) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(pressureBasis,[NumberOfGaussXi,NumberOfGaussXi],err)
      END IF
    ELSE
      CALL CMISSBasis_NumberOfXiSet(pressureBasis,3,err)
      CALL CMISSBasis_InterpolationXiSet(pressureBasis, &
        & [pressureInterpolationType,pressureInterpolationType,pressureInterpolationType],err)
      IF(NumberOfGaussXi>0) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(pressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],err)
      END IF
    END IF
    CALL CMISSBasis_CreateFinish(pressureBasis,err)
  END IF

  !Create a regular generated mesh with the specified dimensions and bases
  CALL CMISSGeneratedMesh_Initialise(generatedMesh,err)
  CALL CMISSGeneratedMesh_CreateStart(generatedMeshUserNumber,region,generatedMesh,err)
  CALL CMISSGeneratedMesh_TypeSet(generatedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  IF(usePressureBasis) THEN
    CALL CMISSGeneratedMesh_BasisSet(generatedMesh,[basis,pressureBasis],err)
  ELSE
    CALL CMISSGeneratedMesh_BasisSet(generatedMesh,[basis],err)
  END IF
  IF(numberGlobalXElements==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(generatedMesh,[width,height],err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(generatedMesh,[numberGlobalXElements,numberGlobalYElements],err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(generatedMesh,[width,height,length],err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(generatedMesh,[numberGlobalXElements,numberGlobalYElements, &
      & numberGlobalZElements],err)
  END IF
  CALL CMISSMesh_Initialise(mesh,err)
  CALL CMISSGeneratedMesh_CreateFinish(generatedMesh,meshUserNumber,mesh,err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(decomposition,err)
  CALL CMISSDecomposition_CreateStart(decompositionUserNumber,mesh,decomposition,err)
  CALL CMISSDecomposition_TypeSet(decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,err)
  CALL CMISSDecomposition_NumberOfDomainsSet(decomposition,numberOfDomains,err)
  CALL CMISSDecomposition_CalculateFacesSet(decomposition,.TRUE.,err)
  CALL CMISSDecomposition_CreateFinish(decomposition,err)

  !Create a geometric field to store the undeformed geometry
  CALL CMISSField_Initialise(geometricField,err)
  CALL CMISSField_CreateStart(fieldGeometryUserNumber,region,geometricField,err)
  CALL CMISSField_MeshDecompositionSet(geometricField,decomposition,err)
  CALL CMISSField_VariableLabelSet(geometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL CMISSField_CreateFinish(geometricField,err)

  !Update the geometric field parameters from the genreated mesh
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create the equations set
  CALL CMISSField_Initialise(equationsSetField,err)
  CALL CMISSEquationsSet_CreateStart(equationSetUserNumber,region,geometricField, &
    & CMISS_EQUATIONS_SET_ELASTICITY_CLASS,CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
    & fieldEquationsSetUserNumber,equationsSetField,equationsSet,err)
  CALL CMISSEquationsSet_CreateFinish(equationsSet,err)

  !Create the dependent field
  CALL CMISSField_Initialise(dependentField,err)
  CALL CMISSEquationsSet_DependentCreateStart(equationsSet,fieldDependentUserNumber,dependentField,err)
  CALL CMISSField_VariableLabelSet(dependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  IF(usePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL CMISSField_ComponentInterpolationSet(dependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL CMISSField_ComponentInterpolationSet(dependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMISS_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL CMISSField_ComponentMeshComponentSet(dependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,2,err)
    CALL CMISSField_ComponentMeshComponentSet(dependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,2,err)
  END IF
  CALL CMISSEquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the material field
  CALL CMISSField_Initialise(materialField,err)
  CALL CMISSEquationsSet_MaterialsCreateStart(equationsSet,fieldMaterialUserNumber,materialField,err)
  CALL CMISSField_VariableLabelSet(materialField,CMISS_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Set Mooney-Rivlin constants c10 and c01.
  CALL CMISSField_ComponentValuesInitialise(materialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,err)
  CALL CMISSField_ComponentValuesInitialise(materialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,6.0_CMISSDP,err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(equations,err)
  CALL CMISSEquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL CMISSEquations_SparsityTypeSet(equations,CMISS_EQUATIONS_SPARSE_MATRICES,err)
  CALL CMISSEquations_OutputTypeSet(equations,CMISS_EQUATIONS_NO_OUTPUT,err)
  CALL CMISSEquationsSet_EquationsCreateFinish(equationsSet,err)

  !Initialise dependent field from undeformed geometry and set the intial hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(geometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(geometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(geometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,err)
  CALL CMISSField_ComponentValuesInitialise(dependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSDP, &
    & err)

  !Define the problem
  CALL CMISSProblem_Initialise(problem,err)
  CALL CMISSProblem_CreateStart(problemUserNumber,problem,err)
  CALL CMISSProblem_SpecificationSet(problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,err)
  CALL CMISSProblem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(problem,err)
  CALL CMISSControlLoop_Initialise(controlLoop,err)
  CALL CMISSProblem_ControlLoopGet(problem,CMISS_CONTROL_LOOP_NODE,controlLoop,err)
  CALL CMISSControlLoop_MaximumIterationsSet(controlLoop,numberLoadSteps,err)
  CALL CMISSProblem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(solver,err)
  CALL CMISSSolver_Initialise(linearSolver,err)
  CALL CMISSProblem_SolversCreateStart(problem,err)
  CALL CMISSProblem_SolverGet(problem,CMISS_CONTROL_LOOP_NODE,1,solver,err)
  CALL CMISSSolver_OutputTypeSet(solver,CMISS_SOLVER_PROGRESS_OUTPUT,err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  CALL CMISSSolver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(solver,1.0e-12_CMISSDP,err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(solver,1.0e-12_CMISSDP,err)
  CALL CMISSSolver_LinearTypeSet(linearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL CMISSProblem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(solver,err)
  CALL CMISSSolverEquations_Initialise(solverEquations,err)
  CALL CMISSProblem_SolverEquationsCreateStart(problem,err)
  CALL CMISSProblem_SolverGet(problem,CMISS_CONTROL_LOOP_NODE,1,solver,err)
  CALL CMISSSolver_SolverEquationsGet(solver,solverEquations,err)
  CALL CMISSSolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL CMISSProblem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions
  CALL CMISSNodes_Initialise(nodes,err)
  CALL CMISSRegion_NodesGet(region,nodes,err)

  CALL CMISSBoundaryConditions_Initialise(boundaryConditions,err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  CALL CMISSNodes_NumberOfNodesGet(nodes,numberOfNodes,err)
  DO nodeIdx=1,numberOfNodes
    CALL CMISSNodes_UserNumberGet(nodes,nodeIdx,nodeNumber,err)
    CALL CMISSDecomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      DO componentIdx=1,3
        CALL CMISSField_ParameterSetGetNode( &
          & geometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
          & 1,1,nodeNumber,componentIdx,nodePos(componentIdx),err)
      END DO

      IF(ABS(nodePos(1)) < eps) THEN
        ! Left face
        ! Fix x value at 0.0
        CALL CMISSBoundaryConditions_SetNode(boundaryConditions,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,err)
      ELSE IF(ABS(nodePos(1) - width) < eps) THEN
        ! Right face
        ! Apply stress using incremented point Neumann conditions
        ! The stress value will be integrated to calculate nodal forces
        CALL CMISSBoundaryConditions_SetNode(boundaryConditions,dependentField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & CMISS_BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED,stress,err)
      ELSE IF((ABS(nodePos(2) - height) < eps) .OR. (ABS(nodePos(2)) < eps) .OR. &
          & (ABS(nodePos(3) - length) < eps) .OR. (ABS(nodePos(3)) < eps)) THEN
        ! Apply integrated only free conditions elsewhere to avoid integrating the
        ! Neumann point conditions across the faces except for the right side
        CALL CMISSBoundaryConditions_SetNode(boundaryConditions,dependentField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & CMISS_BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY,0.0_CMISSDP,err)
      END IF

      IF(ABS(nodePos(2)) < eps) THEN
        ! Front face
        ! Fix y value at 0.0
        CALL CMISSBoundaryConditions_SetNode(boundaryConditions,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
          & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,err)
      END IF

      IF(ABS(nodePos(3)) < eps) THEN
        ! Bottom face
        ! Fix z value at 0.0
        CALL CMISSBoundaryConditions_SetNode(boundaryConditions,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
          & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,err)
      END IF
    END IF
  END DO

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL CMISSProblem_Solve(problem,err)

  !Measure total force applied just to check stress application
  !Have to measure at opposite side due to way boundary conditions work
  cnForce = 0.0_CMISSDP
  DO nodeIdx=1,numberOfNodes
    CALL CMISSNodes_UserNumberGet(nodes,nodeIdx,nodeNumber,err)
    CALL CMISSDecomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      DO componentIdx=1,3
        CALL CMISSField_ParameterSetGetNode( &
          & geometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
          & 1,1,nodeNumber,componentIdx,nodePos(componentIdx),err)
      END DO

      IF(ABS(nodePos(1)) < eps) THEN
        CALL CMISSField_ParameterSetGetNode( &
          & dependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
          & 1,1,nodeNumber,1,force,err)
        cnForce = cnForce + force
      END IF
    END IF
  END DO

  !Last node is opposite origin, use it to get deformed dimensions
  CALL CMISSNodes_UserNumberGet(nodes,numberOfNodes,nodeNumber,err)
  CALL CMISSDecomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
  IF(nodeDomain==computationalNodeNumber) THEN
    DO componentIdx=1,3
      CALL CMISSField_ParameterSetGetNode( &
        & dependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
        & 1,1,nodeNumber,componentIdx,nodePos(componentIdx),err)
    END DO
    CALL MPI_BCAST(nodePos,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpiError)
  END IF

  IF(computationalNodeNumber==0) THEN
    CALL MPI_REDUCE(cnForce, totalForce, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpiError)
    resultStress=totalForce / (nodePos(2) * nodePos(3))
    !Stress integration is done over deformed geometry so totalForce will depend
    !on the change in surface area on the right face.
    WRITE(*,*) "Total force:", totalForce
    WRITE(*,*) "Stress:", resultStress
  END IF

  !Output solution
  CALL CMISSFields_Initialise(fields,err)
  CALL CMISSFields_Create(region,fields,err)
  CALL CMISSFields_NodesExport(fields,"uniaxial_stress","FORTRAN",err)
  CALL CMISSFields_ElementsExport(fields,"uniaxial_stress","FORTRAN",err)
  CALL CMISSFields_Finalise(fields,err)

  CALL CMISSFinalise(err)

  IF(ABS(resultStress + stress) < 1.0e-2) THEN
    WRITE(*,'(A)') "Program successfully completed."
  ELSE
    WRITE(*,'(A)') "Measured stress does not match applied stress."
    STOP 1
  END IF

  STOP

END PROGRAM UniAxialStressExample
