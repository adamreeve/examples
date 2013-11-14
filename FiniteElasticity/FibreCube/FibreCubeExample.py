#!/usr/bin/env python
#> \file
#> \author Adam Reeve
#> \brief This is an example script to solve uniaxial extension of a cube with fibres at 60 degrees from the x axis
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s):
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> \example FiniteElasticity/FibreCubeExample.py
## Example script to solve uniaxial extension of a cube with fibres at 60 degrees from the x axis
#<

#> Main script
""" Uniaxial extension of a cube with an anisotropic constitutive
    relation and fibre angle of 60 degrees.
    Matches example 526 from old CMISS
"""
from math import pi
from opencmiss import iron


# Problem parameters
height = 1.0
width = 1.0
length = 1.0

strain = 0.1
numIncrements = 1

# Constant fibre angles in radians
fibreAngles = [
        60.0 * pi / 180.0,
        0.0,
        0.0]

# Constitutive setup
constitutiveRelation = iron.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_EXPONENTIAL
# Note that W = (c1 / 2) e^Q, old CMISS uses W = c1 e^Q
constitutiveParameters = [2.0, 5.0, 10.0, 0.0, 5.0]

# Elements and interpolation
numberGlobalElements = [1, 1, 1]
interpolations = [
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
nodalPressure = False
# Number of xi directions
numberOfXi = 3
# Number of Gauss points used for integration, depends on
# maximum interpolation degree
numberOfGaussXi = 2

# User numbers for OpenCMISS objects
coordinateSystemUserNumber = 1
regionUserNumber = 1
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
(geometricFieldUserNumber,
    fibreFieldUserNumber,
    materialFieldUserNumber,
    dependentFieldUserNumber,
    equationsSetFieldUserNumber,
    deformedFieldUserNumber) = range(1, 7)
equationsSetUserNumber = 1
problemUserNumber = 1

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.LabelSet("Region")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Define bases
bases = []
for basisUserNumber, interpolation in enumerate(interpolations, 1):
    basis = iron.Basis()
    basis.CreateStart(basisUserNumber)
    basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    basis.numberOfXi = numberOfXi
    basis.interpolationXi = [interpolation] * numberOfXi
    basis.quadratureNumberOfGaussXi = [numberOfGaussXi] * numberOfXi
    basis.CreateFinish()
    bases.append(basis)

# Start the creation of a generated mesh in the region
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = bases
generatedMesh.extent = [width, height, length]
generatedMesh.numberOfElements = numberGlobalElements
mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.NumberOfVariablesSet(1)
fibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 3)
for component in range(1, 4):
    fibreField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
            iron.FieldInterpolationTypes.CONSTANT)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")
fibreField.CreateFinish()

# Set fibre angles
for component, angle in enumerate(fibreAngles, 1):
    fibreField.ComponentValuesInitialise(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
        component, angle)

# Create the equations set and equations set field
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
    iron.EquationsSetTypes.FINITE_ELASTICITY,
    constitutiveRelation]
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create the material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

# Set Mooney-Rivlin constants c10 and c01 respectively.
for component, value in enumerate(constitutiveParameters, 1):
    materialField.ComponentValuesInitialise(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
        component, value)

# Create the dependent field for storing the solution
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
if nodalPressure:
    # Set pressure to last mesh component
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 4,
            iron.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.DELUDELN, 4,
            iron.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentMeshComponentSet(
            iron.FieldVariableTypes.U, 4, len(bases))
    dependentField.ComponentMeshComponentSet(
            iron.FieldVariableTypes.DELUDELN, 4, len(bases))
else:
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 4,
            iron.FieldInterpolationTypes.ELEMENT_BASED)
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.DELUDELN, 4,
            iron.FieldInterpolationTypes.ELEMENT_BASED)
equationsSet.DependentCreateFinish()

# Initialise dependent field from undeformed geometry
for component in [1, 2, 3]:
    geometricField.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component,
        dependentField, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component)
# Set initial hydrostatic pressure
dependentField.ComponentValuesInitialise(
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4, -14.0)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.ProblemSubtypes.NONE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.MaximumIterationsSet(numIncrements)
problem.ControlLoopCreateFinish()

# Create problem solver
nonLinearSolver = iron.Solver()
linearSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
nonLinearSolver.NewtonJacobianCalculationTypeSet(
        iron.JacobianCalculationTypes.EQUATIONS)
nonLinearSolver.NewtonRelativeToleranceSet(1.0e-10)
nonLinearSolver.NewtonAbsoluteToleranceSet(1.0e-10)
nonLinearSolver.NewtonSolutionToleranceSet(1.0e-10)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = iron.LinearSolverTypes.DIRECT
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Get nodes in faces
nodes = iron.Nodes()
region.NodesGet(nodes)
domainNodes = set()
leftNodes = set()
rightNodes = set()
frontNodes = set()
backNodes = set()
topNodes = set()
baseNodes = set()

for node in range(1, nodes.numberOfNodes + 1):
    nodeDomain = decomposition.NodeDomainGet(node, 1)
    if nodeDomain != computationalNodeNumber:
        continue
    domainNodes.add(node)

    # get x, y and z positions at node
    position = [
        geometricField.ParameterSetGetNodeDP(
            iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES,
            1, 1, node, component)
        for component in [1, 2, 3]]
    tol = 1.0e-10
    if abs(position[0]) < tol:
        leftNodes.add(node)
    elif abs(position[0] - width) < tol:
        rightNodes.add(node)
    if abs(position[1]) < tol:
        frontNodes.add(node)
    elif abs(position[1] - length) < tol:
        backNodes.add(node)
    if abs(position[2]) < tol:
        baseNodes.add(node)
    elif abs(position[2] - height) < tol:
        topNodes.add(node)

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

# Set x=0 nodes to no x displacment in x
for node in leftNodes.intersection(domainNodes):
    derivative, version = 1, 1
    component = 1
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U,
            version, derivative, node, component,
            iron.BoundaryConditionsTypes.FIXED, 0.0)

# Set x=width nodes to 10% x displacement
for node in rightNodes.intersection(domainNodes):
    derivative, version = 1, 1
    component = 1
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U,
            version, derivative, node, component,
            iron.BoundaryConditionsTypes.FIXED_INCREMENTED,
            (1.0 + strain) * width)

# Set y=0 nodes to no y displacement
for node in frontNodes.intersection(domainNodes):
    derivative, version = 1, 1
    component = 2
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U,
            version, derivative, node, component,
            iron.BoundaryConditionsTypes.FIXED, 0.0)

# Set z=0 nodes to no z displacement
for node in baseNodes.intersection(domainNodes):
    derivative, version = 1, 1
    component = 3
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U,
            version, derivative, node, component,
            iron.BoundaryConditionsTypes.FIXED, 0.0)

solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Create deformed geometry field, as cmgui doesn't like displaying
# deformed fibres from the dependent field
deformedField = iron.Field()
deformedField.CreateStart(deformedFieldUserNumber, region)
deformedField.MeshDecompositionSet(decomposition)
deformedField.TypeSet(iron.FieldTypes.GEOMETRIC)
deformedField.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
for component in [1, 2, 3]:
    deformedField.ComponentMeshComponentSet(
            iron.FieldVariableTypes.U, component, 1)
deformedField.CreateFinish()
# Copy deformed geometry
for component in [1, 2, 3]:
    dependentField.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component,
        deformedField, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component)

# Export results to exnode/exelem files
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("FibreCube", "FORTRAN")
fields.ElementsExport("FibreCube", "FORTRAN")
fields.Finalise()
