!> \file
!> \author Tim Wu
!> \brief This is an example program to solve 3D data points projecting onto 2D cartesian elements using openCMISS calls.
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
!> Contributor(s): Code based on the examples by Kumar Mithraratne and Prasad Babarenda Gamage
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

!> Main program
PROGRAM DataProjection1DRectangularCartesian

  USE MPI
  USE OPENCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Program parameters
  INTEGER(CMFEIntg),PARAMETER :: BasisUserNumber=1  
  INTEGER(CMFEIntg),PARAMETER :: CoordinateSystemDimension=3
  INTEGER(CMFEIntg),PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMFEIntg),PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMFEIntg),PARAMETER :: FieldUserNumber=1  
  INTEGER(CMFEIntg),PARAMETER :: MeshUserNumber=1
  INTEGER(CMFEIntg),PARAMETER :: RegionUserNumber=1

  REAL(CMFEDP), PARAMETER :: CoordinateSystemOrigin(3)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/)  
  !Program types

  !Program variables   
  INTEGER(CMFEIntg) :: MeshComponentNumber=1
  INTEGER(CMFEIntg) :: NumberOfDataPoints
  INTEGER(CMFEIntg) :: MeshDimensions=2
  INTEGER(CMFEIntg) :: MeshNumberOfElements
  INTEGER(CMFEIntg) :: MeshNumberOfComponents=1
  INTEGER(CMFEIntg) :: NumberOfDomains=1 !NumberOfDomains=2 for parallel processing, need to set up MPI
  INTEGER(CMFEIntg) :: NumberOfNodes
  INTEGER(CMFEIntg) :: NumberOfXi=2
  INTEGER(CMFEIntg) :: BasisInterpolation(2)=(/CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION/)
  INTEGER(CMFEIntg) :: WorldCoordinateSystemUserNumber
  INTEGER(CMFEIntg) :: WorldRegionUserNumber
  
  INTEGER(CMFEIntg) :: FieldNumberOfVariables=1
  INTEGER(CMFEIntg) :: FieldNumberOfComponents=3 

  INTEGER(CMFEIntg) :: data_point_idx,elem_idx,ver_idx,der_idx,node_idx,comp_idx
    
  REAL(CMFEDP), DIMENSION(5,3) :: DataPointValues !(number_of_data_points,dimension)
  REAL(CMFEDP), DIMENSION(5) :: DataPointProjectionDistance !(number_of_data_points)
  INTEGER(CMFEIntg), DIMENSION(5) :: DataPointProjectionElementNumber !(number_of_data_points)
  INTEGER(CMFEIntg), DIMENSION(5) :: DataPointProjectionExitTag !(number_of_data_points)
  REAL(CMFEDP), DIMENSION(5,2) :: DataPointProjectionXi !(number_of_data_points,MeshDimensions)  
  INTEGER(CMFEIntg), DIMENSION(4,4) :: ElementUserNodes  
  REAL(CMFEDP), DIMENSION(4,9,3) :: FieldValues
        
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS and MPI variables
  INTEGER(CMFEIntg) :: Err
  INTEGER(CMFEIntg) :: NUMBER_GLOBAL_X_ELEMENTS=1 !<number of elements on x axis
  INTEGER(CMFEIntg) :: NUMBER_GLOBAL_Y_ELEMENTS=1 !<number of elements on y axis
  INTEGER(CMFEIntg) :: NUMBER_GLOBAL_Z_ELEMENTS=1 !<number of elements on z axis  
  INTEGER(CMFEIntg) :: NUMBER_OF_DOMAINS=1      
  INTEGER(CMFEIntg) :: MPI_IERROR  
  
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
    
  !Define data points
  DataPointValues(1,:)=(/5.0_CMFEDP,3.0_CMFEDP,3.0_CMFEDP/)
  DataPointValues(2,:)=(/6.0_CMFEDP,12.8_CMFEDP,-10.0_CMFEDP/)  
  DataPointValues(3,:)=(/16.0_CMFEDP,6.9_CMFEDP,-5.0_CMFEDP/)  
  DataPointValues(4,:)=(/21.0_CMFEDP,11.0_CMFEDP,20.0_CMFEDP/)
  DataPointValues(5,:)=(/24.0_CMFEDP,21.5_CMFEDP,10.0_CMFEDP/)
  NumberOfDataPoints=SIZE(DataPointValues,1)
  !Define element connectivities
  ElementUserNodes(1,:)=(/1,2,4,5/)
  ElementUserNodes(2,:)=(/2,3,5,6/)
  ElementUserNodes(3,:)=(/4,5,7,8/)
  ElementUserNodes(4,:)=(/5,6,8,9/)
  MeshNumberOfElements=SIZE(ElementUserNodes,1)     
  !Define nodal fields
  FieldValues(1,1,:)=(/0.0_CMFEDP,0.0_CMFEDP,10.0_CMFEDP/) !no der, node 1
  FieldValues(2,1,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 1
  FieldValues(3,1,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 1  
  FieldValues(4,1,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 1    
  
  FieldValues(1,2,:)=(/10.0_CMFEDP,0.0_CMFEDP,-5.0_CMFEDP/) !no der, node 2
  FieldValues(2,2,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 2
  FieldValues(3,2,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 2  
  FieldValues(4,2,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 2    
  
  FieldValues(1,3,:)=(/20.0_CMFEDP,0.0_CMFEDP,10.0_CMFEDP/) !no der, node 3
  FieldValues(2,3,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 3
  FieldValues(3,3,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 3  
  FieldValues(4,3,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 3   
  
  FieldValues(1,4,:)=(/0.0_CMFEDP,10.0_CMFEDP,-5.0_CMFEDP/) !no der, node 4
  FieldValues(2,4,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 4
  FieldValues(3,4,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 4  
  FieldValues(4,4,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 4 
  
  FieldValues(1,5,:)=(/10.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !no der, node 5
  FieldValues(2,5,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 5
  FieldValues(3,5,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 5  
  FieldValues(4,5,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 5   
  
  FieldValues(1,6,:)=(/20.0_CMFEDP,10.0_CMFEDP,-5.0_CMFEDP/) !no der, node 6
  FieldValues(2,6,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 6  
  FieldValues(3,6,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 6  
  FieldValues(4,6,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 6   
  
  FieldValues(1,7,:)=(/0.0_CMFEDP,20.0_CMFEDP,10.0_CMFEDP/) !no der, node 7
  FieldValues(2,7,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 7
  FieldValues(3,7,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 7  
  FieldValues(4,7,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 7
    
  FieldValues(1,8,:)=(/10.0_CMFEDP,20.0_CMFEDP,-5.0_CMFEDP/) !no der, node 8
  FieldValues(2,8,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 8
  FieldValues(3,8,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 8  
  FieldValues(4,8,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 8   
    
  FieldValues(1,9,:)=(/20.0_CMFEDP,20.0_CMFEDP,10.0_CMFEDP/) !no der, node 9
  FieldValues(2,9,:)=(/10.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 der, node 9
  FieldValues(3,9,:)=(/0.0_CMFEDP,10.0_CMFEDP,0.0_CMFEDP/) !s2 der, node 9  
  FieldValues(4,9,:)=(/0.0_CMFEDP,0.0_CMFEDP,0.0_CMFEDP/) !s1 s2 der, node 9
  NumberOfNodes=SIZE(FieldValues,2)
  
  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,Err)
  !Broadcast the number of Elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !=========================================================================================================================
  !Create RC coordinate system
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,Err)
  CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystemUserNumber,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemUserNumber,CoordinateSystemDimension,Err)
  CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystemUserNumber,CoordinateSystemOrigin,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemUserNumber,Err) 

  !=========================================================================================================================
  !Create Region and set CS to newly created 3D RC CS
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegionUserNumber,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionUserNumber,CoordinateSystemUserNumber,Err)
  CALL cmfe_Region_CreateFinish(RegionUserNumber,Err)
    
  !=========================================================================================================================
  !Create Data Points and set the values
  CALL cmfe_DataPoints_CreateStart(RegionUserNumber,SIZE(DataPointValues,1),Err)
  DO data_point_idx=1,NumberOfDataPoints
    CALL cmfe_DataPoints_ValuesSet(RegionUserNumber,data_point_idx,DataPointValues(data_point_idx,:),Err)     
  ENDDO
  CALL cmfe_DataPoints_CreateFinish(RegionUserNumber,Err)  
  !=========================================================================================================================
  !Define basis function - 1D cubic hermite
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Err)
  CALL cmfe_Basis_TypeSet(BasisUserNumber,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(BasisUserNumber,NumberOfXi,Err)
  CALL cmfe_Basis_InterpolationXiSet(BasisUserNumber,BasisInterpolation,Err)
  CALL cmfe_Basis_CreateFinish(BasisUserNumber,Err)  
  !=========================================================================================================================
  !Create a mesh
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,RegionUserNumber,MeshDimensions,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfComponents,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfElements,Err)
  !define nodes for the mesh
  CALL cmfe_Nodes_CreateStart(RegionUserNumber,NumberOfNodes,Err)
  CALL cmfe_Nodes_CreateFinish(RegionUserNumber,Err)  
  !define elements for the mesh
  CALL cmfe_MeshElements_CreateStart(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber,Err)
  Do elem_idx=1,MeshNumberOfElements
    CALL cmfe_MeshElements_NodesSet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,elem_idx,ElementUserNodes(elem_idx,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(RegionUserNumber,MeshUserNumber,MeshComponentNumber,Err)
  CALL cmfe_Mesh_CreateFinish(RegionUserNumber,MeshUserNumber,Err)
  !=========================================================================================================================
  !Create a mesh decomposition 
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,RegionUserNumber,MeshUserNumber,Err)
  CALL cmfe_Decomposition_TypeSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  
  !=========================================================================================================================
  !Create a field to put the geometry
  CALL cmfe_Field_CreateStart(FieldUserNumber,RegionUserNumber,Err)
  CALL cmfe_Field_MeshDecompositionSet(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  CALL cmfe_Field_TypeSet(RegionUserNumber,FieldUserNumber,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(RegionUserNumber,FieldUserNumber,FieldNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(RegionUserNumber,FieldUserNumber,CMFE_FIELD_U_VARIABLE_TYPE,FieldNumberOfComponents,Err)
  !DO xi=1,NumberOfXi
  CALL cmfe_Field_ComponentMeshComponentSet(RegionUserNumber,FieldUserNumber,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  !ENDDO !xi    
  CALL cmfe_Field_CreateFinish(RegionUserNumber,FieldUserNumber,Err)
  !node 1
  DO der_idx=1,SIZE(FieldValues,1)
    DO node_idx=1,SIZE(FieldValues,2)
      DO comp_idx=1,SIZE(FieldValues,3)
        CALL cmfe_Field_ParameterSetUpdateNode(RegionUserNumber,FieldUserNumber,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE, &
          & ver_idx,der_idx,node_idx,comp_idx,FieldValues(der_idx,node_idx,comp_idx),Err)
      ENDDO
    ENDDO
  ENDDO
  
  !=========================================================================================================================
  !Create a data projection
  
  CALL cmfe_DataProjection_CreateStart(RegionUserNumber,FieldUserNumber,RegionUserNumber,Err)
  CALL cmfe_DataProjection_CreateFinish(RegionUserNumber,Err)
  
  !=========================================================================================================================
  !Start data projection
  CALL cmfe_DataProjection_Evaluate(RegionUserNumber,Err)
  
  !=========================================================================================================================
  !Retrieve projection results
  DO data_point_idx=1,NumberOfDataPoints
    CALL cmfe_DataPoints_ProjectionDistanceGet(RegionUserNumber,data_point_idx,DataPointProjectionDistance(data_point_idx),Err)
    CALL cmfe_DataPoints_ProjectionElementNumberGet(RegionUserNumber,data_point_idx,DataPointProjectionElementNumber( &
      & data_point_idx),Err)
    CALL cmfe_DataPoints_ProjectionExitTagGet(RegionUserNumber,data_point_idx,DataPointProjectionExitTag(data_point_idx),Err)
    CALL cmfe_DataPoints_ProjectionXiGet(RegionUserNumber,data_point_idx,DataPointProjectionXi(data_point_idx,:),Err)
  ENDDO
  
  !=========================================================================================================================
  !Destroy used types
  CALL cmfe_DataProjection_Destroy(RegionUserNumber,Err)
  CALL cmfe_DataPoints_Destroy(RegionUserNumber,Err)
    
  CALL cmfe_Region_Destroy(RegionUserNumber,Err)
  CALL cmfe_CoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)  
  
  !=========================================================================================================================
  !Finishing program
  CALL cmfe_Finalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  STOP  
  
END PROGRAM DataProjection1DRectangularCartesian
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
